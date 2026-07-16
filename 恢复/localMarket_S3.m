function res = localMarket_S3(c,b,a,D,init_rho,pg_max,pg_min,init_pi, ...
    init_pes,init_pg,init_pb,init_ps,init_pesc,errTol,inner_cv_tol_kW, ...
    pi_0,pi_max,pi_min,record_trace,beta_qre,~,~,qre_z_cap,qre_max_backoffs, ...
    inner_momentum,max_inner_iter,h0,min_inner_iter,inner_avg_start_iter, ...
    qre_certificate_enabled,agg_certificate_enabled)

A = max_inner_iter + 2;
lambda_step = nan(A,1);
N = length(b);

validateattributes(init_rho, {'numeric'}, {'scalar','real','finite','positive'});
validateattributes(inner_cv_tol_kW, {'numeric'}, {'scalar','real','finite','nonnegative'});
validateattributes(h0, {'numeric'}, {'scalar','real','finite','nonnegative','integer'});
validateattributes(min_inner_iter, {'numeric'}, {'scalar','real','finite','integer','nonnegative'});
validateattributes(inner_avg_start_iter, {'numeric'}, {'scalar','real','finite','positive','integer'});
if nargin < 31
    qre_certificate_enabled = false;
    agg_certificate_enabled = false;
end
if qre_certificate_enabled || agg_certificate_enabled
    error('localMarket_S3 is currently the unchecked_qre profile; enable certificates only after restoring the certified implementation.');
end
rho = init_rho;
alpha = rho * a;

%% 加异步
delay = randi([0,h0],N,1);
count = zeros(N,1);

h = 1;

val_lambda = zeros(A,1);
val_lambda(1) = init_pi;
val_lambda_au = zeros(A,1);
val_lambda_au(1) = init_pi;
val_lambda_used = nan(A,1);
val_lambda_used(1) = init_pi;

val_pes = zeros(N,A);
val_pes(:,1) = init_pes;

val_pg = zeros(N,A);
val_pg(:,1) = init_pg;

val_pb = zeros(N,A);
val_pb(:,1) = init_pb;

val_ps = zeros(N,A);
val_ps(:,1) = init_ps;

val_pesc = zeros(A,1);
val_pesc(1) = init_pesc;

val_violation = zeros(A,1);
val_violation(1) = abs(val_pesc(1) - sum(init_pes)); % 约束违反单位：kW

val_obj = nan(A,1);
val_obj(1) = euLagrangianObj(val_pg(:,1), val_pb(:,1), ...
    val_ps(:,1), val_pes(:,1), init_pi);

fast_pg  = init_pg;
fast_pes = init_pes;
fast_pb  = init_pb;
fast_ps  = init_ps;

sum_pes = zeros(N,1);
sum_pg  = zeros(N,1);
sum_pb  = zeros(N,1);
sum_ps  = zeros(N,1);
sum_pesc = 0;
avg_cnt = 0;

res.success = 1;
res.stop_reason = '';
qre_calls = 0;
qre_perturbed_accepts = 0;
qre_backoff_steps = 0;
qre_fallbacks = 0;
qre_boundary_hits = 0;
qre_max_gap = NaN;

while true

    for j = 1:N
        if count(j) < delay(j)
            count(j) = count(j) + 1;
            fast_pes(j) = val_pes(j,h);
            fast_pg(j)  = val_pg(j,h);
            fast_pb(j)  = val_pb(j,h);
            fast_ps(j)  = val_ps(j,h);
        else
            [fast_pg(j), fast_pes(j), fast_pb(j), fast_ps(j)] = ...
                solveBranch(j, val_lambda(h));
            count(j) = 0;
            delay(j) = randi([0,h0]);
        end

    end

    val_pes(:,h+1) = fast_pes;
    val_pg(:,h+1)  = fast_pg;
    val_pb(:,h+1)  = fast_pb;
    val_ps(:,h+1)  = fast_ps;

    val_pesc(h+1) = (pi_0 - val_lambda(h)) / a;

    if h >= inner_avg_start_iter
        avg_cnt = avg_cnt + 1;
        sum_pes = sum_pes + fast_pes;
        sum_pg  = sum_pg  + fast_pg;
        sum_pb  = sum_pb  + fast_pb;
        sum_ps  = sum_ps  + fast_ps;
        sum_pesc = sum_pesc + val_pesc(h+1);

        bar_pes = sum_pes / avg_cnt;
        bar_pg  = sum_pg  / avg_cnt;
        bar_pb  = sum_pb  / avg_cnt;
        bar_ps  = sum_ps  / avg_cnt;
        bar_pesc_raw = sum_pesc / avg_cnt; % Agg 后缀平均输出
    else
        % burn-in 不进入正式平均；此时输出当前瞬时状态。
        bar_pes = fast_pes;
        bar_pg  = fast_pg;
        bar_pb  = fast_pb;
        bar_ps  = fast_ps;
        bar_pesc_raw = val_pesc(h+1);
    end
    % 对照恢复版：正式上传为 EU 后缀平均输出之和；raw Agg 仍保留用于真实 CV。
    bar_pesc = sum(bar_pes);
    val_lambda_au(h+1) = val_lambda(h) + ...
        alpha * (val_pesc(h+1) - sum(val_pes(:,h+1)));
    val_lambda(h+1) = val_lambda_au(h+1) + inner_momentum * ...
        (val_lambda_au(h+1) - val_lambda_au(h));

    if record_trace
        val_violation(h+1) = abs(bar_pesc_raw - sum(bar_pes)); % 单位：kW
        % EU 输出由本步更新前的 lambda 求得；RE 必须使用同一固定 lambda。
        val_lambda_used(h+1) = val_lambda(h);
        val_obj(h+1) = euLagrangianObj(fast_pg, fast_pb, fast_ps, ...
            fast_pes, val_lambda_used(h+1));
    end

    lambda_step(h) = abs(val_lambda(h+1) - val_lambda(h));
    inner_cv_bar_now_kW = abs(bar_pesc_raw - sum(bar_pes));

    % 价格稳定与正式遍历输出的真实内层 CV 必须同时满足。
    if h >= min_inner_iter && lambda_step(h) <= errTol && ...
            inner_cv_bar_now_kW <= inner_cv_tol_kW
        res.stop_reason = 'price_step_and_average_cv';
        break;
    end

    if h >= max_inner_iter
        res.success = 0;
        res.stop_reason = 'max_inner_iter_unchecked_qre';
        break;
    end

    h = h + 1;

end

if record_trace
    trace_end = h + 1;
    res.trace_h = (0:h)';
    res.trace_violation = abs(val_violation(1:trace_end));
    res.trace_obj = val_obj(1:trace_end);
    res.trace_lambda = val_lambda(1:trace_end);
    res.trace_lambda_used = val_lambda_used(1:trace_end);
else
    res.trace_h = [];
    res.trace_violation = [];
    res.trace_obj = [];
    res.trace_lambda = [];
    res.trace_lambda_used = [];
end

res.pi = val_lambda(h+1);

res.pes = bar_pes;
res.pg  = bar_pg;
res.pb  = bar_pb;
res.ps  = bar_ps;

res.pesc = bar_pesc;
res.pesc_agg_raw = bar_pesc_raw;
res.pes_sum = sum(bar_pes);
res.inner_cv_kW = abs(bar_pesc_raw - res.pes_sum);
res.inner_avg_cnt = avg_cnt;
res.inner_avg_start_iter = inner_avg_start_iter;
res.xj2 = sum(res.pes .^ 2);
res.iter = h;
res.price_updates = h;
res.qre_calls = qre_calls;
res.qre_perturbed_accepts = qre_perturbed_accepts;
res.qre_backoff_steps = qre_backoff_steps;
res.qre_fallbacks = qre_fallbacks;
res.qre_boundary_hits = qre_boundary_hits;
res.qre_max_gap = qre_max_gap;
res.qre_epsilon = NaN;
res.qre_certificate_enabled = false;
res.agg_obj_ref = NaN;
res.agg_obj_raw = NaN;
res.agg_obj_gap = NaN;
res.agg_epsilon = NaN;
res.agg_certificate_enabled = false;
res.agg_objective_certified = false;
res.agg_a2_approx_certified = false;
res.agg_exact_meta = struct('problem',NaN,'info','disabled_for_unchecked_qre_profile', ...
    'solver','none','objective_source','not_evaluated');

    function [pg, pes, pb, ps] = solveBranch(j, lambda)
        [pg, pes, pb, ps, cert] = boundedQREBranch(c(j), b(j), a, D(j), ...
            pg_min(j), pg_max(j), lambda, pi_max, pi_min, beta_qre(j), ...
            NaN, qre_z_cap, qre_max_backoffs);
        qre_calls = qre_calls + 1;
        qre_perturbed_accepts = qre_perturbed_accepts + cert.perturbed_accepted;
        qre_backoff_steps = qre_backoff_steps + cert.backoff_steps;
        qre_fallbacks = qre_fallbacks + cert.used_fallback;
        qre_boundary_hits = qre_boundary_hits + cert.boundary_hit;
        if isfinite(cert.objective_gap)
            qre_max_gap = max(qre_max_gap, cert.objective_gap);
        end
    end

    function obj = euLagrangianObj(pg_v, pb_v, ps_v, pes_v, lambda_v)
        obj = sum(0.5 .* c(:) .* pg_v(:).^2 ...
            + b(:) .* pg_v(:) ...
            + pi_max .* pb_v(:) ...
            - pi_min .* ps_v(:) ...
            + 0.5 .* a .* pes_v(:).^2 ...
            - lambda_v .* pes_v(:));

    end

end
