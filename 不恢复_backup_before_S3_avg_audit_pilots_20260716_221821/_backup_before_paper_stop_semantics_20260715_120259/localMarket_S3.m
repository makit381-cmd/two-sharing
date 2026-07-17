function res = localMarket_S3(c,b,a,D,init_rho,pg_max,pg_min,init_pi, ...
    init_pes,init_pg,init_pb,init_ps,init_pesc,errTol,inner_cv_tol_kW, ...
    pi_0,pi_max,pi_min,record_trace,beta_qre,~,~,qre_z_cap,qre_max_backoffs, ...
    inner_momentum,max_inner_iter,h0,min_inner_iter,inner_avg_start_iter, ...
    qre_certificate_enabled,agg_certificate_enabled,progress_print_enabled, ...
    progress_print_inner_every,outer_iter,community_id)

% S3 内层只保留当前 EU 状态与遍历平均；不再保存每一步全 EU 历史矩阵。
% 瞬时 raw 残差只用于 lambda 更新；正式 CV/RE 使用遍历平均输出。

N = length(b);

validateattributes(init_rho, {'numeric'}, {'scalar','real','finite','positive'});
validateattributes(inner_cv_tol_kW, {'numeric'}, {'scalar','real','finite','nonnegative'});
validateattributes(h0, {'numeric'}, {'scalar','real','finite','nonnegative','integer'});
validateattributes(min_inner_iter, {'numeric'}, {'scalar','real','finite','integer','nonnegative'});
validateattributes(inner_avg_start_iter, {'numeric'}, {'scalar','real','finite','positive','integer'});
if nargin < 35
    progress_print_enabled = false;
    progress_print_inner_every = inf;
    outer_iter = NaN;
    community_id = NaN;
end
if qre_certificate_enabled || agg_certificate_enabled
    error('localMarket_S3 is currently the unchecked_qre profile; enable certificates only after restoring the certified implementation.');
end

rho = init_rho;
alpha = rho * a;
delay = randi([0,h0],N,1);
count = zeros(N,1);

% 当前状态；异步 held 的 EU 直接保留其上一状态，无需全历史数组。
lambda = init_pi;
lambda_aux_prev = init_pi;
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
bar_pes = init_pes;
bar_pg = init_pg;
bar_pb = init_pb;
bar_ps = init_ps;
bar_pesc_raw = init_pesc;

% 仅选中的社区保存 RE/CV 轨迹；这些是一维数组，按需扩容。
if record_trace
    trace_capacity = min(max_inner_iter + 1, 128);
    trace_h = nan(trace_capacity,1);
    trace_violation = nan(trace_capacity,1);
    trace_obj = nan(trace_capacity,1);
    trace_lambda = nan(trace_capacity,1);
    trace_lambda_used = nan(trace_capacity,1);
    trace_h(1) = 0;
    trace_violation(1) = abs(init_pesc - sum(init_pes));
    trace_obj(1) = communityPrimalObj(init_pg, init_pb, init_ps, init_pes, init_pesc);
    trace_lambda(1) = init_pi;
    trace_lambda_used(1) = init_pi;
    trace_count = 1;
end

res.success = 1;
res.stop_reason = '';
h = 1;

if progress_print_enabled
    fprintf('[S3 inner start] k=%d, community=%d, EU=%d, CV limit=%.3f kW\n', ...
        outer_iter, community_id, N, inner_cv_tol_kW);
end

while true
    lambda_used_now = lambda;
    for j = 1:N
        if count(j) < delay(j)
            count(j) = count(j) + 1;
        else
            [fast_pg(j), fast_pes(j), fast_pb(j), fast_ps(j)] = ...
                solveBranch(j, lambda_used_now);
            count(j) = 0;
            delay(j) = randi([0,h0]);
        end
    end

    pesc_now = (pi_0 - lambda_used_now) / a;

    if h >= inner_avg_start_iter
        avg_cnt = avg_cnt + 1;
        sum_pes = sum_pes + fast_pes;
        sum_pg  = sum_pg  + fast_pg;
        sum_pb  = sum_pb  + fast_pb;
        sum_ps  = sum_ps  + fast_ps;
        sum_pesc = sum_pesc + pesc_now;
        bar_pes = sum_pes / avg_cnt;
        bar_pg  = sum_pg  / avg_cnt;
        bar_pb  = sum_pb  / avg_cnt;
        bar_ps  = sum_ps  / avg_cnt;
        bar_pesc_raw = sum_pesc / avg_cnt;
    else
        bar_pes = fast_pes;
        bar_pg  = fast_pg;
        bar_pb  = fast_pb;
        bar_ps  = fast_ps;
        bar_pesc_raw = pesc_now;
    end

    % Algorithm 2 原始 Agg 输出：不以 EU 输出和覆盖。
    bar_pesc = bar_pesc_raw;
    lambda_aux = lambda_used_now + alpha * (pesc_now - sum(fast_pes));
    lambda_next = lambda_aux + inner_momentum * (lambda_aux - lambda_aux_prev);
    lambda_step = abs(lambda_next - lambda_used_now);
    inner_cv_bar_now_kW = abs(bar_pesc_raw - sum(bar_pes));

    if record_trace
        trace_count = trace_count + 1;
        if trace_count > trace_capacity
            trace_capacity = growHistoryCapacity(trace_capacity, trace_count, max_inner_iter + 1);
            trace_h(trace_capacity,1) = NaN;
            trace_violation(trace_capacity,1) = NaN;
            trace_obj(trace_capacity,1) = NaN;
            trace_lambda(trace_capacity,1) = NaN;
            trace_lambda_used(trace_capacity,1) = NaN;
        end
        trace_h(trace_count) = h;
        trace_violation(trace_count) = inner_cv_bar_now_kW;
        trace_obj(trace_count) = communityPrimalObj(bar_pg, bar_pb, bar_ps, bar_pes, bar_pesc_raw);
        trace_lambda(trace_count) = lambda_next;
        trace_lambda_used(trace_count) = lambda_used_now;
    end

    if progress_print_enabled && (h == 1 || mod(h, progress_print_inner_every) == 0)
        fprintf('[S3 inner] k=%d, community=%d, h=%d, CVavg=%.3f kW, limit=%.3f kW, lambdaStep=%.3e\n', ...
            outer_iter, community_id, h, inner_cv_bar_now_kW, inner_cv_tol_kW, lambda_step);
    end

    lambda_aux_prev = lambda_aux;
    lambda = lambda_next;

    if h >= min_inner_iter && lambda_step <= errTol && ...
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

if progress_print_enabled
    fprintf('[S3 inner end] k=%d, community=%d, h=%d, CVavg=%.3f kW, reason=%s\n', ...
        outer_iter, community_id, h, inner_cv_bar_now_kW, res.stop_reason);
end

if record_trace
    res.trace_h = trace_h(1:trace_count);
    res.trace_violation = trace_violation(1:trace_count);
    res.trace_obj = trace_obj(1:trace_count);
    res.trace_lambda = trace_lambda(1:trace_count);
    res.trace_lambda_used = trace_lambda_used(1:trace_count);
else
    res.trace_h = [];
    res.trace_violation = [];
    res.trace_obj = [];
    res.trace_lambda = [];
    res.trace_lambda_used = [];
end

res.pi = lambda;
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
res.iter = h;
res.price_updates = h;

    function [pg, pes, pb, ps] = solveBranch(j, lambda_v)
        [pg, pes, pb, ps] = boundedQREBranch(c(j), b(j), a, D(j), ...
            pg_min(j), pg_max(j), lambda_v, pi_max, pi_min, beta_qre(j), ...
            NaN, qre_z_cap, qre_max_backoffs);
    end

    function obj = communityPrimalObj(pg_v, pb_v, ps_v, pes_v, pesc_v)
        obj = sum(0.5 .* c(:) .* pg_v(:).^2 ...
            + b(:) .* pg_v(:) ...
            + pi_max .* pb_v(:) ...
            - pi_min .* ps_v(:) ...
            + 0.5 .* a .* pes_v(:).^2) ...
            + 0.5 .* a .* pesc_v.^2 - pi_0 .* pesc_v;
    end
end
