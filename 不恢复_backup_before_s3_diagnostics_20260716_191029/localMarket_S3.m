function res = localMarket_S3(c,b,a,D,init_rho,pg_max,pg_min,init_pi, ...
    init_pes,init_pg,init_pb,init_ps,init_pesc,errTol,inner_cv_tol_kW, ...
    pi_0,pi_max,pi_min,record_trace,beta_qre,qre_epsilon,agg_epsilon_i, ...
    qre_z_cap,qre_backoff_factor,qre_max_backoffs,inner_momentum, ...
    max_inner_iter,h0,min_inner_iter,inner_avg_start_iter,inner_cv_stop_enabled, ...
    qre_certificate_enabled,agg_certificate_enabled,agg_cert_tol, ...
    stable_inner_window,qre_audit_enabled,qre_audit_rate, ...
    qre_audit_trace_community_full,qre_audit_seed,agg_gap_diagnostic_enabled, ...
    qre_noise_enabled,progress_print_enabled,progress_print_inner_every, ...
    outer_iter,community_id)
%LOCALMARKET_S3 Asynchronous held EU response with suffix-ergodic output.
% Formal large-scale runs use fast projected QRE plus posterior audit.

N = length(b);
validateattributes(init_rho, {'numeric'}, {'scalar','real','finite','positive'});
validateattributes(inner_cv_tol_kW, {'numeric'}, {'scalar','real','finite','nonnegative'});
validateattributes(qre_epsilon, {'numeric'}, {'scalar','real','finite','nonnegative'});
validateattributes(agg_epsilon_i, {'numeric'}, {'scalar','real','finite','nonnegative'});
validateattributes(qre_audit_rate, {'numeric'}, {'scalar','real','finite','nonnegative','<=',1});
validateattributes(h0, {'numeric'}, {'scalar','real','finite','nonnegative','integer'});
validateattributes(min_inner_iter, {'numeric'}, {'scalar','real','finite','integer','nonnegative'});
validateattributes(inner_avg_start_iter, {'numeric'}, {'scalar','real','finite','positive','integer'});
validateattributes(stable_inner_window, {'numeric'}, {'scalar','real','finite','positive','integer'});
if nargin < 43
    progress_print_enabled = false;
    progress_print_inner_every = inf;
    outer_iter = NaN;
    community_id = NaN;
end

rho = init_rho;
alpha = rho * a;
delay = randi([0,h0],N,1);
count = zeros(N,1);

lambda = init_pi;
lambda_aux_prev = init_pi;
fast_pg = init_pg;
fast_pes = init_pes;
fast_pb = init_pb;
fast_ps = init_ps;

sum_pes = zeros(N,1);
sum_pg = zeros(N,1);
sum_pb = zeros(N,1);
sum_ps = zeros(N,1);
sum_pesc = 0;
avg_cnt = 0;
bar_pes = init_pes;
bar_pg = init_pg;
bar_pb = init_pb;
bar_ps = init_ps;
bar_pesc_raw = init_pesc;

[obj_comm_exact,obj_comm_meta] = exactLowerObjFast( ...
    c,b,a,D,pg_max,pg_min,pi_0,pi_max,pi_min);

% Audit statistics are bounded: no EU-by-iteration history is retained.
qre_audit_count = 0;
qre_audit_gap_sum = 0;
qre_audit_gap_max = 0;
qre_audit_gap_exceed_count = 0;
qre_audit_branch_switch_count = 0;
qre_audit_response_count = 0;
qre_audit_gap_samples = zeros(0,1);
qre_audit_sample_capacity = 10000;
qre_response_count = 0;
qre_gap_max = NaN;
qre_gap_mean = NaN;
qre_gap_p95 = NaN;
qre_branch_switch_count = 0;
qre_branch_switch_rate = NaN;
audit_period = max(1,round(1 / max(qre_audit_rate,eps)));
audit_full_community = qre_audit_trace_community_full && ...
    logical(record_trace);
audit_offset = mod(qre_audit_seed + max(community_id,0),audit_period);

stable_price_count = 0;

if record_trace
    trace_capacity = min(max_inner_iter + 1,128);
    trace_h = nan(trace_capacity,1);
    trace_violation = nan(trace_capacity,1);
    trace_obj = nan(trace_capacity,1);
    trace_lambda = nan(trace_capacity,1);
    trace_lambda_used = nan(trace_capacity,1);
    trace_h(1) = 0;
    trace_violation(1) = abs(init_pesc - sum(init_pes));
    trace_obj(1) = communityPrimalObj(init_pg,init_pb,init_ps,init_pes,init_pesc);
    trace_lambda(1) = init_pi;
    trace_lambda_used(1) = init_pi;
    trace_count = 1;
end

res.success = 1;
res.stop_reason = '';
h = 1;

if progress_print_enabled
    fprintf('[S3 inner start] k=%d, community=%d, EU=%d, CV limit=%.3f kW\n', ...
        outer_iter,community_id,N,inner_cv_tol_kW);
end

while true
    lambda_used_now = lambda;
    for j = 1:N
        if count(j) < delay(j)
            count(j) = count(j) + 1;
        else
            if qre_certificate_enabled
                qre = solveCertifiedQREResponse( ...
                    c(j),b(j),a,D(j),pg_min(j),pg_max(j),lambda_used_now, ...
                    pi_max,pi_min,beta_qre(j),qre_epsilon,qre_z_cap, ...
                    qre_backoff_factor,qre_max_backoffs);
            else
                qre = solveFastQREResponse( ...
                    c(j),b(j),a,D(j),pg_min(j),pg_max(j),lambda_used_now, ...
                    pi_max,pi_min,beta_qre(j),qre_z_cap,qre_noise_enabled);
            end

            fast_pg(j) = qre.pg;
            fast_pes(j) = qre.pes;
            fast_pb(j) = qre.pb;
            fast_ps(j) = qre.ps;
            qre_response_count = qre_response_count + 1;

            audit_due = qre_audit_enabled && ...
                (audit_full_community || ...
                mod(qre_response_count + audit_offset,audit_period) == 0);
            if audit_due
                exact = solveExactEUResponse( ...
                    c(j),b(j),a,D(j),pg_min(j),pg_max(j),lambda_used_now, ...
                    pi_max,pi_min);
                gap_raw = qre.obj - exact.obj;
                cert_tol = 1e-10 * max(1,abs(exact.obj));
                if gap_raw < -cert_tol
                    error('localMarket_S3:objectiveInconsistency', ...
                        'Exact EU response or objective implementation is inconsistent.');
                end
                gap = max(gap_raw,0);
                qre_audit_count = qre_audit_count + 1;
                qre_audit_gap_sum = qre_audit_gap_sum + gap;
                qre_audit_gap_max = max(qre_audit_gap_max,gap);
                qre_audit_gap_exceed_count = qre_audit_gap_exceed_count + ...
                    double(gap > qre_epsilon + cert_tol);
                qre_audit_branch_switch_count = qre_audit_branch_switch_count + ...
                    double(qre.branch ~= exact.branch);
                qre_audit_response_count = qre_audit_response_count + 1;
                if numel(qre_audit_gap_samples) < qre_audit_sample_capacity
                    qre_audit_gap_samples(end+1,1) = gap;
                else
                    replace_idx = mod(qre_audit_count - 1,qre_audit_sample_capacity) + 1;
                    qre_audit_gap_samples(replace_idx) = gap;
                end
            end

            count(j) = 0;
            delay(j) = randi([0,h0]);
        end
    end

    % Raw Agg is used for the instantaneous lambda update and is never
    % replaced by the EU sum.
    pesc_now = (pi_0 - lambda_used_now) / a;
    if h >= inner_avg_start_iter
        avg_cnt = avg_cnt + 1;
        sum_pes = sum_pes + fast_pes;
        sum_pg = sum_pg + fast_pg;
        sum_pb = sum_pb + fast_pb;
        sum_ps = sum_ps + fast_ps;
        sum_pesc = sum_pesc + pesc_now;
        bar_pes = sum_pes / avg_cnt;
        bar_pg = sum_pg / avg_cnt;
        bar_pb = sum_pb / avg_cnt;
        bar_ps = sum_ps / avg_cnt;
        bar_pesc_raw = sum_pesc / avg_cnt;
    else
        bar_pes = fast_pes;
        bar_pg = fast_pg;
        bar_pb = fast_pb;
        bar_ps = fast_ps;
        bar_pesc_raw = pesc_now;
    end

    lambda_aux = lambda_used_now + alpha * (pesc_now - sum(fast_pes));
    lambda_next = lambda_aux + inner_momentum * (lambda_aux - lambda_aux_prev);
    lambda_step = abs(lambda_next - lambda_used_now);
    inner_cv_bar_now_kW = abs(bar_pesc_raw - sum(bar_pes));

    if agg_gap_diagnostic_enabled
        obj_comm_out = communityPrimalObj(bar_pg,bar_pb,bar_ps,bar_pes,bar_pesc_raw);
        raw_agg_output_objective_gap = max(obj_comm_out - obj_comm_exact,0);
    else
        obj_comm_out = NaN;
        raw_agg_output_objective_gap = NaN;
    end

    if lambda_step <= errTol
        stable_price_count = stable_price_count + 1;
    else
        stable_price_count = 0;
    end

    if record_trace
        trace_count = trace_count + 1;
        if trace_count > trace_capacity
            trace_capacity = growHistoryCapacity(trace_capacity,trace_count,max_inner_iter + 1);
            trace_h(trace_capacity,1) = NaN;
            trace_violation(trace_capacity,1) = NaN;
            trace_obj(trace_capacity,1) = NaN;
            trace_lambda(trace_capacity,1) = NaN;
            trace_lambda_used(trace_capacity,1) = NaN;
        end
        trace_h(trace_count) = h;
        trace_violation(trace_count) = inner_cv_bar_now_kW;
        trace_obj(trace_count) = obj_comm_out;
        trace_lambda(trace_count) = lambda_next;
        trace_lambda_used(trace_count) = lambda_used_now;
    end

    if progress_print_enabled && (h == 1 || mod(h,progress_print_inner_every) == 0)
        fprintf('[S3 inner] k=%d, community=%d, h=%d, CV=%.3f/%.3f kW, ', ...
            outer_iter,community_id,h,inner_cv_bar_now_kW,inner_cv_tol_kW);
        fprintf('rawAggGap=%.4g, lambdaStep=%.3e, auditMax=%.4g\n', ...
            raw_agg_output_objective_gap,lambda_step,qre_audit_gap_max);
    end

    lambda_aux_prev = lambda_aux;
    lambda = lambda_next;

    inner_cv_pass = (~logical(inner_cv_stop_enabled)) || ...
        (inner_cv_bar_now_kW <= inner_cv_tol_kW);
    stable_pass = stable_price_count >= stable_inner_window;
    formal_sample = h >= inner_avg_start_iter;
    if h >= min_inner_iter && formal_sample && stable_pass && inner_cv_pass
        res.stop_reason = 'stable_price_and_relative_cv';
        break;
    end
    if isfinite(max_inner_iter) && h >= max_inner_iter
        res.success = 0;
        res.stop_reason = 'max_inner_iter';
        break;
    end
    h = h + 1;
end

if progress_print_enabled
    fprintf('[S3 inner end] k=%d, community=%d, h=%d, CV=%.3f, reason=%s\n', ...
        outer_iter,community_id,h,inner_cv_bar_now_kW,res.stop_reason);
end

if record_trace
    res.trace_h = trace_h(1:trace_count);
    res.trace_violation = trace_violation(1:trace_count);
    res.trace_obj = trace_obj(1:trace_count);
    res.trace_lambda = trace_lambda(1:trace_count);
    res.trace_lambda_used = trace_lambda_used(1:trace_count);
    res.trace_price_update_count = res.trace_h;
    res.trace_update_event = [false; true(max(numel(res.trace_h)-1,0),1)];
else
    res.trace_h = [];
    res.trace_violation = [];
    res.trace_obj = [];
    res.trace_lambda = [];
    res.trace_lambda_used = [];
    res.trace_price_update_count = [];
    res.trace_update_event = [];
end

if isempty(qre_audit_gap_samples)
    qre_gap_p95 = NaN;
else
    qre_gap_p95 = prctile(qre_audit_gap_samples,95);
end
qre_gap_max = qre_audit_gap_max;
qre_gap_mean = qre_audit_gap_sum / max(qre_audit_count,1);
qre_branch_switch_count = qre_audit_branch_switch_count;
qre_branch_switch_rate = qre_branch_switch_count / max(qre_audit_count,1);

res.pi = lambda;
res.pes = bar_pes;
res.pg = bar_pg;
res.pb = bar_pb;
res.ps = bar_ps;
res.pesc = bar_pesc_raw;
res.pesc_agg_raw = bar_pesc_raw;
res.pes_sum = sum(bar_pes);
res.inner_cv_kW = abs(bar_pesc_raw - res.pes_sum);
res.inner_avg_cnt = avg_cnt;
res.inner_avg_start_iter = inner_avg_start_iter;
res.formal_output_type = 'suffix_ergodic_average';
res.qre_mode = ternary(qre_certificate_enabled,'certified','fast_qre_with_posterior_audit');
res.qre_certificate_enabled = logical(qre_certificate_enabled);
res.qre_audit_enabled = logical(qre_audit_enabled);
res.qre_audit_rate = qre_audit_rate;
res.qre_audit_seed = qre_audit_seed;
res.qre_audit_count = qre_audit_count;
res.qre_gap_max = qre_gap_max;
res.qre_gap_mean = qre_gap_mean;
res.qre_gap_p95 = qre_gap_p95;
res.qre_gap_exceed_count = qre_audit_gap_exceed_count;
res.qre_gap_exceed_rate = qre_audit_gap_exceed_count / max(qre_audit_count,1);
res.qre_branch_switch_count = qre_branch_switch_count;
res.qre_branch_switch_rate = qre_branch_switch_rate;
res.qre_response_count = qre_response_count;
res.qre_audit_response_count = qre_audit_response_count;
res.qre_all_certificate_pass = qre_audit_gap_exceed_count == 0;
res.qre_gamma_min = NaN;
res.qre_gamma_mean = NaN;
res.qre_backoff_total = 0;
res.qre_backoff_max = 0;
res.qre_fallback_count = 0;
res.agg_certificate_enabled = logical(agg_certificate_enabled);
res.agg_certificate_pass = NaN;
res.agg_gap = raw_agg_output_objective_gap;
res.agg_epsilon = agg_epsilon_i;
res.agg_cert_tol = agg_cert_tol;
res.raw_agg_output_objective_gap = raw_agg_output_objective_gap;
res.raw_agg_reference_diagnostic = raw_agg_output_objective_gap;
res.obj_comm_exact = obj_comm_exact;
res.obj_comm_out = obj_comm_out;
res.obj_comm_exact_meta = obj_comm_meta;
res.agg_gap_diagnostic_enabled = logical(agg_gap_diagnostic_enabled);
res.stable_inner_window = stable_inner_window;
res.stable_price_count = stable_price_count;
res.inner_cv_stop_enabled = logical(inner_cv_stop_enabled);
res.iter = h;
res.price_updates = h;

    function obj = communityPrimalObj(pg_v,pb_v,ps_v,pes_v,pesc_v)
        obj = sum(0.5 .* c(:) .* pg_v(:).^2 + b(:) .* pg_v(:) ...
            + pi_max .* pb_v(:) - pi_min .* ps_v(:) ...
            + 0.5 .* a .* pes_v(:).^2) ...
            + 0.5 .* a .* pesc_v.^2 - pi_0 .* pesc_v;
    end

    function value = ternary(condition,true_value,false_value)
        if condition
            value = true_value;
        else
            value = false_value;
        end
    end
end
