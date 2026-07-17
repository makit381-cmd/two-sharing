function res = localMarket_S3(c,b,a,D,init_rho,pg_max,pg_min,init_pi, ...
    init_pes,init_pg,init_pb,init_ps,init_pesc,errTol,inner_cv_tol_kW, ...
    pi_0,pi_max,pi_min,record_trace,beta_qre,qre_epsilon,agg_epsilon_i, ...
    qre_z_cap,qre_backoff_factor,qre_max_backoffs,inner_momentum, ...
    max_inner_iter,h0,min_inner_iter,inner_avg_start_iter,inner_cv_stop_enabled, ...
    qre_certificate_enabled,agg_certificate_enabled,agg_cert_tol, ...
    stable_inner_window,qre_audit_enabled,qre_audit_rate, ...
    qre_audit_trace_community_full,qre_audit_seed,agg_gap_diagnostic_enabled, ...
    qre_noise_enabled,diagnostic_record_every,exact_sync_diagnostic_every, ...
    rolling_window,progress_print_enabled,progress_print_inner_every, ...
    outer_iter,community_id,rng_config,avg_config,local_call_count)
%LOCALMARKET_S3 Asynchronous held EU response with suffix-ergodic output.
% Formal large-scale runs use fast projected QRE plus posterior audit.

N = length(b);
validateattributes(init_rho, {'numeric'}, {'scalar','real','finite','positive'});
validateattributes(inner_cv_tol_kW, {'numeric'}, {'scalar','real','finite','nonnegative'});
validateattributes(qre_epsilon, {'numeric'}, {'scalar','real','finite','nonnegative'});
validateattributes(agg_epsilon_i, {'numeric'}, {'scalar','real','finite','nonnegative'});
validateattributes(qre_audit_rate, {'numeric'}, {'scalar','real','finite','nonnegative','<=',1});
validateattributes(diagnostic_record_every, {'numeric'}, {'scalar','real','finite','positive','integer'});
validateattributes(exact_sync_diagnostic_every, {'numeric'}, {'scalar','real','finite','positive','integer'});
validateattributes(rolling_window, {'numeric'}, {'scalar','real','finite','positive','integer'});
validateattributes(h0, {'numeric'}, {'scalar','real','finite','nonnegative','integer'});
validateattributes(min_inner_iter, {'numeric'}, {'scalar','real','finite','integer','nonnegative'});
validateattributes(inner_avg_start_iter, {'numeric'}, {'scalar','real','finite','positive','integer'});
validateattributes(stable_inner_window, {'numeric'}, {'scalar','real','finite','positive','integer'});
if nargin < 49 || isempty(rng_config)
    rng_config = struct();
end
if nargin < 50 || isempty(avg_config)
    avg_config = struct();
end
if nargin < 51 || isempty(local_call_count)
    local_call_count = 1;
end
if nargin < 43
    progress_print_enabled = false;
    progress_print_inner_every = inf;
    outer_iter = NaN;
    community_id = NaN;
end

if ~isfield(rng_config,'inner_delay_seed'), rng_config.inner_delay_seed = 31001; end
if ~isfield(rng_config,'qre_noise_seed'), rng_config.qre_noise_seed = 31002; end
if ~isfield(rng_config,'qre_audit_seed'), rng_config.qre_audit_seed = 31003; end
if ~isfield(avg_config,'policy'), avg_config.policy = 'fixed_inner_start'; end
if ~isfield(avg_config,'policy_version'), avg_config.policy_version = 'A2_dynamic_stable_start_v1'; end
if ~isfield(avg_config,'start_cv_factor'), avg_config.start_cv_factor = 5; end
if ~isfield(avg_config,'start_price_factor'), avg_config.start_price_factor = 10; end
if ~isfield(avg_config,'start_stable_window'), avg_config.start_stable_window = 20; end
if ~isfield(avg_config,'min_price_updates'), avg_config.min_price_updates = 100; end
if ~isfield(avg_config,'min_samples'), avg_config.min_samples = 200; end
if ~isfield(avg_config,'formal_cv_stable_window'), avg_config.formal_cv_stable_window = 5; end

validateattributes(avg_config.start_cv_factor, {'numeric'}, {'scalar','real','finite','nonnegative'});
validateattributes(avg_config.start_price_factor, {'numeric'}, {'scalar','real','finite','nonnegative'});
validateattributes(avg_config.start_stable_window, {'numeric'}, {'scalar','real','finite','positive','integer'});
validateattributes(avg_config.min_price_updates, {'numeric'}, {'scalar','real','finite','nonnegative','integer'});
validateattributes(avg_config.min_samples, {'numeric'}, {'scalar','real','finite','nonnegative','integer'});
validateattributes(avg_config.formal_cv_stable_window, {'numeric'}, ...
    {'scalar','real','finite','positive','integer'});
validateattributes(local_call_count, {'numeric'}, ...
    {'scalar','real','finite','positive','integer'});
use_dynamic_average = strcmp(avg_config.policy,'dynamic_stable_start');

community_key = 0;
if isfinite(community_id), community_key = max(0,floor(community_id)); end
inner_delay_seed = deriveDeterministicSeed( ...
    rng_config.inner_delay_seed,community_key,local_call_count,1);
qre_noise_seed = deriveDeterministicSeed( ...
    rng_config.qre_noise_seed,community_key,local_call_count,2);
qre_audit_seed_effective = deriveDeterministicSeed( ...
    rng_config.qre_audit_seed,community_key,local_call_count,3);
delay_stream = RandStream('mt19937ar','Seed',inner_delay_seed);
qre_stream = RandStream('mt19937ar','Seed',qre_noise_seed);
audit_stream = RandStream('mt19937ar','Seed',qre_audit_seed_effective);
hash_horizon = 10000;
[delay_sequence_hash,update_sequence_hash] = computeScheduleHashes( ...
    inner_delay_seed,h0,N,hash_horizon);

diagnostic_record_every = max(1,diagnostic_record_every);
exact_sync_diagnostic_every = max(1,exact_sync_diagnostic_every);
rolling_window = max(1,rolling_window);

rho = init_rho;
alpha = rho * a;
delay = randi(delay_stream,[0,h0],N,1);
delay_hash_state = 2166136261;
for j = 1:N
    delay_hash_state = updateHash(delay_hash_state,delay(j) + 1);
end
update_sequence_hash_state = 2166136261;
count = zeros(N,1);
age_since_last_update = zeros(N,1);
update_count = zeros(N,1);
max_observed_age = 0;
sum_observed_age = 0;
age_observation_count = 0;

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
price_update_count = 0;
avg_start_stable_count = 0;
averaging_started = false;
formal_average_available = false;
inner_avg_start_wall_iter = NaN;
inner_avg_start_price_update = NaN;
formal_cv_kW = NaN;
formal_cv_pass = false;
formal_cv_stable_pass = false;
recent_lambda_steps = nan(max(1,stable_inner_window),1);
recent_step_count = 0;
recent_step_position = 0;
recent_formal_cv_kW = nan(avg_config.formal_cv_stable_window,1);
recent_formal_cv_count = 0;
recent_formal_cv_position = 0;
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
qre_noise_hash_state = 2166136261;
qre_gap_max = NaN;
qre_gap_mean = NaN;
qre_gap_p95 = NaN;
qre_branch_switch_count = 0;
qre_branch_switch_rate = NaN;
audit_full_community = qre_audit_trace_community_full && ...
    logical(record_trace);

lambda_exact_result = struct('lambda_exact',NaN,'residual_exact',NaN, ...
    'success',false,'iter',0,'lower_bound',NaN,'upper_bound',NaN);
if record_trace
    lambda_exact_result = solveExactCommunityLambda( ...
        c,b,a,D,pg_min,pg_max,pi_0,pi_max,pi_min,init_pi);
end

rolling_pesc_buffer = zeros(rolling_window,1);
rolling_pes_sum_buffer = zeros(rolling_window,1);
rolling_count = 0;
rolling_pos = 0;
last_exact_sync_cv = NaN;
last_exact_sync_h = NaN;
trace_exact_sync_h = zeros(0,1);
trace_exact_sync_cv_kW = zeros(0,1);

stable_price_count = 0;

if record_trace
    trace_capacity = max(128,min(ceil(max_inner_iter / diagnostic_record_every) + 2,10000));
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
    trace_lambda_step = nan(trace_capacity,1);
    trace_cv_held_kW = nan(trace_capacity,1);
    trace_cv_exact_sync_kW = nan(trace_capacity,1);
    trace_cv_avg_kW = nan(trace_capacity,1);
    trace_cv_rolling_kW = nan(trace_capacity,1);
    trace_lambda_error = nan(trace_capacity,1);
    trace_max_age = nan(trace_capacity,1);
    trace_mean_age = nan(trace_capacity,1);
    trace_averaging_started = false(trace_capacity,1);
    trace_avg_cnt = zeros(trace_capacity,1);
    trace_avg_start_stable_count = zeros(trace_capacity,1);
    trace_formal_cv_stable_pass = false(trace_capacity,1);
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
    age_since_last_update = age_since_last_update + 1;
    lambda_used_now = lambda;
    for j = 1:N
        if count(j) < delay(j)
            count(j) = count(j) + 1;
        else
            if qre_certificate_enabled
                qre = solveCertifiedQREResponse( ...
                    c(j),b(j),a,D(j),pg_min(j),pg_max(j),lambda_used_now, ...
                    pi_max,pi_min,beta_qre(j),qre_epsilon,qre_z_cap, ...
                    qre_backoff_factor,qre_max_backoffs,qre_stream);
            else
                qre = solveFastQREResponse( ...
                    c(j),b(j),a,D(j),pg_min(j),pg_max(j),lambda_used_now, ...
                    pi_max,pi_min,beta_qre(j),qre_z_cap,qre_noise_enabled,qre_stream);
            end

            fast_pg(j) = qre.pg;
            fast_pes(j) = qre.pes;
            fast_pb(j) = qre.pb;
            fast_ps(j) = qre.ps;
            qre_response_count = qre_response_count + 1;
            qre_noise_hash_state = updateHash(qre_noise_hash_state, ...
                round((qre.z_pg_raw + qre_z_cap) * 1e9));
            qre_noise_hash_state = updateHash(qre_noise_hash_state, ...
                round((qre.z_pes_raw + qre_z_cap) * 1e9));

            audit_due = qre_audit_enabled && ...
                (audit_full_community || rand(audit_stream,1) <= qre_audit_rate);
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
                audit_exceed = gap > qre_epsilon + cert_tol;
                qre_audit_gap_exceed_count = qre_audit_gap_exceed_count + ...
                    double(audit_exceed);
                if audit_exceed && qre_audit_gap_exceed_count == 1
                    warning('localMarket_S3:qreAuditExceed', ...
                        'Posterior QRE audit exceeded qre_epsilon; formal response was not modified.');
                end
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
            delay(j) = randi(delay_stream,[0,h0]);
            delay_hash_state = updateHash(delay_hash_state,delay(j) + 1);
            age_since_last_update(j) = 0;
            update_count(j) = update_count(j) + 1;
            update_sequence_hash_state = updateHash(update_sequence_hash_state,h);
            update_sequence_hash_state = updateHash(update_sequence_hash_state,j);
        end
    end

    max_observed_age = max(max_observed_age,max(age_since_last_update));
    sum_observed_age = sum_observed_age + sum(age_since_last_update);
    age_observation_count = age_observation_count + N;

    % Raw Agg is used for the instantaneous lambda update and is never
    % replaced by the EU sum.
    pesc_now = (pi_0 - lambda_used_now) / a;
    residual_now = pesc_now - sum(fast_pes);
    residual_held_kW = residual_now;
    cv_held_kW = abs(residual_held_kW);
    lambda_aux = lambda_used_now + alpha * residual_now;
    lambda_next = lambda_aux + inner_momentum * (lambda_aux - lambda_aux_prev);
    lambda_step = abs(lambda_next - lambda_used_now);

    % S3 has one real price update per wall slot. Waiting is not represented
    % by an artificial zero step in this inner routine.
    price_update_count = price_update_count + 1;
    recent_step_count = min(recent_step_count + 1,numel(recent_lambda_steps));
    recent_step_position = mod(recent_step_position,numel(recent_lambda_steps)) + 1;
    recent_lambda_steps(recent_step_position) = lambda_step;
    if lambda_step <= errTol
        stable_price_count = stable_price_count + 1;
    else
        stable_price_count = 0;
    end

    average_active_before_slot = averaging_started;
    if use_dynamic_average
        held_cv_pass_for_avg = cv_held_kW <= ...
            avg_config.start_cv_factor * inner_cv_tol_kW;
        price_step_pass_for_avg = lambda_step <= ...
            avg_config.start_price_factor * errTol;
        if price_update_count >= avg_config.min_price_updates && ...
                held_cv_pass_for_avg && price_step_pass_for_avg
            avg_start_stable_count = avg_start_stable_count + 1;
        else
            avg_start_stable_count = 0;
        end
        if ~averaging_started && ...
                price_update_count >= avg_config.min_price_updates && ...
                avg_start_stable_count >= avg_config.start_stable_window
            averaging_started = true;
            inner_avg_start_wall_iter = h;
            inner_avg_start_price_update = price_update_count;
            sum_pes = zeros(N,1);
            sum_pg = zeros(N,1);
            sum_pb = zeros(N,1);
            sum_ps = zeros(N,1);
            sum_pesc = 0;
            avg_cnt = 0;
            formal_average_available = false;
        end
    else
        if ~averaging_started && h >= inner_avg_start_iter
            averaging_started = true;
            inner_avg_start_wall_iter = h;
            inner_avg_start_price_update = price_update_count;
        end
    end

    if average_active_before_slot
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
    formal_average_available = averaging_started && ...
        avg_cnt >= avg_config.min_samples;
    residual_avg_kW = bar_pesc_raw - sum(bar_pes);
    cv_avg_kW = abs(residual_avg_kW);
    inner_cv_bar_now_kW = cv_avg_kW;
    if formal_average_available
        formal_cv_kW = cv_avg_kW;
        recent_formal_cv_count = min(recent_formal_cv_count + 1, ...
            numel(recent_formal_cv_kW));
        recent_formal_cv_position = mod(recent_formal_cv_position, ...
            numel(recent_formal_cv_kW)) + 1;
        recent_formal_cv_kW(recent_formal_cv_position) = formal_cv_kW;
    else
        formal_cv_kW = NaN;
    end
    formal_cv_pass = formal_average_available && formal_cv_kW <= inner_cv_tol_kW;
    formal_cv_stable_pass = recent_formal_cv_count >= ...
        avg_config.formal_cv_stable_window && ...
        all(recent_formal_cv_kW <= inner_cv_tol_kW);

    rolling_pos = mod(rolling_pos,rolling_window) + 1;
    rolling_pesc_buffer(rolling_pos) = pesc_now;
    rolling_pes_sum_buffer(rolling_pos) = sum(fast_pes);
    rolling_count = min(rolling_count + 1,rolling_window);
    rolling_pesc_avg = sum(rolling_pesc_buffer(1:rolling_count)) / rolling_count;
    rolling_pes_sum_avg = sum(rolling_pes_sum_buffer(1:rolling_count)) / rolling_count;
    cv_rolling_kW = abs(rolling_pesc_avg - rolling_pes_sum_avg);

    if record_trace && (h == 1 || mod(h,exact_sync_diagnostic_every) == 0)
        pes_exact_sync = zeros(N,1);
        for j = 1:N
            exact_sol = solveExactEUResponse( ...
                c(j),b(j),a,D(j),pg_min(j),pg_max(j), ...
                lambda_used_now,pi_max,pi_min);
            pes_exact_sync(j) = exact_sol.pes;
        end
        residual_exact_sync_kW = pesc_now - sum(pes_exact_sync);
        last_exact_sync_cv = abs(residual_exact_sync_kW);
        last_exact_sync_h = h;
        trace_exact_sync_h(end+1,1) = h;
        trace_exact_sync_cv_kW(end+1,1) = last_exact_sync_cv;
    end

    if agg_gap_diagnostic_enabled
        obj_comm_out = communityPrimalObj(bar_pg,bar_pb,bar_ps,bar_pes,bar_pesc_raw);
        raw_agg_output_objective_gap = max(obj_comm_out - obj_comm_exact,0);
    else
        obj_comm_out = NaN;
        raw_agg_output_objective_gap = NaN;
    end

    if record_trace && (h == 1 || mod(h,diagnostic_record_every) == 0)
        trace_count = trace_count + 1;
        if trace_count > trace_capacity
            trace_capacity = growHistoryCapacity(trace_capacity,trace_count,max_inner_iter + 1);
            trace_h(trace_capacity,1) = NaN;
            trace_violation(trace_capacity,1) = NaN;
            trace_obj(trace_capacity,1) = NaN;
            trace_lambda(trace_capacity,1) = NaN;
            trace_lambda_used(trace_capacity,1) = NaN;
            trace_lambda_step(trace_capacity,1) = NaN;
            trace_cv_held_kW(trace_capacity,1) = NaN;
            trace_cv_exact_sync_kW(trace_capacity,1) = NaN;
            trace_cv_avg_kW(trace_capacity,1) = NaN;
            trace_cv_rolling_kW(trace_capacity,1) = NaN;
            trace_lambda_error(trace_capacity,1) = NaN;
            trace_max_age(trace_capacity,1) = NaN;
            trace_mean_age(trace_capacity,1) = NaN;
            trace_averaging_started(trace_capacity,1) = false;
            trace_avg_cnt(trace_capacity,1) = NaN;
            trace_avg_start_stable_count(trace_capacity,1) = NaN;
            trace_formal_cv_stable_pass(trace_capacity,1) = false;
        end
        trace_h(trace_count) = h;
        trace_violation(trace_count) = inner_cv_bar_now_kW;
        trace_obj(trace_count) = obj_comm_out;
        trace_lambda(trace_count) = lambda_used_now;
        trace_lambda_used(trace_count) = lambda_used_now;
        trace_lambda_step(trace_count) = lambda_step;
        trace_cv_held_kW(trace_count) = cv_held_kW;
        trace_cv_exact_sync_kW(trace_count) = last_exact_sync_cv;
        trace_cv_avg_kW(trace_count) = cv_avg_kW;
        trace_cv_rolling_kW(trace_count) = cv_rolling_kW;
        if lambda_exact_result.success
            trace_lambda_error(trace_count) = abs(lambda_used_now - lambda_exact_result.lambda_exact);
        end
        trace_max_age(trace_count) = max(age_since_last_update);
        trace_mean_age(trace_count) = mean(age_since_last_update);
        trace_averaging_started(trace_count) = averaging_started;
        trace_avg_cnt(trace_count) = avg_cnt;
        trace_avg_start_stable_count(trace_count) = avg_start_stable_count;
        trace_formal_cv_stable_pass(trace_count) = formal_cv_stable_pass;
    end

    if progress_print_enabled && (h == 1 || mod(h,progress_print_inner_every) == 0)
        fprintf('[S3 inner] k=%d, community=%d, h=%d, CV=%.3f/%.3f kW, ', ...
            outer_iter,community_id,h,inner_cv_bar_now_kW,inner_cv_tol_kW);
        fprintf('rawAggGap=%.4g, lambdaStep=%.3e, auditMax=%.4g\n', ...
            raw_agg_output_objective_gap,lambda_step,qre_audit_gap_max);
    end

    lambda_aux_prev = lambda_aux;
    lambda = lambda_next;

    stable_pass = recent_step_count >= stable_inner_window && ...
        all(recent_lambda_steps <= errTol);
    if use_dynamic_average
        if h >= min_inner_iter && formal_average_available && ...
                stable_pass && formal_cv_pass && formal_cv_stable_pass
            res.stop_reason = 'dynamic_average_price_step_and_relative_cv';
            break;
        end
    else
        inner_cv_pass = (~logical(inner_cv_stop_enabled)) || ...
            (inner_cv_bar_now_kW <= inner_cv_tol_kW);
        formal_sample = h >= inner_avg_start_iter;
        if h >= min_inner_iter && formal_sample && stable_pass && inner_cv_pass
            res.stop_reason = 'stable_price_and_relative_cv';
            break;
        end
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
    res.trace_lambda_step = trace_lambda_step(1:trace_count);
    res.trace_cv_held_kW = trace_cv_held_kW(1:trace_count);
    res.trace_cv_exact_sync_kW = trace_cv_exact_sync_kW(1:trace_count);
    res.trace_cv_avg_kW = trace_cv_avg_kW(1:trace_count);
    res.trace_cv_rolling_kW = trace_cv_rolling_kW(1:trace_count);
    res.trace_lambda_error = trace_lambda_error(1:trace_count);
    res.trace_max_age = trace_max_age(1:trace_count);
    res.trace_mean_age = trace_mean_age(1:trace_count);
    res.trace_averaging_started = trace_averaging_started(1:trace_count);
    res.trace_avg_cnt = trace_avg_cnt(1:trace_count);
    res.trace_avg_start_stable_count = trace_avg_start_stable_count(1:trace_count);
    res.trace_formal_cv_stable_pass = trace_formal_cv_stable_pass(1:trace_count);
    res.trace_price_update_count = res.trace_h;
    res.trace_update_event = [false; true(max(numel(res.trace_h)-1,0),1)];
else
    res.trace_h = [];
    res.trace_violation = [];
    res.trace_obj = [];
    res.trace_lambda = [];
    res.trace_lambda_used = [];
    res.trace_lambda_step = [];
    res.trace_cv_held_kW = [];
    res.trace_cv_exact_sync_kW = [];
    res.trace_cv_rolling_kW = [];
    res.trace_cv_avg_kW = [];
    res.trace_lambda_error = [];
    res.trace_max_age = [];
    res.trace_mean_age = [];
    res.trace_averaging_started = [];
    res.trace_avg_cnt = [];
    res.trace_avg_start_stable_count = [];
    res.trace_formal_cv_stable_pass = [];
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
res.residual_held_final_kW = residual_held_kW;
res.cv_held_final_kW = cv_held_kW;
res.residual_avg_final_kW = residual_avg_kW;
res.cv_avg_final_kW = cv_avg_kW;
res.cv_rolling_final_kW = cv_rolling_kW;
res.cv_exact_sync_final_kW = last_exact_sync_cv;
res.trace_exact_sync_h = trace_exact_sync_h;
res.trace_exact_sync_cv_kW = trace_exact_sync_cv_kW;
res.lambda_exact = lambda_exact_result.lambda_exact;
res.lambda_exact_result = lambda_exact_result;
res.final_lambda_error = abs(lambda - lambda_exact_result.lambda_exact);
res.inner_avg_policy = avg_config.policy;
res.inner_avg_policy_version = avg_config.policy_version;
res.averaging_started = averaging_started;
res.inner_avg_start_wall_iter = inner_avg_start_wall_iter;
res.inner_avg_start_price_update = inner_avg_start_price_update;
res.inner_avg_start_cv_factor = avg_config.start_cv_factor;
res.inner_avg_start_price_factor = avg_config.start_price_factor;
res.inner_avg_start_stable_window = avg_config.start_stable_window;
res.inner_avg_min_price_updates = avg_config.min_price_updates;
res.inner_avg_min_samples = avg_config.min_samples;
res.inner_formal_cv_stable_window = avg_config.formal_cv_stable_window;
res.formal_average_available = formal_average_available;
res.final_formal_cv_kW = formal_cv_kW;
res.formal_cv_pass = formal_cv_pass;
res.formal_cv_stable_pass = formal_cv_stable_pass;
res.final_held_cv_kW = cv_held_kW;
res.final_lambda_step = lambda_step;
res.price_update_count = price_update_count;
res.wall_iter = h;
res.avg_start_stable_count = avg_start_stable_count;
res.inner_delay_seed = inner_delay_seed;
res.qre_noise_seed = qre_noise_seed;
res.qre_audit_seed = qre_audit_seed_effective;
res.qre_noise_sequence_hash = sprintf('%.0f',qre_noise_hash_state);
res.delay_sequence_hash = delay_sequence_hash;
res.update_sequence_hash = update_sequence_hash;
res.executed_delay_sequence_hash = sprintf('%.0f',delay_hash_state);
res.executed_update_sequence_hash = sprintf('%.0f',update_sequence_hash_state);
res.delay_hash_horizon = hash_horizon;
res.community_id = community_id;
res.local_call_count = local_call_count;
res.seed_derivation_policy = 'deriveDeterministicSeed(base_seed,community_id,local_call_count,stream_tag)';
res.warm_start_initialized = true;
res.initialization_semantics = 'fast_pg/pes/pb/ps and lambda use init_* warm start; pesc raw is recomputed from pi_0 and lambda';
res.max_observed_EU_age = max_observed_age;
res.mean_observed_EU_age = sum_observed_age / max(age_observation_count,1);
res.min_update_count = min(update_count);
res.max_update_count = max(update_count);
res.all_EU_updated = all(update_count > 0);
res.age_bound_h0 = h0;
res.age_bound_pass_h0 = max_observed_age <= h0;
res.age_bound_pass_h0_plus_1 = max_observed_age <= h0 + 1;
res.diagnostic_record_every = diagnostic_record_every;
res.exact_sync_diagnostic_every = exact_sync_diagnostic_every;
res.rolling_window = rolling_window;
res.avg_start_policy = avg_config.policy;
res.iter = h;
res.price_updates = h;

    function hash_state = updateHash(hash_state,value)
        hash_state = mod(hash_state * 16777619 + double(value),4294967291);
    end

    function [delay_hash,update_hash] = computeScheduleHashes(seed,h0_local,N_local,horizon)
        preview_stream = RandStream('mt19937ar','Seed',seed);
        preview_delay = randi(preview_stream,[0,h0_local],N_local,1);
        preview_count = zeros(N_local,1);
        delay_state = 2166136261;
        update_state = 2166136261;
        for jj = 1:N_local
            delay_state = updateHash(delay_state,preview_delay(jj) + 1);
        end
        for hh = 1:horizon
            for jj = 1:N_local
                if preview_count(jj) < preview_delay(jj)
                    preview_count(jj) = preview_count(jj) + 1;
                else
                    preview_count(jj) = 0;
                    preview_delay(jj) = randi(preview_stream,[0,h0_local]);
                    delay_state = updateHash(delay_state,preview_delay(jj) + 1);
                    update_state = updateHash(update_state,hh);
                    update_state = updateHash(update_state,jj);
                end
            end
        end
        delay_hash = sprintf('%.0f',delay_state);
        update_hash = sprintf('%.0f',update_state);
    end

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
