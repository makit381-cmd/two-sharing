function res = localMarket_S2(c,b,a,D,init_rho,pg_max,pg_min,init_pi, ...
    init_pes,init_pg,init_pb,init_ps,init_pesc,errTol,inner_cv_tol_kW, ...
    pi_0,pi_max,pi_min,record_trace,beta_qre,~,~,qre_z_cap,~, ...
    inner_momentum,max_inner_iter,h0,min_inner_iter,inner_avg_start_iter,inner_cv_stop_enabled, ...
    qre_certificate_enabled,agg_certificate_enabled,progress_print_enabled, ...
    progress_print_inner_every,outer_iter,community_id,rng_config,avg_config,local_call_count)

% S2 内层：与 S3 相同的 raw Agg/QRE/动态记录写法；不调用公共算法函数。
% 瞬时 raw 残差只用于 lambda 更新；正式 CV/RE 使用遍历平均输出。

N = length(b);

validateattributes(init_rho, {'numeric'}, {'scalar','real','finite','positive'});
validateattributes(inner_cv_tol_kW, {'numeric'}, {'scalar','real','finite','nonnegative'});
validateattributes(h0, {'numeric'}, {'scalar','real','finite','nonnegative','integer'});
validateattributes(min_inner_iter, {'numeric'}, {'scalar','real','finite','integer','nonnegative'});
validateattributes(inner_avg_start_iter, {'numeric'}, {'scalar','real','finite','positive','integer'});
validateattributes(inner_cv_stop_enabled, {'logical','numeric'}, {'scalar','real','finite'});
if nargin < 35
    progress_print_enabled = false;
    progress_print_inner_every = inf;
    outer_iter = NaN;
    community_id = NaN;
end
if nargin < 36 || isempty(rng_config), rng_config = struct(); end
if nargin < 37 || isempty(avg_config), avg_config = struct(); end
if nargin < 38 || isempty(local_call_count), local_call_count = 1; end
if ~isfield(rng_config,'inner_delay_seed'), rng_config.inner_delay_seed = 31001; end
if ~isfield(rng_config,'qre_noise_seed'), rng_config.qre_noise_seed = 31002; end
if ~isfield(avg_config,'policy'), avg_config.policy = 'A2_dynamic_stable_start'; end
if ~isfield(avg_config,'policy_version'), avg_config.policy_version = 'A2_dynamic_stable_start_v1'; end
if ~isfield(avg_config,'start_cv_factor'), avg_config.start_cv_factor = 5; end
if ~isfield(avg_config,'start_price_factor'), avg_config.start_price_factor = 10; end
if ~isfield(avg_config,'start_stable_window'), avg_config.start_stable_window = 20; end
if ~isfield(avg_config,'min_price_updates'), avg_config.min_price_updates = 100; end
if ~isfield(avg_config,'min_samples'), avg_config.min_samples = 200; end
if ~isfield(avg_config,'formal_cv_stable_window'), avg_config.formal_cv_stable_window = 5; end
if ~isfield(avg_config,'price_stable_window'), avg_config.price_stable_window = 20; end
if ~isfield(rng_config,'qre_noise_enabled'), rng_config.qre_noise_enabled = true; end
validateattributes(local_call_count,{'numeric'},{'scalar','real','finite','positive','integer'});
validateattributes(avg_config.start_stable_window,{'numeric'},{'scalar','real','finite','positive','integer'});
validateattributes(avg_config.formal_cv_stable_window,{'numeric'},{'scalar','real','finite','positive','integer'});
validateattributes(avg_config.price_stable_window,{'numeric'},{'scalar','real','finite','positive','integer'});
if qre_certificate_enabled || agg_certificate_enabled
    error('localMarket_S2 is currently the unchecked_qre profile; enable certificates only after restoring the certified implementation.');
end

rho = init_rho;
alpha = rho * a;
community_key = 0;
if isfinite(community_id), community_key = max(0,floor(community_id)); end
inner_delay_seed = deriveDeterministicSeed(rng_config.inner_delay_seed, ...
    community_key,local_call_count,1);
qre_noise_seed = deriveDeterministicSeed(rng_config.qre_noise_seed, ...
    community_key,local_call_count,2);
delay_stream = RandStream('mt19937ar','Seed',inner_delay_seed);
qre_stream = RandStream('mt19937ar','Seed',qre_noise_seed);
delay = randi(delay_stream,[0,h0],N,1); % 与 S3 相同的 keyed 延迟机制。
count = zeros(N,1);
ready = false(N,1);
wall_iter = 0;
price_update_count = 0;
eu_response_count = 0;

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
averaging_started = false;
avg_start_stable_count = 0;
inner_avg_start_wall_iter = NaN;
inner_avg_start_price_update = NaN;
recent_lambda_steps = nan(avg_config.price_stable_window,1);
recent_lambda_step_count = 0;
recent_lambda_step_position = 0;
recent_formal_cv_kW = nan(avg_config.formal_cv_stable_window,1);
recent_formal_cv_count = 0;
recent_formal_cv_position = 0;
formal_cv_stable_pass = false;
formal_cv_kW = NaN;
stable_price_pass = false;
bar_pes = init_pes;
bar_pg = init_pg;
bar_pb = init_pb;
bar_ps = init_ps;
bar_pesc_raw = init_pesc;
tol_F = 1e-8;

% 仅选中的社区保存 RE/CV 轨迹；这些是一维数组，按需扩容。
if record_trace
    trace_capacity = min(max_inner_iter + 1, 128);
    trace_h = nan(trace_capacity,1);
    trace_violation = nan(trace_capacity,1);
    trace_obj = nan(trace_capacity,1);
    trace_lambda = nan(trace_capacity,1);
    trace_lambda_used = nan(trace_capacity,1);
    trace_lambda_step = nan(trace_capacity,1);
    trace_price_update_count = nan(trace_capacity,1);
    trace_update_event = false(trace_capacity,1);
    trace_avg_cnt = zeros(trace_capacity,1);
    trace_h(1) = 0;
    trace_violation(1) = abs(init_pesc - sum(init_pes));
    trace_obj(1) = communityPrimalObj(init_pg, init_pb, init_ps, init_pes, init_pesc);
    trace_lambda(1) = init_pi;
    trace_lambda_used(1) = init_pi;
    trace_lambda_step(1) = 0;
    trace_price_update_count(1) = 0;
    trace_avg_cnt(1) = 0;
    trace_count = 1;
end

res.success = 1;
res.stop_reason = '';
h = 1;

if progress_print_enabled
    fprintf('[S2 inner start] k=%d, community=%d, EU=%d, CV limit=%.3f kW\n', ...
        outer_iter, community_id, N, inner_cv_tol_kW);
end

while true
    wall_iter = wall_iter + 1;
    lambda_cycle = lambda;
    lambda_used_now = lambda_cycle;
    for j = 1:N
        if ~ready(j)
        if count(j) < delay(j)
            count(j) = count(j) + 1;
        else
            % 未认证 QRE 的解析响应直接内联：保留原有两次 randn 的顺序、
            % 投影和买家/卖家/内部三分支，避免每次 EU 更新的两层函数调用。
            if rng_config.qre_noise_enabled
                z_pg = max(min(randn(qre_stream,1), qre_z_cap), -qre_z_cap);
                z_pes = max(min(randn(qre_stream,1), qre_z_cap), -qre_z_cap);
            else
                z_pg = 0;
                z_pes = 0;
            end
            c_j = c(j);
            b_j = b(j);
            D_j = D(j);
            pg_min_j = pg_min(j);
            pg_max_j = pg_max(j);
            beta_j = beta_qre(j);

            z_middle = (sqrt(c_j) * z_pg + sqrt(a) * z_pes) / sqrt(c_j + a);
            pg_buyer = min(max((pi_max - b_j) / c_j + ...
                sqrt(1 / (beta_j * c_j)) * z_pg, pg_min_j), pg_max_j);
            pes_buyer = (lambda_used_now - pi_max) / a + ...
                sqrt(1 / (beta_j * a)) * z_pes;
            F_buyer = D_j + pes_buyer - pg_buyer;

            pg_seller = min(max((pi_min - b_j) / c_j + ...
                sqrt(1 / (beta_j * c_j)) * z_pg, pg_min_j), pg_max_j);
            pes_seller = (lambda_used_now - pi_min) / a + ...
                sqrt(1 / (beta_j * a)) * z_pes;
            F_seller = D_j + pes_seller - pg_seller;

            if F_buyer > tol_F
                fast_pg(j) = pg_buyer;
                fast_pes(j) = pes_buyer;
                fast_pb(j) = F_buyer;
                fast_ps(j) = 0;
            elseif F_seller < -tol_F
                fast_pg(j) = pg_seller;
                fast_pes(j) = pes_seller;
                fast_pb(j) = 0;
                fast_ps(j) = -F_seller;
            else
                pes_unclipped = (lambda_used_now - b_j - c_j * D_j) / (c_j + a) + ...
                    sqrt(1 / (beta_j * (c_j + a))) * z_middle;
                fast_pes(j) = min(max(pes_unclipped, pg_min_j - D_j), pg_max_j - D_j);
                fast_pg(j) = D_j + fast_pes(j);
                fast_pb(j) = 0;
                fast_ps(j) = 0;
            end
            ready(j) = true;
        end
        end
    end
    % 同步等待：仅在本轮所有 EU 均用同一 lambda 响应后更新 lambda。
    if ~all(ready)
        if record_trace
            trace_count = trace_count + 1;
            if trace_count > trace_capacity
                trace_capacity = growHistoryCapacity(trace_capacity, trace_count, Inf);
                trace_h(trace_capacity,1) = NaN;
                trace_violation(trace_capacity,1) = NaN;
                trace_obj(trace_capacity,1) = NaN;
                trace_lambda(trace_capacity,1) = NaN;
                trace_lambda_used(trace_capacity,1) = NaN;
                trace_lambda_step(trace_capacity,1) = NaN;
                trace_price_update_count(trace_capacity,1) = NaN;
                trace_update_event(trace_capacity,1) = false;
                trace_avg_cnt(trace_capacity,1) = 0;
            end
            trace_h(trace_count) = wall_iter;
            trace_violation(trace_count) = abs(bar_pesc_raw - sum(bar_pes));
            trace_obj(trace_count) = communityPrimalObj(bar_pg, bar_pb, bar_ps, bar_pes, bar_pesc_raw);
            trace_lambda(trace_count) = lambda_cycle;
            trace_lambda_used(trace_count) = lambda_cycle;
            trace_lambda_step(trace_count) = 0;
            trace_price_update_count(trace_count) = price_update_count;
            trace_avg_cnt(trace_count) = avg_cnt;
        end
        continue;
    end

    pesc_now = (pi_0 - lambda_used_now) / a;

    price_update_count = price_update_count + 1;
    h = price_update_count;
    lambda_aux = lambda_cycle + alpha * (pesc_now - sum(fast_pes));
    lambda_next = lambda_aux + inner_momentum * (lambda_aux - lambda_aux_prev);
    lambda_step = abs(lambda_next - lambda_cycle);
    sync_cv_kW = abs(pesc_now - sum(fast_pes));
    averaging_started_before = averaging_started;
    start_cv_pass = sync_cv_kW <= avg_config.start_cv_factor * inner_cv_tol_kW;
    start_price_pass = isfinite(errTol) && ...
        lambda_step <= avg_config.start_price_factor * errTol;
    if start_cv_pass && start_price_pass
        avg_start_stable_count = avg_start_stable_count + 1;
    else
        avg_start_stable_count = 0;
    end
    if ~averaging_started && ...
            price_update_count >= avg_config.min_price_updates && ...
            avg_start_stable_count >= avg_config.start_stable_window
        averaging_started = true;
        inner_avg_start_wall_iter = wall_iter;
        inner_avg_start_price_update = price_update_count;
    end

    if averaging_started_before
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
    inner_cv_bar_now_kW = abs(bar_pesc_raw - sum(bar_pes));
    recent_lambda_step_count = min(recent_lambda_step_count + 1, numel(recent_lambda_steps));
    recent_lambda_step_position = mod(recent_lambda_step_position, numel(recent_lambda_steps)) + 1;
    recent_lambda_steps(recent_lambda_step_position) = lambda_step;
    stable_price_pass = recent_lambda_step_count >= numel(recent_lambda_steps) && ...
        all(recent_lambda_steps <= errTol);
    if averaging_started
        recent_formal_cv_count = min(recent_formal_cv_count + 1, numel(recent_formal_cv_kW));
        recent_formal_cv_position = mod(recent_formal_cv_position, numel(recent_formal_cv_kW)) + 1;
        recent_formal_cv_kW(recent_formal_cv_position) = inner_cv_bar_now_kW;
    end
    formal_cv_kW = inner_cv_bar_now_kW;
    formal_cv_stable_pass = averaging_started && ...
        avg_cnt >= avg_config.min_samples && ...
        recent_formal_cv_count >= numel(recent_formal_cv_kW) && ...
        all(recent_formal_cv_kW <= inner_cv_tol_kW);
    eu_response_count = eu_response_count + N;

    if record_trace
        trace_count = trace_count + 1;
        if trace_count > trace_capacity
            trace_capacity = growHistoryCapacity(trace_capacity, trace_count, Inf);
            trace_h(trace_capacity,1) = NaN;
            trace_violation(trace_capacity,1) = NaN;
            trace_obj(trace_capacity,1) = NaN;
            trace_lambda(trace_capacity,1) = NaN;
            trace_lambda_used(trace_capacity,1) = NaN;
            trace_lambda_step(trace_capacity,1) = NaN;
            trace_price_update_count(trace_capacity,1) = NaN;
            trace_update_event(trace_capacity,1) = false;
            trace_avg_cnt(trace_capacity,1) = 0;
        end
        trace_h(trace_count) = wall_iter;
        trace_violation(trace_count) = inner_cv_bar_now_kW;
        trace_obj(trace_count) = communityPrimalObj(bar_pg, bar_pb, bar_ps, bar_pes, bar_pesc_raw);
        trace_lambda(trace_count) = lambda_next;
        trace_lambda_used(trace_count) = lambda_cycle;
        trace_lambda_step(trace_count) = lambda_step;
        trace_price_update_count(trace_count) = price_update_count;
        trace_update_event(trace_count) = true;
        trace_avg_cnt(trace_count) = avg_cnt;
    end

    if progress_print_enabled && (h == 1 || mod(h, progress_print_inner_every) == 0)
        if h < inner_avg_start_iter
            cv_label = 'CVnow';
        else
            cv_label = 'CVavg';
        end
        fprintf('[S2 inner] k=%d, community=%d, h=%d, %s=%.3f kW, diagnostic limit=%.3f kW, lambdaStep=%.3e\n', ...
            outer_iter, community_id, h, cv_label, inner_cv_bar_now_kW, inner_cv_tol_kW, lambda_step);
    end

    lambda_aux_prev = lambda_aux;
    lambda = lambda_next;
    ready(:) = false;
    count(:) = 0;
    delay = randi(delay_stream,[0,h0],N,1);

    inner_cv_pass = (~logical(inner_cv_stop_enabled)) || ...
        (formal_cv_kW <= inner_cv_tol_kW);
    if price_update_count >= max(min_inner_iter,1) && ...
            averaging_started && avg_cnt >= avg_config.min_samples && ...
            stable_price_pass && formal_cv_stable_pass && inner_cv_pass
        if logical(inner_cv_stop_enabled)
            res.stop_reason = 'dynamic_average_sync_price_step_and_relative_cv';
        else
            res.stop_reason = 'price_step_only';
        end
        break;
    end
    if isfinite(max_inner_iter) && price_update_count >= max_inner_iter
        res.success = 0;
        res.stop_reason = 'max_inner_iter_unchecked_qre';
        break;
    end
    h = h + 1;
end

if progress_print_enabled
    fprintf('[S2 inner end] k=%d, community=%d, h=%d, CVavg=%.3f kW, reason=%s\n', ...
        outer_iter, community_id, h, inner_cv_bar_now_kW, res.stop_reason);
end

if record_trace
    res.trace_h = trace_h(1:trace_count);
    res.trace_violation = trace_violation(1:trace_count);
    res.trace_obj = trace_obj(1:trace_count);
    res.trace_lambda = trace_lambda(1:trace_count);
    res.trace_lambda_used = trace_lambda_used(1:trace_count);
    res.trace_lambda_step = trace_lambda_step(1:trace_count);
    res.trace_price_update_count = trace_price_update_count(1:trace_count);
    res.trace_update_event = trace_update_event(1:trace_count);
    res.trace_avg_cnt = trace_avg_cnt(1:trace_count);
else
    res.trace_h = [];
    res.trace_violation = [];
    res.trace_obj = [];
    res.trace_lambda = [];
    res.trace_lambda_used = [];
    res.trace_lambda_step = [];
    res.trace_price_update_count = [];
    res.trace_update_event = [];
    res.trace_avg_cnt = [];
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
res.inner_avg_policy = avg_config.policy;
res.inner_avg_policy_version = avg_config.policy_version;
res.averaging_started = averaging_started;
res.inner_avg_start_wall_iter = inner_avg_start_wall_iter;
res.inner_avg_start_price_update = inner_avg_start_price_update;
res.formal_cv_kW = formal_cv_kW;
res.formal_cv_stable_pass = formal_cv_stable_pass;
res.stable_price_pass = stable_price_pass;
res.price_update_count = price_update_count;
res.sync_cycle_count = price_update_count;
res.wall_clock_count = wall_iter;
res.eu_response_count = eu_response_count;
res.inner_delay_seed = inner_delay_seed;
res.qre_noise_seed = qre_noise_seed;
res.community_id = community_id;
res.local_call_count = local_call_count;
res.seed_derivation_policy = ...
    'deriveDeterministicSeed(base_seed,community_id,local_call_count,stream_tag)';
res.inner_cv_stop_enabled = logical(inner_cv_stop_enabled);
res.iter = h;
res.wall_iter = wall_iter;
res.price_updates = price_update_count;

    function obj = communityPrimalObj(pg_v, pb_v, ps_v, pes_v, pesc_v)
        obj = sum(0.5 .* c(:) .* pg_v(:).^2 ...
            + b(:) .* pg_v(:) ...
            + pi_max .* pb_v(:) ...
            - pi_min .* ps_v(:) ...
            + 0.5 .* a .* pes_v(:).^2) ...
            + 0.5 .* a .* pesc_v.^2 - pi_0 .* pesc_v;
    end
end
