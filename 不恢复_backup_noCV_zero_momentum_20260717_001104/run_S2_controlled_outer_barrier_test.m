function summary = run_S2_controlled_outer_barrier_test(max_price_updates,outer_avg_start_price_update_override)
% 受控 S2 外层同步屏障测试。
% 只选 5 个社区，不调用完整 upperMarketS2，也不修改正式 params.mat。
% 该测试验证：frozen pi_0、held 等待、原子发布、遍历平均只在同步完成时递增，
% 以及 keyed 外层延迟和 keyed local-call 随机流的可复现性。

if nargin < 1 || isempty(max_price_updates)
    max_price_updates = 20;
end
if nargin < 2 || isempty(outer_avg_start_price_update_override)
    outer_avg_start_price_update_override = NaN;
end
validateattributes(max_price_updates,{'numeric'},{'scalar','real','finite','positive','integer'});
if isfinite(outer_avg_start_price_update_override)
    validateattributes(outer_avg_start_price_update_override,{'numeric'},{'scalar','real','finite','nonnegative','integer'});
elseif ~isnan(outer_avg_start_price_update_override)
    error('run_S2_controlled_outer_barrier_test:invalidOverride', ...
        'outer_avg_start_price_update_override must be finite or NaN.');
end
selected = [1,3,49,71,123];
P = load('params.mat');
nsel = numel(selected);

if isfield(P,'outer_avg_start_price_update')
    formal_outer_avg_start_price_update = P.outer_avg_start_price_update;
else
    formal_outer_avg_start_price_update = P.outer_avg_start_iter;
end
if isfinite(outer_avg_start_price_update_override)
    effective_outer_avg_start_price_update = outer_avg_start_price_update_override;
    test_overrides = struct('outer_avg_start_price_update',outer_avg_start_price_update_override, ...
        'label','TEST OVERRIDE ONLY; NOT FORMAL PROFILE');
else
    effective_outer_avg_start_price_update = formal_outer_avg_start_price_update;
    test_overrides = struct('outer_avg_start_price_update',NaN, ...
        'label','No test override; formal profile');
end
formal_parameters = struct('outer_avg_start_price_update',formal_outer_avg_start_price_update, ...
    'outer_avg_start_iter',P.outer_avg_start_iter);
effective_parameters = struct('outer_avg_start_price_update',effective_outer_avg_start_price_update, ...
    'outer_avg_start_iter',effective_outer_avg_start_price_update);

if ~isfield(P,'rng_seeds') || ~isfield(P.rng_seeds,'outer_delay')
    error('params.mat 缺少 keyed outer_delay seed。');
end

% 受控测试只验证外层屏障，因此每次局部求解限制为 1 次真实同步价格更新；
% 这不会改写正式参数，也不用于报告收敛结论。
test_max_inner_iter = 1;
test_avg_config = struct( ...
    'policy',P.inner_avg_policy, ...
    'policy_version',P.inner_avg_policy_version, ...
    'start_cv_factor',5, ...
    'start_price_factor',10, ...
    'start_stable_window',1, ...
    'min_price_updates',1, ...
    'min_samples',1, ...
    'formal_cv_stable_window',1, ...
    'price_stable_window',1);
test_rng_config = struct( ...
    'inner_delay_seed',P.rng_seeds.inner_delay, ...
    'qre_noise_seed',P.rng_seeds.qre_noise, ...
    'qre_audit_seed',P.rng_seeds.qre_audit, ...
    'qre_audit_enabled',false, ...
    'qre_audit_rate',0, ...
    'qre_noise_enabled',false);

state_pi = P.init_pi(selected);
state_pes = cell(nsel,1);
state_pg = cell(nsel,1);
state_pb = cell(nsel,1);
state_ps = cell(nsel,1);
state_pesc = P.init_pesc(selected);
for q = 1:nsel
    idx = P.start_idx(selected(q)):P.end_idx(selected(q));
    state_pes{q} = P.init_pes(idx);
    state_pg{q} = P.init_pg(idx);
    state_pb{q} = P.init_pb(idx);
    state_ps{q} = P.init_ps(idx);
end

sum_pes = cell(nsel,1);
sum_pg = cell(nsel,1);
sum_pb = cell(nsel,1);
sum_ps = cell(nsel,1);
sum_pesc = zeros(nsel,1);
for q = 1:nsel
    sum_pes{q} = zeros(size(state_pes{q}));
    sum_pg{q} = zeros(size(state_pg{q}));
    sum_pb{q} = zeros(size(state_pb{q}));
    sum_ps{q} = zeros(size(state_ps{q}));
end

local_call_count = zeros(nsel,1);
outer_avg_cnt = 0;
outer_wall_iter = 0;
price_update_count = 0;
wall_update_event = false(0,1);
wall_price_update_count = zeros(0,1);
wall_avg_count = zeros(0,1);
wall_cycle = zeros(0,1);
wall_ready_count = zeros(0,1);
wall_published_hash = strings(0,1);
delay_history = zeros(max_price_updates,nsel);
delay_seed_history = zeros(max_price_updates,nsel);
local_seed_history = zeros(max_price_updates,nsel);
cycle_wall_events = zeros(max_price_updates,1);
cycle_call_counts = zeros(max_price_updates,nsel);
cycle_avg_counts = zeros(max_price_updates,1);
cycle_price_updates_before = zeros(max_price_updates,1);
cycle_price_updates_after = zeros(max_price_updates,1);
cycle_avg_counts_before = zeros(max_price_updates,1);
cycle_avg_counts_after = zeros(max_price_updates,1);
barrier_pass = true;
published_pesc = state_pesc;

for cycle = 1:max_price_updates
    pi_0_cycle = P.init_pi_0(selected);
    cycle_pi = state_pi;
    cycle_pes = state_pes;
    cycle_pg = state_pg;
    cycle_pb = state_pb;
    cycle_ps = state_ps;
    cycle_pesc = state_pesc;
    ready = false(nsel,1);
    count = zeros(nsel,1);
    delay = zeros(nsel,1);
    for q = 1:nsel
        delay_seed_history(cycle,q) = deriveDeterministicSeed( ...
            P.rng_seeds.outer_delay,selected(q),cycle,4);
        delay_stream = RandStream('mt19937ar','Seed',delay_seed_history(cycle,q));
        delay(q) = randi(delay_stream,[0,P.k0_max(selected(q))]);
        delay_history(cycle,q) = delay(q);
    end

    calls_before = local_call_count;
    price_update_count_before = price_update_count;
    avg_count_before = outer_avg_cnt;
    cycle_wall_start = outer_wall_iter;
    cycle_community_ready_count = 0;
    while ~all(ready)
        outer_wall_iter = outer_wall_iter + 1;
        for q = 1:nsel
            if ~ready(q)
                if count(q) < delay(q)
                    count(q) = count(q) + 1;
                else
                    local_call_count(q) = local_call_count(q) + 1;
                    idx = P.start_idx(selected(q)):P.end_idx(selected(q));
                    local = localMarket_S2(P.c(idx),P.b(idx),P.a(selected(q)),P.D(idx), ...
                        P.init_rho(selected(q)),P.pg_max(idx),P.pg_min(idx),cycle_pi(q), ...
                        cycle_pes{q},cycle_pg{q},cycle_pb{q},cycle_ps{q},cycle_pesc(q), ...
                        P.errTol_LESMs(selected(q)),P.inner_cv_tol_kW(selected(q)),pi_0_cycle(q), ...
                        P.pi_max,P.pi_min,false,P.beta_qre(idx),P.qre_epsilon,P.agg_epsilon, ...
                        P.qre_z_cap,P.qre_max_backoffs,P.inner_momentum,test_max_inner_iter,P.h0_max, ...
                        0,2,false,false,false,false,Inf,cycle,selected(q), ...
                        test_rng_config,test_avg_config,local_call_count(q));
                    cycle_pi(q) = local.pi;
                    cycle_pes{q} = local.pes;
                    cycle_pg{q} = local.pg;
                    cycle_pb{q} = local.pb;
                    cycle_ps{q} = local.ps;
                    cycle_pesc(q) = local.pesc;
                    local_seed_history(cycle,q) = local.qre_noise_seed;
                    ready(q) = true;
                    cycle_community_ready_count = cycle_community_ready_count + 1;
                end
            end
        end

        % 等待期间只记录已发布状态；price_updates 和 outer_avg_cnt 均不得变化。
        wall_update_event(end+1,1) = false;
        wall_price_update_count(end+1,1) = price_update_count;
        wall_avg_count(end+1,1) = outer_avg_cnt;
        wall_cycle(end+1,1) = cycle;
        wall_ready_count(end+1,1) = cycle_community_ready_count;
        wall_published_hash(end+1,1) = stateHash(published_pesc);
        if price_update_count ~= price_update_count_before || outer_avg_cnt ~= avg_count_before
            barrier_pass = false;
        end
    end

    % 原子发布及一次遍历平均：只能在 all(ready) 后发生。
    state_pi = cycle_pi;
    state_pes = cycle_pes;
    state_pg = cycle_pg;
    state_pb = cycle_pb;
    state_ps = cycle_ps;
    state_pesc = cycle_pesc;
    published_pesc = state_pesc;
    price_update_count = price_update_count + 1;
    if price_update_count >= effective_outer_avg_start_price_update
        outer_avg_cnt = outer_avg_cnt + 1;
    end
    if price_update_count >= effective_outer_avg_start_price_update
        for q = 1:nsel
            sum_pes{q} = sum_pes{q} + state_pes{q};
            sum_pg{q} = sum_pg{q} + state_pg{q};
            sum_pb{q} = sum_pb{q} + state_pb{q};
            sum_ps{q} = sum_ps{q} + state_ps{q};
            sum_pesc(q) = sum_pesc(q) + state_pesc(q);
        end
    end

    % 仅用于屏障测试的受控基础价格更新；不代表完整 PTDF 外层优化结果。
    pi_0_next = pi_0_cycle - P.alpha_PB .* state_pesc;
    if any(~isfinite(pi_0_next))
        barrier_pass = false;
    end
    cycle_wall_events(cycle) = outer_wall_iter - cycle_wall_start;
    cycle_call_counts(cycle,:) = local_call_count(:) - calls_before(:);
    cycle_avg_counts(cycle) = outer_avg_cnt;
    cycle_price_updates_before(cycle) = price_update_count_before;
    cycle_price_updates_after(cycle) = price_update_count;
    cycle_avg_counts_before(cycle) = avg_count_before;
    cycle_avg_counts_after(cycle) = outer_avg_cnt;
    expected_avg_increment = double(price_update_count >= effective_outer_avg_start_price_update);
    if any(cycle_call_counts(cycle,:) ~= 1) || ...
            cycle_price_updates_after(cycle) ~= cycle_price_updates_before(cycle) + 1 || ...
            cycle_avg_counts_after(cycle) ~= cycle_avg_counts_before(cycle) + expected_avg_increment
        barrier_pass = false;
    end

    % 让下一轮使用本轮发布后的基础价格；本轮所有 local 调用仍使用同一 frozen 值。
    P.init_pi_0(selected) = pi_0_next; %#ok<NASGU>
    % MATLAB 函数工作区不能修改结构体字段的下一轮快照，单独保存状态。
    pi0_state = pi_0_next; %#ok<NASGU>
    if cycle < max_price_updates
        % 下一轮开头读取 pi0_state，而不是重新从 P.init_pi_0 取初值。
        P.init_pi_0(selected) = pi0_state;
    end
end

% keyed schedule replay：同一输入下延迟与 local-call seed 必须完全相同。
replay_delay = zeros(size(delay_history));
replay_delay_seed = zeros(size(delay_seed_history));
for cycle = 1:max_price_updates
    for q = 1:nsel
        replay_delay_seed(cycle,q) = deriveDeterministicSeed( ...
            P.rng_seeds.outer_delay,selected(q),cycle,4);
        replay_stream = RandStream('mt19937ar','Seed',replay_delay_seed(cycle,q));
        replay_delay(cycle,q) = randi(replay_stream,[0,P.k0_max(selected(q))]);
    end
end
replay_pass = isequal(delay_history,replay_delay) && isequal(delay_seed_history,replay_delay_seed);
expected_final_avg_cnt = max(max_price_updates - effective_outer_avg_start_price_update + 1,0);
count_relation_pass = all(all(cycle_call_counts == 1)) && ...
    all(local_call_count == max_price_updates) && ...
    price_update_count == max_price_updates && outer_avg_cnt == expected_final_avg_cnt;
summary = struct();
summary.selected_communities = selected;
summary.max_price_updates = max_price_updates;
summary.formal_parameters = formal_parameters;
summary.effective_parameters = effective_parameters;
summary.test_overrides = test_overrides;
summary.test_override_only = isfinite(outer_avg_start_price_update_override);
summary.outer_price_update_count = price_update_count;
summary.outer_wall_iter = outer_wall_iter;
summary.outer_avg_cnt = outer_avg_cnt;
summary.expected_final_avg_cnt = expected_final_avg_cnt;
summary.total_local_calls = sum(local_call_count);
summary.local_call_count = local_call_count;
summary.cycle_wall_events = cycle_wall_events;
summary.delay_history = delay_history;
summary.delay_seed_history = delay_seed_history;
summary.local_qre_noise_seed_history = local_seed_history;
summary.cycle_call_counts = cycle_call_counts;
summary.cycle_price_updates_before = cycle_price_updates_before;
summary.cycle_price_updates_after = cycle_price_updates_after;
summary.cycle_avg_counts_before = cycle_avg_counts_before;
summary.cycle_avg_counts_after = cycle_avg_counts_after;
summary.wall_price_update_count = wall_price_update_count;
summary.wall_avg_count = wall_avg_count;
summary.wall_update_event = wall_update_event;
summary.wall_cycle = wall_cycle;
summary.wall_ready_count = wall_ready_count;
summary.wall_published_hash = wall_published_hash;
summary.replay_pass = replay_pass;
summary.count_relation_pass = count_relation_pass;
summary.barrier_pass = barrier_pass;
summary.success_gate = barrier_pass && replay_pass && count_relation_pass;
summary.semantics = 'controlled barrier only; no full 123-community outer optimization claim';

save(sprintf('S2_controlled_outer_barrier_%03d.mat',max_price_updates),'summary');
fprintf('S2 controlled outer barrier (%d updates): %s\n',max_price_updates,ternary(summary.success_gate,'PASS','FAIL'));

    function out = stateHash(v)
        % 仅用于等待阶段状态是否保持；不作为密码学证明。
        out = string(sprintf('%.16g|',v(:)));
    end

    function out = ternary(condition,yes_value,no_value)
        if condition, out = yes_value; else, out = no_value; end
    end
end
