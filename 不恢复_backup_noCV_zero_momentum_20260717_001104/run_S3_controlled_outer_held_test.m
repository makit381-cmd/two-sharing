function summary = run_S3_controlled_outer_held_test(controlled_outer_price_updates)
% S3 held 外层受控回归：只验证五社区事件/计数/随机流，不调用昂贵 local solver。
% 这是语义回归，不是完整 S3 收敛或 PTDF 优化结果。

if nargin < 1 || isempty(controlled_outer_price_updates)
    controlled_outer_price_updates = 20;
end
validateattributes(controlled_outer_price_updates,{'numeric'}, ...
    {'scalar','real','finite','positive','integer'});
P = load('params.mat');
selected = [1,3,49,71,123];
formal_start = P.outer_avg_start_price_update;
formal_parameters = struct('outer_avg_start_price_update',formal_start, ...
    'outer_avg_start_iter',P.outer_avg_start_iter, ...
    'h0_max',P.h0_max,'k0_max',P.k0_max(selected));
test_overrides = struct('outer_avg_start_price_update',NaN, ...
    'local_solver','stub_only','label','CONTROLLED REGRESSION ONLY; NOT FORMAL S3');

case_off = simulateHeldCase(false);
case_on = simulateHeldCase(true);
schedule_pass = isequal(case_off.active_mask,case_on.active_mask) && ...
    isequal(case_off.delay_history,case_on.delay_history) && ...
    sameWithNaN(case_off.qre_noise_seed_history,case_on.qre_noise_seed_history);
held_pass = case_on.held_preservation_pass;
count_pass = case_on.outer_wall_iter == controlled_outer_price_updates && ...
    case_on.outer_price_update_count == controlled_outer_price_updates && ...
    all(case_on.local_call_count >= 0) && ...
    case_on.outer_avg_cnt == max(controlled_outer_price_updates - formal_start + 1,0);
finite_pass = all(isfinite(case_on.pb_cv_kW)) && all(isfinite(case_on.outer_price_step));
summary = struct();
summary.selected_communities = selected;
summary.controlled_outer_price_updates = controlled_outer_price_updates;
summary.formal_parameters = formal_parameters;
summary.effective_parameters = formal_parameters;
summary.test_overrides = test_overrides;
summary.audit_switch_lifecycle_pass = schedule_pass;
summary.held_state_pass = held_pass;
summary.count_semantics_pass = count_pass;
summary.finite_pass = finite_pass;
summary.success_gate = schedule_pass && held_pass && count_pass && finite_pass;
summary.audit_off = case_off;
summary.audit_on = case_on;
summary.table = table((1:controlled_outer_price_updates)', ...
    case_on.active_count,case_on.held_count,case_on.local_call_total, ...
    case_on.outer_avg_cnt_history,case_on.pb_cv_kW,case_on.line_cv_kW, ...
    'VariableNames',{'update','active_communities','held_communities', ...
    'local_calls','outer_avg_cnt','PB_CV_kW','line_CV_kW'});
summary.semantics = 'S3 held controlled regression only; each wall slot performs one price update using active plus held states; no full local solver claim';
save(sprintf('S3_controlled_outer_held_%03d.mat',controlled_outer_price_updates),'summary');
fprintf('S3 controlled held regression (%d updates): %s\n', ...
    controlled_outer_price_updates,ternary(summary.success_gate,'PASS','FAIL'));

    function out = simulateHeldCase(audit_enabled)
        nsel = numel(selected);
        state_pesc = selected(:) .* 1e-3; % 非零 warm-start 状态
        pi0_state = P.init_pi_0(selected);
        delay_remaining = zeros(nsel,1);
        local_call_count = zeros(nsel,1);
        outer_avg_cnt = 0;
        outer_wall_iter = 0;
        outer_price_update_count = 0;
        sum_pesc = zeros(nsel,1);
        active_mask = false(controlled_outer_price_updates,nsel);
        delay_history = nan(controlled_outer_price_updates,nsel);
        qre_noise_seed_history = nan(controlled_outer_price_updates,nsel);
        active_count = zeros(controlled_outer_price_updates,1);
        held_count = zeros(controlled_outer_price_updates,1);
        local_call_total = zeros(controlled_outer_price_updates,1);
        outer_avg_cnt_history = zeros(controlled_outer_price_updates,1);
        pb_cv_kW = zeros(controlled_outer_price_updates,1);
        line_cv_kW = zeros(controlled_outer_price_updates,1); % stub 不构造 PTDF 线路
        outer_price_step = zeros(controlled_outer_price_updates,1);
        held_preservation_pass = true;
        audit_count = zeros(controlled_outer_price_updates,1);
        for wall = 1:controlled_outer_price_updates
            outer_wall_iter = outer_wall_iter + 1;
            previous_state = state_pesc;
            previous_pi0 = pi0_state;
            for q = 1:nsel
                delay_history(wall,q) = delay_remaining(q);
                if delay_remaining(q) > 0
                    delay_remaining(q) = delay_remaining(q) - 1;
                else
                    active_mask(wall,q) = true;
                    local_call_count(q) = local_call_count(q) + 1;
                    noise_seed = deriveDeterministicSeed(P.rng_seeds.qre_noise, ...
                        selected(q),local_call_count(q),2);
                    qre_noise_seed_history(wall,q) = noise_seed;
                    state_pesc(q) = previous_state(q) + ...
                        1e-4 * (1 + mod(double(noise_seed),1000)/1000) + ...
                        1e-5 * previous_pi0(q);
                    delay_seed = deriveDeterministicSeed(P.rng_seeds.outer_delay, ...
                        selected(q),local_call_count(q),4);
                    delay_stream = RandStream('mt19937ar','Seed',delay_seed);
                    delay_remaining(q) = randi(delay_stream,[0,P.k0_max(selected(q))]);
                    if audit_enabled
                        audit_count(wall) = audit_count(wall) + 1;
                        % audit 开关只记录，不参与 delay/QRE noise 生成。
                        deriveDeterministicSeed(P.rng_seeds.qre_audit, ...
                            selected(q),local_call_count(q),3);
                    end
                end
            end
            inactive = ~active_mask(wall,:); %#ok<NASGU>
            if any(~active_mask(wall,:))
                held_preservation_pass = held_preservation_pass && ...
                    all(state_pesc(~active_mask(wall,:)) == previous_state(~active_mask(wall,:)));
            end
            % held 全系统状态每个 wall slot 都参与一次受控基础价格更新。
            outer_price_update_count = outer_price_update_count + 1;
            outer_price_step(wall) = norm(pi0_state - (pi0_state - P.alpha_PB .* state_pesc),inf);
            pi0_state = pi0_state - P.alpha_PB .* state_pesc;
            if outer_price_update_count >= formal_start
                outer_avg_cnt = outer_avg_cnt + 1;
                sum_pesc = sum_pesc + state_pesc;
            end
            active_count(wall) = sum(active_mask(wall,:));
            held_count(wall) = nsel - active_count(wall);
            local_call_total(wall) = sum(local_call_count);
            outer_avg_cnt_history(wall) = outer_avg_cnt;
            pb_cv_kW(wall) = abs(sum(state_pesc));
        end
        out = struct('audit_enabled',audit_enabled, ...
            'active_mask',active_mask,'delay_history',delay_history, ...
            'qre_noise_seed_history',qre_noise_seed_history, ...
            'active_count',active_count,'held_count',held_count, ...
            'local_call_total',local_call_total,'local_call_count',local_call_count, ...
            'outer_avg_cnt_history',outer_avg_cnt_history,'outer_avg_cnt',outer_avg_cnt, ...
            'outer_wall_iter',outer_wall_iter,'outer_price_update_count',outer_price_update_count, ...
            'pb_cv_kW',pb_cv_kW,'line_cv_kW',line_cv_kW, ...
            'outer_price_step',outer_price_step,'audit_count',audit_count, ...
            'held_preservation_pass',held_preservation_pass,'warm_start_pesc',selected(:).*1e-3, ...
            'final_pesc',state_pesc,'final_pi_0',pi0_state);
    end

    function out = ternary(condition,yes_value,no_value)
        if condition, out = yes_value; else, out = no_value; end
    end
    function out = sameWithNaN(x,y)
        out = isequal(size(x),size(y)) && all((x(:) == y(:)) | ...
            (isnan(x(:)) & isnan(y(:))));
    end
end
