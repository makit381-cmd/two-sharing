clc;clear;
tic;
load('params.mat');
params_snapshot = load('params.mat');

pai_il = PTDF;%具体的量化网络约束系数    123*122
finite_line_mask = isfinite(F_l);
% %% 准备
% tic;tok;
state.success = 1;
outer_stop_reason = 'running';
stable_outer_price_count = 0;
k = 1;
val_pi_PB(1,k) = init_pi_PB;
val_pi_l_pos(:,k) = init_pi_l_pos;
val_pi_l_neg(:,k) = init_pi_l_neg;
val_pi_0(:,k) = init_pi_0;
pi = init_pi;
pes = init_pes;
pesc = init_pesc;
pg = init_pg;
pb = init_pb;
ps = init_ps;
pi_0 = init_pi_0;%迭代过程中的不变的能源共享基础价格

res(num_LESMs) = struct('success',[],'pi',[],'pes',[],'pg',[], ...
    'pb',[],'ps',[],'pesc',[],'pesc_agg_raw',[],'pes_sum',[], ...
    'inner_cv_kW',[],'inner_avg_cnt',[],'inner_avg_start_iter',[], ...
    'inner_avg_policy',[],'averaging_started',[],'inner_avg_start_wall_iter',[], ...
    'inner_avg_policy_version',[], ...
    'inner_avg_start_price_update',[],'inner_avg_start_cv_factor',[], ...
    'inner_avg_start_price_factor',[],'inner_avg_start_stable_window',[], ...
    'inner_avg_min_price_updates',[],'inner_avg_min_samples',[], ...
    'inner_formal_cv_stable_window',[], ...
    'formal_average_available',[],'final_formal_cv_kW',[],'formal_cv_pass',[], ...
    'formal_cv_stable_pass',[], ...
    'final_held_cv_kW',[],'final_lambda_step',[],'price_update_count',[],'wall_iter',[], ...
    'inner_delay_seed',[],'qre_noise_seed',[],'qre_audit_seed',[], ...
    'qre_noise_sequence_hash',[],'community_id',[],'local_call_count',[], ...
    'seed_derivation_policy',[], ...
    'delay_sequence_hash',[],'update_sequence_hash',[], ...
    'executed_delay_sequence_hash',[],'executed_update_sequence_hash',[], ...
    'delay_hash_horizon',[],'warm_start_initialized',[], ...
    'initialization_semantics',[], ...
    'formal_output_type',[],'inner_cv_stop_enabled',[],'iter',[],'price_updates',[],'stop_reason',[], ...
    'qre_certificate_enabled',[],'qre_epsilon',[],'qre_gap_max',[],'qre_gap_mean',[], ...
    'qre_gamma_min',[],'qre_gamma_mean',[],'qre_backoff_total',[],'qre_backoff_max',[], ...
    'qre_fallback_count',[],'qre_branch_switch_count',[],'qre_branch_switch_rate',[], ...
    'qre_all_certificate_pass',[],'qre_response_count',[],'qre_mode',[], ...
    'qre_audit_enabled',[],'qre_audit_rate',[],'qre_audit_seed',[],'qre_audit_count',[], ...
    'qre_gap_p95',[],'qre_gap_exceed_count',[],'qre_gap_exceed_rate',[], ...
    'qre_audit_response_count',[],'agg_certificate_enabled',[], ...
    'agg_gap',[],'agg_epsilon',[],'agg_certificate_pass',[],'agg_cert_tol',[], ...
    'raw_agg_output_objective_gap',[],'raw_agg_reference_diagnostic',[], ...
    'obj_comm_exact',[],'obj_comm_out',[],'obj_comm_exact_meta',[], ...
    'trace_violation',[],'trace_h',[],'trace_obj',[],'trace_lambda',[],'trace_lambda_used',[], ...
    'trace_lambda_step',[],'trace_cv_held_kW',[],'trace_cv_exact_sync_kW',[], ...
    'trace_cv_avg_kW',[],'trace_cv_rolling_kW',[],'trace_lambda_error',[], ...
    'trace_max_age',[],'trace_mean_age',[],'trace_exact_sync_h',[], ...
    'trace_exact_sync_cv_kW',[],'trace_price_update_count',[],'trace_update_event',[], ...
    'trace_averaging_started',[],'trace_avg_cnt',[],'trace_avg_start_stable_count',[], ...
    'trace_formal_cv_stable_pass',[], ...
    'residual_held_final_kW',[],'cv_held_final_kW',[],'residual_avg_final_kW',[], ...
    'cv_avg_final_kW',[],'cv_rolling_final_kW',[],'cv_exact_sync_final_kW',[], ...
    'lambda_exact',[],'lambda_exact_result',[],'final_lambda_error',[], ...
    'max_observed_EU_age',[],'mean_observed_EU_age',[],'min_update_count',[], ...
    'max_update_count',[],'all_EU_updated',[],'age_bound_h0',[], ...
    'age_bound_pass_h0',[],'age_bound_pass_h0_plus_1',[],'diagnostic_record_every',[], ...
    'exact_sync_diagnostic_every',[],'rolling_window',[],'avg_start_policy',[]);
val_pi_l_pos_au(:,k) = init_pi_l_pos;
val_pi_l_neg_au(:,k) = init_pi_l_neg;
val_pi_PB_au(:,k) = init_pi_PB;

outer_rng_seed = rng_seeds.outer_delay;
outer_delay_stream = RandStream('mt19937ar','Seed',outer_rng_seed);
s3_rng_config = struct('inner_delay_seed',rng_seeds.inner_delay, ...
    'qre_noise_seed',rng_seeds.qre_noise,'qre_audit_seed',rng_seeds.qre_audit);
s3_avg_config = struct('policy',inner_avg_policy, ...
    'policy_version',inner_avg_policy_version, ...
    'start_cv_factor',inner_avg_start_cv_factor, ...
    'start_price_factor',inner_avg_start_price_factor, ...
    'start_stable_window',inner_avg_start_stable_window, ...
    'min_price_updates',inner_avg_min_price_updates, ...
    'min_samples',inner_avg_min_samples);
s3_avg_config.formal_cv_stable_window = inner_formal_cv_stable_window;
k0 = k0_max(:);
if numel(k0) ~= num_LESMs || any(k0 < 0) || any(k0 ~= floor(k0))
    error('params.mat k0_max must be a nonnegative integer delay bound per community.');
end

delay = zeros(num_LESMs,1);
count = zeros(num_LESMs,1);
community_local_call_count = zeros(num_LESMs,1);
for i = 1:num_LESMs
    delay(i) = randi(outer_delay_stream,[0, k0(i)]);%deta
end
%% 开始迭代
outer_history_capacity = min(max_outer_iter + 2, 128);
outer_price_step = nan(outer_history_capacity, 1);
% 线路诊断只用于定位停止缓慢的来源，不参与任何价格更新或停止判断。
finite_line_idx = find(finite_line_mask);
num_finite_lines = numel(finite_line_idx);
line_inst_ratio_history = nan(num_finite_lines, outer_history_capacity);
line_avg_ratio_history = nan(num_finite_lines, outer_history_capacity);
line_dual_pos_history = nan(num_finite_lines, outer_history_capacity);
line_dual_neg_history = nan(num_finite_lines, outer_history_capacity);
over_trace = 0;
inner_k = struct('h',[], ...
    'val_violation',[], ...
    'val_obj',[], ...
    'obj_ref',[], ...
    'RE_obj',[]);
total_inner_iter = 0;
total_inner_price_updates = 0;
num_local_call = 0;
max_inner_wall_iter_observed = 0;
max_inner_price_updates_observed = 0;
total_local_time = 0;

%% 上层遍历平均：只恢复原变量，不参与价格更新
sum_outer_pes  = zeros(num_prosumers,1);
sum_outer_pg   = zeros(num_prosumers,1);
sum_outer_pb   = zeros(num_prosumers,1);
sum_outer_ps   = zeros(num_prosumers,1);
sum_outer_pesc = zeros(num_LESMs,1);
outer_avg_cnt = 0;
val_totcost_bar = nan(outer_history_capacity, 1);
val_outer_violation_bar = nan(outer_history_capacity, 1);
pb_cv_history_kW = nan(outer_history_capacity, 1);
max_line_cv_history_kW = nan(outer_history_capacity, 1);
componentwise_cv_history = nan(outer_history_capacity, 1);

while true
    if k + 1 > outer_history_capacity
        nextCapacity = growHistoryCapacity(outer_history_capacity, k + 1, max_outer_iter + 2);
        outer_price_step(nextCapacity,1) = NaN;
        val_totcost_bar(nextCapacity,1) = NaN;
        val_outer_violation_bar(nextCapacity,1) = NaN;
        pb_cv_history_kW(nextCapacity,1) = NaN;
        max_line_cv_history_kW(nextCapacity,1) = NaN;
        componentwise_cv_history(nextCapacity,1) = NaN;
        line_inst_ratio_history(:,nextCapacity) = NaN;
        line_avg_ratio_history(:,nextCapacity) = NaN;
        line_dual_pos_history(:,nextCapacity) = NaN;
        line_dual_neg_history(:,nextCapacity) = NaN;
        outer_history_capacity = nextCapacity;
    end

    beta_k = outer_momentum;
    if progress_print_enabled && (k == 1 || mod(k, progress_print_outer_every) == 0)
        fprintf('[S3 outer start] k=%d, completed local calls=%d\n', k, num_local_call);
    end
    for i = 1:num_LESMs
        record_inner_or_not = ismember(k,record_k) && (i == record_inner_comm) && (count(i) == delay(i)) && ~over_trace;%等于1：记录；等于0：不记录
        if count(i) < delay(i)
            count(i) = count(i) + 1;
        else
            local_tic = tic;
            community_local_call_count(i) = community_local_call_count(i) + 1;
            res(i) = localMarket_S3(c(start_idx(i):end_idx(i)), ...
                b(start_idx(i):end_idx(i)), a(i), ...
                D(start_idx(i):end_idx(i)), init_rho(i), ...
                pg_max(start_idx(i):end_idx(i)), pg_min(start_idx(i):end_idx(i)), ...
                pi(i), ...
                pes(start_idx(i):end_idx(i)), ...
                pg(start_idx(i):end_idx(i)), ...
                pb(start_idx(i):end_idx(i)), ...
                ps(start_idx(i):end_idx(i)), ...
                pesc(i), ...
                errTol_LESMs(i), inner_cv_tol_kW(i), pi_0(i), pi_max, pi_min, record_inner_or_not, ...
                beta_qre(start_idx(i):end_idx(i)), qre_epsilon, agg_epsilon_i(i), qre_z_cap, ...
                qre_backoff_factor, qre_max_backoffs, inner_momentum, max_inner_iter, h0_max, ...
                min_inner_iter, inner_avg_start_iter, inner_cv_stop_enabled, ...
                qre_certificate_enabled, agg_certificate_enabled, agg_cert_tol, stable_inner_window, ...
                qre_audit_enabled, qre_audit_rate, qre_audit_trace_community_full, qre_audit_seed, ...
                agg_gap_diagnostic_enabled, qre_noise_enabled, diagnostic_record_every, ...
                exact_sync_diagnostic_every, rolling_window, progress_print_enabled, ...
                progress_print_inner_every, k, i,s3_rng_config,s3_avg_config, ...
                community_local_call_count(i));

            if ~res(i).success
                error('upperMarketS3:localFailure', ...
                    'Community %d failed: %s', i, res(i).stop_reason);
            end
            if qre_certificate_enabled && ~res(i).qre_all_certificate_pass
                error('upperMarketS3:qreCertificateFailure', ...
                    'Community %d returned an uncertified QRE response.', i);
            end
            if agg_certificate_enabled && ~res(i).agg_certificate_pass
                error('upperMarketS3:aggregateCertificateFailure', ...
                    'Community %d exceeded its aggregate error budget.', i);
            end

            local_elapsed = toc(local_tic);
            total_inner_iter = total_inner_iter + res(i).iter;
            total_inner_price_updates = total_inner_price_updates + res(i).price_updates;
            num_local_call = num_local_call + 1;
            max_inner_wall_iter_observed = max(max_inner_wall_iter_observed, res(i).iter);
            max_inner_price_updates_observed = max(max_inner_price_updates_observed, res(i).price_updates);
            total_local_time = total_local_time + local_elapsed;
            % 输出
            pi(i) = res(i).pi;
            pes(start_idx(i):end_idx(i)) = res(i).pes;
            pg(start_idx(i):end_idx(i))  = res(i).pg;
            pb(start_idx(i):end_idx(i))  = res(i).pb;
            ps(start_idx(i):end_idx(i))  = res(i).ps;
            pesc(i) = res(i).pesc;

            count(i) = 0;
            delay(i) = randi([0,k0(i)]);
            % 记录时，RE 与 CV 均使用同一份遍历平均社区输出。
            if record_inner_or_not
                [obj_ref_scalar, obj_ref_meta] = exactLowerObjFast( ...
                        c(start_idx(i):end_idx(i)), ...
                        b(start_idx(i):end_idx(i)), ...
                        a(i), ...
                        D(start_idx(i):end_idx(i)), ...
                        pg_max(start_idx(i):end_idx(i)), ...
                        pg_min(start_idx(i):end_idx(i)), ...
                        pi_0(i), pi_max, pi_min);

                obj_trace = res(i).trace_obj(:);
                violation_trace = res(i).trace_violation(:);
                h_trace = res(i).trace_h(:);
                obj_ref = repmat(obj_ref_scalar, size(obj_trace));

                obj_scale = max(abs(obj_ref_scalar), eps);
                RE_obj_trace = 100 * abs(obj_trace - obj_ref) ./ obj_scale;

                inner_k.h = h_trace;
                inner_k.price_update_count = res(i).trace_price_update_count(:);
                inner_k.update_event = res(i).trace_update_event(:);
                inner_k.val_violation = violation_trace;
                inner_k.val_obj = obj_trace;
                inner_k.obj_ref = obj_ref;
                inner_k.RE_obj = RE_obj_trace;
                inner_k.CV_inner_avg_kW = violation_trace;
                inner_k.RE_inner_avg = RE_obj_trace;
                inner_k.trace_metric = 'ergodic_community_primal_objective';
                inner_k.reference_metric = 'exact_community_primal_optimum_at_fixed_pi0';
                inner_k.obj_ref_meta = obj_ref_meta;
                 inner_k.qre_error_diagnostic = struct( ...
                     'qre_gap_max',res(i).qre_gap_max, ...
                     'qre_gap_mean',res(i).qre_gap_mean, ...
                     'qre_gamma_min',res(i).qre_gamma_min, ...
                     'qre_gamma_mean',res(i).qre_gamma_mean, ...
                     'qre_backoff_total',res(i).qre_backoff_total, ...
                     'qre_backoff_max',res(i).qre_backoff_max, ...
                     'qre_fallback_count',res(i).qre_fallback_count, ...
                     'qre_branch_switch_count',res(i).qre_branch_switch_count);
                inner_k.record_outer_iter = k;
                inner_k.record_comm = i;
                inner_k.final_RE_obj = RE_obj_trace(end);
                inner_k.final_val_violation = violation_trace(end);
                inner_k.final_RE_inner_avg = RE_obj_trace(end);
                inner_k.final_CV_inner_avg_kW = violation_trace(end);
                inner_k.record_h = h_trace(end);
                inner_k.record_h_used = h_trace(end);
                inner_k.record_index = numel(h_trace);
                inner_k.record_source = 'community_primal_trace_on_global_wall_clock';
                inner_k.record_extrapolated_by_hold = false;
                inner_k.record_RE_obj = inner_k.final_RE_obj;
                inner_k.record_val_violation = inner_k.final_val_violation;
                inner_k.own_final_RE_obj = inner_k.final_RE_obj;
                inner_k.own_final_val_violation = inner_k.final_val_violation;

                over_trace = 1;

            end
        end
    end

    %% HDEM式上层运行平均：从第一个外层更新起计入输出。
    if k >= outer_avg_start_iter
        outer_avg_cnt = outer_avg_cnt + 1;

        sum_outer_pes  = sum_outer_pes  + pes;
        sum_outer_pg   = sum_outer_pg   + pg;
        sum_outer_pb   = sum_outer_pb   + pb;
        sum_outer_ps   = sum_outer_ps   + ps;
        sum_outer_pesc = sum_outer_pesc + pesc;

        bar_outer_pes  = sum_outer_pes  / outer_avg_cnt;
        bar_outer_pg   = sum_outer_pg   / outer_avg_cnt;
        bar_outer_pb   = sum_outer_pb   / outer_avg_cnt;
        bar_outer_ps   = sum_outer_ps   / outer_avg_cnt;
        bar_outer_pesc = sum_outer_pesc / outer_avg_cnt;
    else
        bar_outer_pes = pes;
        bar_outer_pg = pg;
        bar_outer_pb = pb;
        bar_outer_ps = ps;
        bar_outer_pesc = pesc;
    end

    val_totcost_bar(k+1) = sum(0.5 .* c .* bar_outer_pg.^2 ...
        + b .* bar_outer_pg ...
        + pi_max .* bar_outer_pb ...
        - pi_min .* bar_outer_ps ...
        + 0.5 .* a_extend .* bar_outer_pes.^2) ...
        + sum(0.5 .* a .* bar_outer_pesc.^2);

    PB_bar_now_kW = abs(sum(bar_outer_pesc));
    LINE_pos_bar_now_kW = max(pai_il' * bar_outer_pesc - F_l, 0);
    LINE_neg_bar_now_kW = max(-pai_il' * bar_outer_pesc - F_l, 0);
    val_outer_violation_bar(k+1) = norm([PB_bar_now_kW; LINE_pos_bar_now_kW; LINE_neg_bar_now_kW], 2);
    line_violation_bar_kW = max(LINE_pos_bar_now_kW, LINE_neg_bar_now_kW);
    % 与遍历平均诊断并列记录当前 raw Agg 上传的瞬时线路违反。
    LINE_pos_now_kW = max(pai_il' * pesc - F_l, 0);
    LINE_neg_now_kW = max(-pai_il' * pesc - F_l, 0);
    line_violation_now_kW = max(LINE_pos_now_kW, LINE_neg_now_kW);
    line_inst_ratio = line_violation_now_kW(finite_line_mask) ./ F_l(finite_line_mask);
    line_avg_ratio = line_violation_bar_kW(finite_line_mask) ./ F_l(finite_line_mask);
    line_inst_ratio_history(:,k) = line_inst_ratio;
    line_avg_ratio_history(:,k) = line_avg_ratio;
    pb_cv_history_kW(k) = PB_bar_now_kW;
    max_line_cv_history_kW(k) = max(line_violation_bar_kW(finite_line_mask));
    componentwise_cv_history(k) = max([PB_bar_now_kW / system_balance_scale_kW; line_avg_ratio]);
    outer_pb_cv_pass = PB_bar_now_kW <= outer_pb_cv_tol_kW;
    outer_line_cv_pass = all(line_violation_bar_kW(finite_line_mask) <= ...
        outer_line_cv_tol_kW(finite_line_mask));
    if progress_print_enabled && (k == 1 || mod(k, progress_print_outer_every) == 0)
        max_line_ratio_now = max(line_violation_bar_kW(finite_line_mask) ./ F_l(finite_line_mask));
        fprintf('[S3 outer] k=%d, PBavg=%.3f/%.3f kW, maxLineCV=%.4f, L2diag=%.3f kW\n', ...
            k, PB_bar_now_kW, outer_pb_cv_tol_kW, max_line_ratio_now, val_outer_violation_bar(k+1));
        [max_inst_ratio, max_inst_pos] = max(line_inst_ratio);
        [max_avg_ratio, max_avg_pos] = max(line_avg_ratio);
        fprintf('[S3 line] k=%d, inst(line %d)=%.4f, avg(line %d)=%.4f\n', ...
            k, finite_line_idx(max_inst_pos), max_inst_ratio, ...
            finite_line_idx(max_avg_pos), max_avg_ratio);
    end

    %% 上层运营者更新基础电价 using instantaneous held raw Agg
    val_pi_PB_au(1,k+1) = val_pi_PB(1,k) - alpha_PB * sum(pesc);
    val_pi_PB(1,k+1) = val_pi_PB_au(1,k+1);

    %计算pi_l——这个是不等式约束限制的价格，所以必须投影（122行）
    raw_pi_l_pos = val_pi_l_pos(:,k) - alpha_l * (pai_il' * pesc - F_l);
    raw_pi_l_neg = val_pi_l_neg(:,k) - alpha_l * (-pai_il' * pesc - F_l);
    val_pi_l_pos_au(:,k+1) = min(raw_pi_l_pos, 0);
    val_pi_l_neg_au(:,k+1) = min(raw_pi_l_neg, 0);
    val_pi_l_pos(:,k+1) = val_pi_l_pos_au(:,k+1);
    val_pi_l_neg(:,k+1) = val_pi_l_neg_au(:,k+1);
    line_dual_pos_history(:,k) = val_pi_l_pos(finite_line_mask,k+1);
    line_dual_neg_history(:,k) = val_pi_l_neg(finite_line_mask,k+1);

    %由pi_l、pi_PB计算pi_0
    val_pi_0(:,k+1) = val_pi_PB(1,k+1) + pai_il * (val_pi_l_pos(:,k+1) - val_pi_l_neg(:,k+1));

    pi_0 = val_pi_0(:,k+1);

    % if norm(val_pi_0(:,k+1) - val_pi_0(:,k),inf) <= errTol_UESM
    % 外层当前约束残差，单位 kW。
    PB_res_kW = abs(sum(pesc));
    LINE_pos_res_kW = max(pai_il' * pesc - F_l, 0);
    LINE_neg_res_kW = max(-pai_il' * pesc - F_l, 0);
    CV_U_now_kW = norm([PB_res_kW; LINE_pos_res_kW; LINE_neg_res_kW], 2);

    % 当前配置下 CV 仅作逐 k 诊断；可通过开关恢复价格步长与平均 CV 双条件。
    outer_price_step(k) = norm(val_pi_0(:,k+1) - val_pi_0(:,k), inf);

    % 不设人为迭代下限；outer_avg_start_iter 仅定义后缀平均开始时刻。
    if outer_price_step(k) <= errTol_UESM
        stable_outer_price_count = stable_outer_price_count + 1;
    else
        stable_outer_price_count = 0;
    end
    outer_cv_pass = (~logical(outer_cv_stop_enabled)) || ...
        (outer_pb_cv_pass && outer_line_cv_pass);
    if outer_cv_l2_stop_enabled
        outer_cv_pass = outer_cv_pass && ...
            (val_outer_violation_bar(k+1) <= outer_cv_l2_tol_kW);
    end
    if k >= max(min_outer_iter, outer_avg_start_iter) && ...
            stable_outer_price_count >= stable_outer_window && outer_cv_pass
        if logical(outer_cv_stop_enabled)
            outer_stop_reason = 'stable_price_and_componentwise_relative_cv';
        else
            outer_stop_reason = 'price_step_only';
        end
        break;
    end

    if isfinite(max_outer_iter) && k >= max_outer_iter
        state.success = 0;
        outer_stop_reason = 'max_outer_iter';
        break;
    end

    k = k + 1;
end
final_gap_kW = zeros(num_LESMs,1);
for i = 1:num_LESMs
    final_gap_kW(i) = abs(bar_outer_pesc(i) - sum(bar_outer_pes(start_idx(i):end_idx(i))));
end

outer = struct('k',[], ...
    'val_violation',[], ...
    'totcost',[], ...
    'totcost_ref',[], ...
    'RE_obj',[], ...
    'pi',[], ...
    'pesc',[], ...
    'pi_0',[], ...
    'final_gap_kW',[] ...
    );
outer.final_gap_kW = final_gap_kW;
%% 计算外层遍历平均点的目标函数值相对误差
totcost_ref = exactUpperObj(a,c,b,pi_max,pi_min,D,F_l,pai_il, ...
    start_idx,end_idx,pg_max,pg_min,a_extend);
totcost_scale = max(abs(totcost_ref), eps);

outer_iter = (1:k)';
totcost_iter = val_totcost_bar(2:k+1);
violation_iter = val_outer_violation_bar(2:k+1);
RE_totcost = 100 * abs(totcost_iter - totcost_ref) / totcost_scale;
%% 计算外层遍历平均点的约束违反，单位 kW。
PB_bar_kW = abs(sum(bar_outer_pesc));
LINE_pos_bar_kW = max(pai_il' * bar_outer_pesc - F_l, 0);
LINE_neg_bar_kW = max(-pai_il' * bar_outer_pesc - F_l, 0);
CV_kW = norm([PB_bar_kW; LINE_pos_bar_kW; LINE_neg_bar_kW], 2);
line_kW = max([0; LINE_pos_bar_kW; LINE_neg_bar_kW]);
line_violation_kW = max(LINE_pos_bar_kW, LINE_neg_bar_kW);
line_cv_ratio = zeros(size(F_l));
line_cv_ratio(finite_line_mask) = line_violation_kW(finite_line_mask) ./ F_l(finite_line_mask);
PB_cv_ratio = PB_bar_kW / system_balance_scale_kW;
componentwise_CV = max([PB_cv_ratio; line_cv_ratio(finite_line_mask)]);
componentwise_CV_pass = PB_bar_kW <= outer_pb_cv_tol_kW && ...
    all(line_violation_kW(finite_line_mask) <= outer_line_cv_tol_kW(finite_line_mask));

%% 外层存档
outer.k = outer_iter;
outer.iter = outer_iter;
outer.val_violation = violation_iter;
outer.totcost = totcost_iter;
outer.totcost_ref = totcost_ref;
outer.RE_totcost = RE_totcost;
outer.CV_l2_history_kW = violation_iter;
% 异步保持的 k 本身就是全局墙钟；补齐与 S0/S2 相同的绘图字段。
outer.wall_k = outer_iter;
outer.RE_totcost_wall = RE_totcost;
outer.CV_l2_history_wall_kW = violation_iter;
outer.PB_history_wall_kW = pb_cv_history_kW(1:k);
outer.max_finite_line_cv_history_wall_kW = max_line_cv_history_kW(1:k);
outer.wall_price_update_count = outer_iter;
outer.wall_update_event = true(k,1);
outer.PB_history_kW = pb_cv_history_kW(1:k);
outer.max_finite_line_cv_history_kW = max_line_cv_history_kW(1:k);
outer.componentwise_CV_history = componentwise_cv_history(1:k);
outer.outer_history_semantics = ['each k is a global wall-clock event: completed communities upload new suffix-average ', ...
    'raw Agg and unfinished communities hold their last upload; RE/CV are diagnostics, not stop criteria'];
outer.final_RE_totcost = RE_totcost(end);
outer.final_val_violation = CV_kW;
outer.record_iter = numel(RE_totcost);
outer.record_iter_used = numel(RE_totcost);
outer.record_source = 'own_final_iter';
outer.record_extrapolated_by_hold = false;
outer.record_RE_totcost = outer.final_RE_totcost;
outer.record_val_violation = outer.final_val_violation;
outer.own_final_RE_totcost = outer.final_RE_totcost;
outer.own_final_val_violation = outer.final_val_violation;
outer.CV_kW = CV_kW;
outer.PB_kW = PB_bar_kW;
outer.line_kW = line_kW;
outer.raw_l2_cv_diagnostic_only = true;
outer.raw_agg_objective_diagnostic = true;
outer.l2_cv_stop_enabled = logical(outer_cv_l2_stop_enabled);
outer.cv_policy = cv_policy;
outer.inner_cv_ratio = inner_cv_ratio;
outer.inner_cv_tol_kW = inner_cv_tol_kW;
outer.community_load_scale_kW = community_load_scale_kW;
outer.outer_cv_ratio = outer_cv_ratio;
outer.system_balance_scale_kW = system_balance_scale_kW;
outer.outer_pb_cv_tol_kW = outer_pb_cv_tol_kW;
outer.outer_line_cv_tol_kW = outer_line_cv_tol_kW;
outer.final_PB_cv_ratio = PB_cv_ratio;
outer.final_line_cv_ratio = line_cv_ratio;
outer.final_componentwise_CV = componentwise_CV;
outer.final_componentwise_CV_pass = componentwise_CV_pass;
outer.final_gap_mean_kW = mean(outer.final_gap_kW,'omitnan');
outer.final_gap_max_kW = max(outer.final_gap_kW);
outer.final_gap_p95_kW = prctile(outer.final_gap_kW,95);
outer.outer_avg_cnt = outer_avg_cnt;
outer.stable_outer_window = stable_outer_window;
outer.stable_outer_price_count = stable_outer_price_count;
outer.outer_cv_l2_tol_kW = outer_cv_l2_tol_kW;
outer.qre_certificate_enabled = qre_certificate_enabled;
outer.run_profile = run_profile;
if qre_certificate_enabled
    outer.qre_mode = 'certified';
else
    outer.qre_mode = 'fast_qre_with_posterior_audit';
end
outer.qre_epsilon = qre_epsilon;
outer.agg_certificate_enabled = agg_certificate_enabled;
outer.agg_epsilon_i = agg_epsilon_i;
outer.agg_cert_tol = agg_cert_tol;
outer.qre_gap_max = max([res.qre_gap_max]);
outer.qre_gap_mean = mean([res.qre_gap_mean],'omitnan');
outer.qre_gap_p95 = max([res.qre_gap_p95]);
outer.qre_gap_exceed_count = sum([res.qre_gap_exceed_count]);
outer.qre_gap_exceed_rate = outer.qre_gap_exceed_count / ...
    max(sum([res.qre_audit_count]),1);
outer.qre_audit_count = sum([res.qre_audit_count]);
outer.qre_audit_response_count = sum([res.qre_audit_response_count]);
outer.qre_audit_enabled = qre_audit_enabled;
outer.qre_audit_rate = qre_audit_rate;
outer.qre_audit_seed = qre_audit_seed;
outer.qre_gamma_min = min([res.qre_gamma_min]);
outer.qre_gamma_mean = mean([res.qre_gamma_mean],'omitnan');
outer.qre_backoff_total = sum([res.qre_backoff_total]);
outer.qre_backoff_max = max([res.qre_backoff_max]);
outer.qre_fallback_count = sum([res.qre_fallback_count]);
outer.qre_branch_switch_count = sum([res.qre_branch_switch_count]);
outer.qre_branch_switch_rate = outer.qre_branch_switch_count / ...
    max(sum([res.qre_response_count]),1);
outer.agg_gap_max = max([res.agg_gap]);
outer.agg_gap_mean = mean([res.agg_gap],'omitnan');
outer.raw_agg_output_objective_gap_max = max([res.raw_agg_output_objective_gap]);
outer.raw_agg_output_objective_gap_mean = mean([res.raw_agg_output_objective_gap],'omitnan');
outer.raw_agg_reference_diagnostic = true;
outer.outer_price_step = outer_price_step(1:k);
outer.finite_line_idx = finite_line_idx;
outer.line_ratio_instantaneous = line_inst_ratio_history(:,1:k);
outer.line_ratio_average = line_avg_ratio_history(:,1:k);
outer.line_dual_pos = line_dual_pos_history(:,1:k);
outer.line_dual_neg = line_dual_neg_history(:,1:k);
outer.line_diagnostic_semantics = ['instantaneous=raw Agg upload at outer k; ', ...
    'average=suffix-average raw Agg upload; dual=post-update multiplier at k'];
outer.min_outer_iter = min_outer_iter;
outer.max_outer_iter = max_outer_iter;
outer.errTol_UESM = errTol_UESM;
outer.inner_cv_diagnostic_reference_kW = inner_cv_tol_kW;
outer.inner_avg_start_iter = inner_avg_start_iter;
outer.inner_cv_stop_enabled = inner_cv_stop_enabled;
outer.outer_cv_stop_enabled = outer_cv_stop_enabled;
outer.stop_reason = outer_stop_reason;
outer.outer_avg_start_iter = outer_avg_start_iter;
outer.random_seed = outer_rng_seed;
outer.state = state;
outer.num_local_call = num_local_call;
outer.community_local_call_count = community_local_call_count;
outer.seed_derivation_policy = ...
    'deriveDeterministicSeed(base_seed,community_id,local_call_count,stream_tag)';
outer.total_inner_iter = total_inner_iter;
outer.total_inner_price_updates = total_inner_price_updates;
outer.max_inner_iter = max_inner_wall_iter_observed;
outer.max_inner_wall_iter = max_inner_wall_iter_observed;
outer.max_inner_price_updates = max_inner_price_updates_observed;
outer.total_local_time = total_local_time;
outer.avg_inner_iter = total_inner_iter / max(num_local_call, 1);
outer.avg_inner_price_updates = total_inner_price_updates / max(num_local_call, 1);
% 价格为最后迭代值；原变量为上层遍历平均输出。
outer.pi = pi(:);
outer.pi_0 = pi_0(:);
outer.pes = bar_outer_pes;
outer.pg = bar_outer_pg;
outer.pb = bar_outer_pb;
outer.ps = bar_outer_ps;
outer.pesc = bar_outer_pesc;


if save_s3_results
    % 用户指定：S3 每次仅覆盖同一对固定文件，避免累积时间戳结果。
    inner_result_file = 'inner_k_data_3.mat';
    outer_result_file = 'outer_data_3.mat';
    save(inner_result_file,'inner_k','params_snapshot','run_tag','result_unit');
    save(outer_result_file,'outer','params_snapshot','run_tag','result_unit');
else
    % 保留 inner_k、outer、params_snapshot 于当前 MATLAB 工作区，供即时查看。
    inner_result_file = '';
    outer_result_file = '';
end
toc;
%% 对照组
% totCost_SS = SS(c,b,pi_max,pi_min,D,pg_max,pg_min);
% totCost_LS = LS(num_LESMs,start_idx,end_idx,c,b,a,D,init_rho,pg_max,pg_min,pi,pes,errTol_LESMs,pi_max,pi_min,a_extend);
% totCost_LO = LO(num_LESMs, start_idx, end_idx, c, b, a, D, pg_max, pg_min, pi_max, pi_min);
% totCost_GS = val_totcost(1,k);
% totCost_GO = GO(c,b,pi_max,pi_min,D,F_l,pai_il,start_idx,end_idx,pg_max,pg_min);
% %生成table
% totcosts_k = [totCost_SS, totCost_LS, totCost_LO, totCost_GS, totCost_GO] / 1000;
% fprintf('\n\n');
% fprintf('Total cost in different conditions.\n');
% fprintf('------------------------------------------------------------\n');
% fprintf('%-16s %-8s %-8s %-8s %-8s %-8s\n', 'Conditions', 'SS', 'LS', 'LO', 'GS', 'GO');
% fprintf('------------------------------------------------------------\n');
% fprintf('%-16s %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f\n', 'Total cost/k$', ...
%     totcosts_k(1), totcosts_k(2), totcosts_k(3), totcosts_k(4), totcosts_k(5));
% fprintf('------------------------------------------------------------\n\n');
% %% 画图：1.几个节点的每次迭代的相对误差；2.几个节点每次迭代的能源共享价格pes
% colors = {'#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30'};
% %1.相对误差
% figure;
% subplot(2,2,1);
% hold on; box on;
% for N = 1:length(num_list)
%     m = num_list(N);
%     valid_len = inner_iter_len(m);
%     graph_RE_pi = RE_pi(m,1:valid_len);
%     plot(graph_RE_pi,'-x','Color',colors{N}, ...
%         'DisplayName',sprintf('L-ESM %d',m),'LineWidth', ...
%         1.5,'MarkerSize', 6,'MarkerIndices', 1:5:length(graph_RE_pi));
% end
% xlabel('iteration');
% ylabel('relative_error / %');
% lgd = legend(subplot(2,2,1), 'Orientation', 'horizontal');
% lgd.Box = 'off';
% lgd.NumColumns = 5;
% set(lgd, 'Position', [0.25, 0.94, 0.5, 0.05]);
% %2.能源共享价格
% subplot(2,2,2);
% hold on; box on;
% for N = 1:length(num_list)
%     m = num_list(N);
%     graph_val_pi = val_pi(m,1:k+1);
%     plot(graph_val_pi,'-x','Color',colors{N}, ...
%         'DisplayName',sprintf('L-ESM %d',m), ...
%         'LineWidth',1.5, ...
%         'MarkerSize',6, ...
%         'MarkerIndices',1:20:length(graph_val_pi));
% end
% xlabel('iteration');
% ylabel('sharing price / ($/kW)');
% %% 画图
% idx = 1:40:k;
% %1.相对误差
% subplot(2,2,3);
% hold on; box on;
% exac_totcost = exac_totalCost(a,c,b,pi_max,pi_min,D,F_l,pai_il,start_idx,end_idx,pg_max,pg_min,a_extend);
% RE_totcost = 100 * abs(val_totcost(1,2:k+1) - exac_totcost) / exac_totcost;
% plot(2:k+1,RE_totcost,'-x','Color','k', ...
%     'LineWidth',1.5, 'MarkerSize', 6,'MarkerIndices',idx);
% % plot(100 * abs(xj2_val_totcost(1,:) - exac_totcost) / exac_totcost,'-x','Color','r', ...
% %         'LineWidth',0.5, 'MarkerSize', 2,'MarkerIndices',idx)
% xlabel('iteration');
% ylabel('relative error / %');
% %2.能源共享基础价格
% subplot(2,2,4);
% hold on; box on;
% for N = 1:length(num_list)
%     m = num_list(N);
%     graph_val_pi_0 = val_pi_0(m,1:k+1);
%     plot(graph_val_pi_0,'-x','Color',colors{N}, ...
%         'DisplayName',sprintf('L-ESM %d',m), ...
%         'LineWidth',1.5, ...
%         'MarkerSize',6, ...
%         'MarkerIndices',1:20:length(graph_val_pi_0));
% end
% xlabel('iteration');
% ylabel('base price / ($/kW)');
% %% 算法1、2结合得到每个社区的情况
% figure;
% x = 1:num_LESMs;
% yyaxis left;
% ba = bar(x,pesc .* 1e-3,1.0);
% ylabel('Share energy / kW');
% yyaxis right;
% p1 = plot(x,pi_0,'--','LineWidth',1.5,'Color','#FF7F11');hold on;
% p2 = plot(x,pi,'-','LineWidth',1.5,'Color','#F6CE71');
% ylabel('Prices / ($/kW)');
% xlabel('Lower-layer market index');
% grid on; box on;
% xlim([x(1)-0.5, x(end)+0.5]);%防止柱子贴边
% yyaxis left;ylim([-4,6]);
% yyaxis right;ylim([0,0.25]);
% lgd = legend([ba p1 p2], {'Shared energy','Base price','Sharing price'}, ...
%     'Location','northoutside','Orientation','horizontal');
% lgd.Box = 'off';%图例外面的框
% ax = gca;
% ax.YAxis(1).Color = 'k';   % 左轴黑色
% ax.YAxis(2).Color = 'k';   % 右轴黑色
% 
% 
% figure;
% plot(1:k, cnt00_outer(1:k), '-x', 'LineWidth', 1.5);
% xlabel('outer iteration k');
% ylabel('number of f1=0 and f2=0');
% grid on;
% box on;
% 


     





