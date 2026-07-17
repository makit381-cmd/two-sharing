clc;clear;
tic;
load('params.mat','a','b','c','D','start_idx','end_idx', ...
    'num_LESMs','num_prosumers','pi_max','pi_min','u','F_l', ...
    'alpha_l','alpha_PB','PTDF','init_rho','h0_max','k0_max', ...
    'errTol_LESMs','errTol_UESM', ...
    'pg_max','pg_min','a_extend','beta_qre','qre_epsilon','agg_epsilon','qre_z_cap','qre_max_backoffs','qre_certificate_enabled','agg_certificate_enabled','outer_momentum','inner_momentum','max_inner_iter', ...
    'min_inner_iter','inner_cv_tol_kW','inner_avg_start_iter', ...
    'max_outer_iter','min_outer_iter','outer_cv_tol_kW', ...
    'outer_avg_start_iter','record_inner_comm','record_k','save_s3_results','run_tag','result_unit','rho_qre','sigma_qre', ...
    'small','big','min_x_scale','rng_seeds', ...
    'init_pb','init_pes','init_pesc','init_pg','init_xj2','init_pi','init_pi_0', ...
    'init_pi_l_neg','init_pi_l_pos','init_pi_PB','init_ps');

params_snapshot = load('params.mat');

pai_il = PTDF;%具体的量化网络约束系数    123*122
% %% 准备
% tic;tok;
state.success = 1;
outer_stop_reason = 'running';
k = 1;
val_pi_PB(1,k) = init_pi_PB;
val_pi_l_pos(:,k) = init_pi_l_pos;
val_pi_l_neg(:,k) = init_pi_l_neg;
val_pi_0(:,k) = init_pi_0;
exac_pi = zeros(num_LESMs,1);
exac_pi_0 = zeros(num_LESMs,1);
RE_pi = zeros(num_LESMs,500);
pi = init_pi;
pes = init_pes;
pesc = init_pesc;
xj2 = init_xj2;
pg = init_pg;
pb = init_pb;
ps = init_ps;
pi_0 = init_pi_0;%迭代过程中的不变的能源共享基础价格

res(num_LESMs) = struct('success',[],'pi',[],'pes',[],'pg',[], ...
    'pb',[],'ps',[],'pesc',[],'pesc_agg_raw',[],'pes_sum',[], ...
    'inner_cv_kW',[],'inner_avg_cnt',[],'inner_avg_start_iter',[],'xj2',[],'iter',[],'price_updates',[],'stop_reason',[], ...
    'trace_violation',[],'trace_h',[],'trace_obj',[],'trace_lambda',[],'trace_lambda_used',[], ...
    'qre_calls',[],'qre_perturbed_accepts',[],'qre_backoff_steps',[], ...
    'qre_fallbacks',[],'qre_boundary_hits',[],'qre_max_gap',[],'qre_epsilon',[],'qre_certificate_enabled',[], ...
    'agg_obj_ref',[],'agg_obj_raw',[],'agg_obj_gap',[],'agg_epsilon',[], ...
    'agg_certificate_enabled',[],'agg_objective_certified',[],'agg_a2_approx_certified',[],'agg_exact_meta',[]);
val_pi_l_pos_au(:,k) = init_pi_l_pos;
val_pi_l_neg_au(:,k) = init_pi_l_neg;
val_pi_PB_au(:,k) = init_pi_PB;

outer_rng_seed = rng_seeds.outer_schedule;
rng(outer_rng_seed, 'twister');
k0 = k0_max(:);
if numel(k0) ~= num_LESMs || any(k0 < 0) || any(k0 ~= floor(k0))
    error('params.mat k0_max must be a nonnegative integer delay bound per community.');
end

delay = zeros(num_LESMs,1);
count = zeros(num_LESMs,1);
for i = 1:num_LESMs
    delay(i) = randi([0, k0(i)]);%deta
end
val_pi = init_pi;
val_pes(:,1) = init_pes;
val_pg(:,1) = pg;
val_pb(:,1) = pb;
val_ps(:,1) = ps;
val_pesc(:,1) = init_pesc;
val_xj2(:,1) = init_xj2;
%% 开始迭代
outer_history_capacity = min(max_outer_iter + 2, 128);
outer_price_step = nan(outer_history_capacity, 1);
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

while true
    if k + 1 > outer_history_capacity
        nextCapacity = growHistoryCapacity(outer_history_capacity, k + 1, max_outer_iter + 2);
        outer_price_step(nextCapacity,1) = NaN;
        val_totcost_bar(nextCapacity,1) = NaN;
        val_outer_violation_bar(nextCapacity,1) = NaN;
        outer_history_capacity = nextCapacity;
    end

    beta_k = outer_momentum;
    for i = 1:num_LESMs
        record_inner_or_not = ismember(k,record_k) && (i == record_inner_comm) && (count(i) == delay(i)) && ~over_trace;%等于1：记录；等于0：不记录
        if count(i) < delay(i)
            count(i) = count(i) + 1;
            pi(i) = val_pi(i,k);
            pes(start_idx(i):end_idx(i)) = val_pes(start_idx(i):end_idx(i),k);
            pg(start_idx(i):end_idx(i))  = val_pg(start_idx(i):end_idx(i),k);
            pb(start_idx(i):end_idx(i))  = val_pb(start_idx(i):end_idx(i),k);
            ps(start_idx(i):end_idx(i))  = val_ps(start_idx(i):end_idx(i),k);
            pesc(i) = val_pesc(i,k);
            xj2(i) = val_xj2(i,k);
        else
            local_tic = tic;
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
                errTol_LESMs(i), inner_cv_tol_kW, pi_0(i), pi_max, pi_min, record_inner_or_not, ...
                beta_qre(start_idx(i):end_idx(i)), qre_epsilon, agg_epsilon, qre_z_cap, qre_max_backoffs, ...
                inner_momentum, max_inner_iter, h0_max, ...
                min_inner_iter, inner_avg_start_iter, qre_certificate_enabled, agg_certificate_enabled);

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
            xj2(i) = res(i).xj2;%xj2=sum(xj^2)

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
                inner_k.val_violation = violation_trace;
                inner_k.val_obj = obj_trace;
                inner_k.obj_ref = obj_ref;
                inner_k.RE_obj = RE_obj_trace;
                inner_k.CV_inner_avg_kW = violation_trace;
                inner_k.RE_inner_avg = RE_obj_trace;
                inner_k.trace_metric = 'ergodic_community_primal_objective';
                inner_k.reference_metric = 'exact_community_primal_optimum_at_fixed_pi0';
                inner_k.obj_ref_meta = obj_ref_meta;
                inner_k.qre_error_diagnostic = 'not_evaluated_in_unchecked_qre_profile';
                inner_k.record_outer_iter = k;
                inner_k.record_comm = i;
                inner_k.final_RE_obj = RE_obj_trace(end);
                inner_k.final_val_violation = violation_trace(end);
                inner_k.final_RE_inner_avg = RE_obj_trace(end);
                inner_k.final_CV_inner_avg_kW = violation_trace(end);
                inner_k.record_h = h_trace(end);
                inner_k.record_h_used = h_trace(end);
                inner_k.record_index = numel(h_trace);
                inner_k.record_source = 'ergodic_community_primal_trace';
                inner_k.record_extrapolated_by_hold = false;
                inner_k.record_RE_obj = inner_k.final_RE_obj;
                inner_k.record_val_violation = inner_k.final_val_violation;
                inner_k.own_final_RE_obj = inner_k.final_RE_obj;
                inner_k.own_final_val_violation = inner_k.final_val_violation;

                over_trace = 1;

            end
        end
        % 更新各个社区的值
        val_pi(i,k+1) = pi(i);
        val_pes(start_idx(i):end_idx(i),k+1) = pes(start_idx(i):end_idx(i));
        val_pg(start_idx(i):end_idx(i),k+1) = pg(start_idx(i):end_idx(i));
        val_pb(start_idx(i):end_idx(i),k+1) = pb(start_idx(i):end_idx(i));
        val_ps(start_idx(i):end_idx(i),k+1) = ps(start_idx(i):end_idx(i));
        val_pesc(i,k+1) = pesc(i);
        val_xj2(i,k+1) = xj2(i);
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

    %% 上层运营者更新基础电价
    %计算pi_PB
    val_pi_PB_au(1,k+1) = val_pi_PB(1,k) - alpha_PB * sum(pesc);
    momentum_pi_PB = beta_k * (val_pi_PB_au(1,k+1) - val_pi_PB_au(1,k));
    val_pi_PB(1,k+1) = val_pi_PB_au(1,k+1) + momentum_pi_PB;

    %计算pi_l——这个是不等式约束限制的价格，所以必须投影（122行）
    raw_pi_l_pos = val_pi_l_pos(:,k) - alpha_l * (pai_il' * pesc - F_l);
    raw_pi_l_neg = val_pi_l_neg(:,k) - alpha_l * (-pai_il' * pesc - F_l);
    val_pi_l_pos_au(:,k+1) = min(raw_pi_l_pos, 0);
    val_pi_l_neg_au(:,k+1) = min(raw_pi_l_neg, 0);
    momentum_pi_l_pos = beta_k * (val_pi_l_pos_au(:,k+1) - val_pi_l_pos_au(:,k));
    momentum_pi_l_neg = beta_k * (val_pi_l_neg_au(:,k+1) - val_pi_l_neg_au(:,k));
    val_pi_l_pos(:,k+1) = val_pi_l_pos_au(:,k+1) + momentum_pi_l_pos;
    val_pi_l_neg(:,k+1) = val_pi_l_neg_au(:,k+1) + momentum_pi_l_neg;

    %由pi_l、pi_PB计算pi_0
    val_pi_0(:,k+1) = val_pi_PB(1,k+1) + pai_il * (val_pi_l_pos(:,k+1) - val_pi_l_neg(:,k+1));

    pi_0 = val_pi_0(:,k+1);

    % if norm(val_pi_0(:,k+1) - val_pi_0(:,k),inf) <= errTol_UESM
    % 外层当前约束残差，单位 kW。
    PB_res_kW = abs(sum(pesc));
    LINE_pos_res_kW = max(pai_il' * pesc - F_l, 0);
    LINE_neg_res_kW = max(-pai_il' * pesc - F_l, 0);
    CV_U_now_kW = norm([PB_res_kW; LINE_pos_res_kW; LINE_neg_res_kW], 2);

    % 价格稳定与正式遍历平均上传量的 PB/线路 CV 必须同时满足。
    outer_price_step(k) = norm(val_pi_0(:,k+1) - val_pi_0(:,k), inf);

    if k >= min_outer_iter && outer_price_step(k) <= errTol_UESM && ...
            val_outer_violation_bar(k+1) <= outer_cv_tol_kW
        outer_stop_reason = 'price_step_and_average_cv';
        break;
    end

    if k >= max_outer_iter
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

%% 外层存档
outer.k = outer_iter;
outer.iter = outer_iter;
outer.val_violation = violation_iter;
outer.totcost = totcost_iter;
outer.totcost_ref = totcost_ref;
outer.RE_totcost = RE_totcost;
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
outer.final_gap_mean_kW = mean(outer.final_gap_kW,'omitnan');
outer.final_gap_max_kW = max(outer.final_gap_kW);
outer.final_gap_p95_kW = prctile(outer.final_gap_kW,95);
outer.outer_avg_cnt = outer_avg_cnt;
outer.outer_price_step = outer_price_step(1:k);
outer.min_outer_iter = min_outer_iter;
outer.max_outer_iter = max_outer_iter;
outer.errTol_UESM = errTol_UESM;
outer.inner_cv_diagnostic_reference_kW = inner_cv_tol_kW;
outer.inner_avg_start_iter = inner_avg_start_iter;
outer.errTol_output_CV_kW = outer_cv_tol_kW;
outer.inner_cv_stop_enabled = true;
outer.outer_cv_stop_enabled = true;
outer.stop_reason = outer_stop_reason;
outer.outer_avg_start_iter = outer_avg_start_iter;
outer.random_seed = outer_rng_seed;
outer.state = state;
outer.num_local_call = num_local_call;
outer.total_inner_iter = total_inner_iter;
outer.total_inner_price_updates = total_inner_price_updates;
outer.max_inner_iter = max_inner_wall_iter_observed;
outer.max_inner_wall_iter = max_inner_wall_iter_observed;
outer.max_inner_price_updates = max_inner_price_updates_observed;
outer.total_local_time = total_local_time;
outer.avg_inner_iter = total_inner_iter / max(num_local_call, 1);
outer.avg_inner_price_updates = total_inner_price_updates / max(num_local_call, 1);
outer.agg_epsilon = NaN;
outer.agg_certificate_enabled = false;
outer.agg_objective_certified_calls = sum([res.agg_objective_certified]);
outer.agg_a2_approx_certified_calls = sum([res.agg_a2_approx_certified]);
outer.agg_certified_call_rate = outer.agg_a2_approx_certified_calls / max(num_local_call, 1);
outer.agg_max_obj_gap = NaN;
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


     





