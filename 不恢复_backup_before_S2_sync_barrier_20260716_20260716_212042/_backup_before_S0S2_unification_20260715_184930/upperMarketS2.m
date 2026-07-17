clc;clear;
tic;
% S2: 同步等待处理外部延迟 + 遍历平均输出
load('params.mat','a','b','c','D','start_idx','end_idx', ...
    'num_LESMs','num_prosumers','pi_max','pi_min','u','F_l', ...
    'alpha_l','alpha_PB','PTDF','init_rho', ...
    'errTol_LESMs','errTol_UESM', ...
    'pg_max','pg_min','a_extend','beta_qre','qre_error_rel','qre_z_cap','outer_momentum','inner_momentum','max_inner_iter','rho_qre','sigma_qre', ...
    'small','big','min_x_scale','rng_seeds', ...
    'init_pb','init_pes','init_pesc','init_pg','init_pi','init_pi_0', ...
    'init_pi_l_neg','init_pi_l_pos','init_pi_PB','init_ps');

pai_il = PTDF;
state.success = 1;
k = 1;
val_pi_PB(1,k) = init_pi_PB;
val_pi_l_pos(:,k) = init_pi_l_pos;
val_pi_l_neg(:,k) = init_pi_l_neg;
val_pi_0(:,k) = init_pi_0;

pi = init_pi;
pes = init_pes;
pesc = init_pesc;
xj2 = zeros(num_LESMs, 1);
pg = init_pg;
pb = init_pb;
ps = init_ps;
pi_0 = init_pi_0;

record_inner_comm = 71;

res(num_LESMs) = struct('success',[],'pi',[],'pes',[],'pg',[], ...
    'pb',[],'ps',[],'pesc',[],'xj2',[],'iter',[],'price_updates',[], ...
    'trace_violation',[],'trace_h',[],'trace_obj',[],'trace_lambda',[]);
val_pi_l_pos_au(:,1) = init_pi_l_pos;
val_pi_l_neg_au(:,1) = init_pi_l_neg;
val_pi_PB_au(1,1) = init_pi_PB;

outer_rng_seed = 999;
rng(outer_rng_seed, 'twister');
n_comm = end_idx - start_idx + 1;
k0_max = 10;
k0 = min(k0_max,ceil(0.5 * sqrt(n_comm)));
delay = zeros(num_LESMs,1);
count = zeros(num_LESMs,1);
for i = 1:num_LESMs
    delay(i) = randi([0, k0(i)]);
    val_xj2(i,1) = sum(init_pes(start_idx(i):end_idx(i)) .^ 2);
end
mask = false(num_LESMs,1);
sync_cnt = 0;
sync_price_step = nan(128,1);

val_pi = init_pi;
val_pes(:,1) = init_pes;
val_pg(:,1) = pg;
val_pb(:,1) = pb;
val_ps(:,1) = ps;
val_pesc(:,1) = init_pesc;
val_xj2(:,1) = xj2;

min_outer_iter = 100;
stable_outer_window = 5;
max_outer_iter = 1000;
errTol_output_CV_MW = 0.01;
outer_avg_start_iter = 1;
outer_history_capacity = min(max_outer_iter + 2, 128);
outer_price_step = nan(outer_history_capacity, 1);
outer_wait_steps = nan(outer_history_capacity, 1);
inner_record_start_iter = min_outer_iter;
record_k = 20:30;
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

%% 上层遍历平均：只在一次完整同步协调完成后计入
sum_outer_pes  = zeros(num_prosumers,1);
sum_outer_pg   = zeros(num_prosumers,1);
sum_outer_pb   = zeros(num_prosumers,1);
sum_outer_ps   = zeros(num_prosumers,1);
sum_outer_pesc = zeros(num_LESMs,1);
outer_avg_cnt = 0;
bar_pes = pes;
bar_pg = pg;
bar_pb = pb;
bar_ps = ps;
bar_pesc = pesc;
val_totcost_iter = nan(outer_history_capacity, 1);
val_outer_violation_iter = nan(outer_history_capacity, 1);

while true
    if k + 1 > outer_history_capacity
        nextCapacity = growHistoryCapacity(outer_history_capacity, k + 1, max_outer_iter + 2);
        sync_price_step(nextCapacity,1) = NaN;
        outer_price_step(nextCapacity,1) = NaN;
        outer_wait_steps(nextCapacity,1) = NaN;
        val_totcost_iter(nextCapacity,1) = NaN;
        val_outer_violation_iter(nextCapacity,1) = NaN;
        outer_history_capacity = nextCapacity;
    end

    beta_k = outer_momentum;
    % Synchronous case: one k is one completed coordination round.  Delays
    % contribute to communication time, but do not create fictitious price steps.
    outer_wait_steps(k) = max(delay) + 1;
    count = delay;
    for i = 1:num_LESMs
        record_inner_or_not = ismember(k,record_k) && (i == record_inner_comm) && (count(i) == delay(i)) && ~over_trace;
        if ~mask(i)
            if count(i) == delay(i)
                local_tic = tic;
                res(i) = localMarket_S2(c(start_idx(i):end_idx(i)), ...
                    b(start_idx(i):end_idx(i)), a(i), ...
                    D(start_idx(i):end_idx(i)), init_rho(i), ...
                    pg_max(start_idx(i):end_idx(i)), pg_min(start_idx(i):end_idx(i)), ...
                    pi(i), ...
                    pes(start_idx(i):end_idx(i)), ...
                    pg(start_idx(i):end_idx(i)), ...
                    pb(start_idx(i):end_idx(i)), ...
                    ps(start_idx(i):end_idx(i)), ...
                    pesc(i), ...
                    errTol_LESMs(i), pi_0(i), pi_max, pi_min, record_inner_or_not, ...
                beta_qre(start_idx(i):end_idx(i)), qre_error_rel, qre_z_cap, inner_momentum, max_inner_iter);

                local_elapsed = toc(local_tic);
            total_inner_iter = total_inner_iter + res(i).iter;
            total_inner_price_updates = total_inner_price_updates + res(i).price_updates;
                num_local_call = num_local_call + 1;
                max_inner_wall_iter_observed = max(max_inner_wall_iter_observed, res(i).iter);
                max_inner_price_updates_observed = max(max_inner_price_updates_observed, res(i).price_updates);
                total_local_time = total_local_time + local_elapsed;

                pi(i) = res(i).pi;
                pes(start_idx(i):end_idx(i)) = res(i).pes;
                pg(start_idx(i):end_idx(i))  = res(i).pg;
                pb(start_idx(i):end_idx(i))  = res(i).pb;
                ps(start_idx(i):end_idx(i))  = res(i).ps;
                pesc(i) = res(i).pesc;
                xj2(i) = res(i).xj2;
                mask(i) = true;

                count(i) = 0;
                delay(i) = randi([0,k0(i)]);

                if record_inner_or_not
                    obj_ref = exactLowerObj( ...
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

                    obj_scale = max(abs(obj_ref), eps);
                    RE_obj_trace = 100 * abs(obj_trace - obj_ref) / obj_scale;

                    inner_k.h = h_trace;
                    inner_k.val_violation = violation_trace;
                    inner_k.val_obj = obj_trace;
                    inner_k.obj_ref = obj_ref;
                    inner_k.RE_obj = RE_obj_trace;
                    inner_k.record_outer_iter = k;
                    inner_k.record_comm = i;

                    record_idx = numel(h_trace);
                    record_h = h_trace(record_idx);
                    inner_k.final_RE_obj = RE_obj_trace(record_idx);
                    inner_k.final_val_violation = violation_trace(record_idx);
                    inner_k.record_h = record_h;
                    inner_k.record_h_used = h_trace(record_idx);
                    inner_k.record_index = record_idx;
                    inner_k.record_source = 'own_final_h';
                    inner_k.record_extrapolated_by_hold = false;
                    inner_k.record_RE_obj = inner_k.final_RE_obj;
                    inner_k.record_val_violation = inner_k.final_val_violation;
                    inner_k.own_final_RE_obj = inner_k.final_RE_obj;
                    inner_k.own_final_val_violation = inner_k.final_val_violation;

                    over_trace = 1;
                end
            else
                count(i) = count(i) + 1;
            end
        end

        val_pi(i,k+1) = pi(i);
        val_pes(start_idx(i):end_idx(i),k+1) = pes(start_idx(i):end_idx(i));
        val_pg(start_idx(i):end_idx(i),k+1) = pg(start_idx(i):end_idx(i));
        val_pb(start_idx(i):end_idx(i),k+1) = pb(start_idx(i):end_idx(i));
        val_ps(start_idx(i):end_idx(i),k+1) = ps(start_idx(i):end_idx(i));
        val_pesc(i,k+1) = pesc(i);
        val_xj2(i,k+1) = xj2(i);
    end

    if all(mask)

        if k >= outer_avg_start_iter
            outer_avg_cnt = outer_avg_cnt + 1;
            sum_outer_pes  = sum_outer_pes  + pes;
            sum_outer_pg   = sum_outer_pg   + pg;
            sum_outer_pb   = sum_outer_pb   + pb;
            sum_outer_ps   = sum_outer_ps   + ps;
            sum_outer_pesc = sum_outer_pesc + pesc;

            bar_pes  = sum_outer_pes  / outer_avg_cnt;
            bar_pg   = sum_outer_pg   / outer_avg_cnt;
            bar_pb   = sum_outer_pb   / outer_avg_cnt;
            bar_ps   = sum_outer_ps   / outer_avg_cnt;
            bar_pesc = sum_outer_pesc / outer_avg_cnt;
        else
            bar_pes = pes;
            bar_pg = pg;
            bar_pb = pb;
            bar_ps = ps;
            bar_pesc = pesc;
        end

        val_pi_PB_au(1,sync_cnt+2) = val_pi_PB(1,sync_cnt+1) - alpha_PB * sum(pesc);
        momentum_pi_PB = beta_k * (val_pi_PB_au(1,sync_cnt+2) - val_pi_PB_au(1,sync_cnt+1));
        val_pi_PB(1,sync_cnt+2) = val_pi_PB_au(1,sync_cnt+2) + momentum_pi_PB;

        raw_pi_l_pos = val_pi_l_pos(:,sync_cnt+1) - alpha_l * (pai_il' * pesc - F_l);
        raw_pi_l_neg = val_pi_l_neg(:,sync_cnt+1) - alpha_l * (-pai_il' * pesc - F_l);
        val_pi_l_pos_au(:,sync_cnt+2) = min(raw_pi_l_pos, 0);
        val_pi_l_neg_au(:,sync_cnt+2) = min(raw_pi_l_neg, 0);
        momentum_pi_l_pos = beta_k * (val_pi_l_pos_au(:,sync_cnt+2) - val_pi_l_pos_au(:,sync_cnt+1));
        momentum_pi_l_neg = beta_k * (val_pi_l_neg_au(:,sync_cnt+2) - val_pi_l_neg_au(:,sync_cnt+1));
        val_pi_l_pos(:,sync_cnt+2) = val_pi_l_pos_au(:,sync_cnt+2) + momentum_pi_l_pos;
        val_pi_l_neg(:,sync_cnt+2) = val_pi_l_neg_au(:,sync_cnt+2) + momentum_pi_l_neg;

        old_pi_0 = pi_0;
        val_pi_0(:,sync_cnt+2) = val_pi_PB(1,sync_cnt+2) + pai_il * (val_pi_l_pos(:,sync_cnt+2) - val_pi_l_neg(:,sync_cnt+2));
        pi_0 = val_pi_0(:,sync_cnt+2);

        sync_cnt = sync_cnt + 1;
        sync_price_step(sync_cnt) = norm(pi_0 - old_pi_0, inf);
        outer_price_step(k) = sync_price_step(sync_cnt);
        mask(:) = false;
    end

    val_totcost_iter(k) = sum(0.5 .* c .* bar_pg.^2 ...
        + b .* bar_pg ...
        + pi_max .* bar_pb ...
        - pi_min .* bar_ps ...
        + 0.5 .* a_extend .* bar_pes.^2) ...
        + sum(0.5 .* a .* bar_pesc.^2);%由于内层传上来的就已经是内层残差为0的结果了，所以就这样

    F_l_MW_bar = F_l / 1e3;
    y_bar_now_MW = bar_pesc / 1e3;
    PB_bar_now_MW = abs(sum(y_bar_now_MW));
    LINE_pos_bar_now_MW = max(pai_il' * y_bar_now_MW - F_l_MW_bar, 0);
    LINE_neg_bar_now_MW = max(-pai_il' * y_bar_now_MW - F_l_MW_bar, 0);
    val_outer_violation_iter(k) = norm([PB_bar_now_MW; LINE_pos_bar_now_MW; LINE_neg_bar_now_MW], 2);

    if k >= min_outer_iter && sync_cnt >= stable_outer_window
        recent_step = sync_price_step(sync_cnt-stable_outer_window+1:sync_cnt);
        if all(recent_step <= errTol_UESM) && ...
                val_outer_violation_iter(k) <= errTol_output_CV_MW
            break;
        end
    end

    if k >= max_outer_iter
        state.success = 0;
        break;
    end

    k = k + 1;
end

output_pes  = bar_pes;
output_pg   = bar_pg;
output_pb   = bar_pb;
output_ps   = bar_ps;
output_pesc = bar_pesc;

final_gap_kW = zeros(num_LESMs,1);
for i = 1:num_LESMs
    final_gap_kW(i) = abs(output_pesc(i) - sum(output_pes(start_idx(i):end_idx(i))));
end

outer = struct('k',[], ...
    'iter',[], ...
    'val_violation',[], ...
    'totcost',[], ...
    'totcost_ref',[], ...
    'RE_totcost',[], ...
    'pi',[], ...
    'pesc',[], ...
    'pi_0',[], ...
    'final_gap_MW',[]);
outer.final_gap_MW = final_gap_kW / 1e3;

totcost_ref = exactUpperObj(a,c,b,pi_max,pi_min,D,F_l,pai_il, ...
    start_idx,end_idx,pg_max,pg_min,a_extend);
totcost_scale = max(abs(totcost_ref), eps);

outer_iter = (1:k)';
totcost_iter = val_totcost_iter(1:k);
violation_iter = val_outer_violation_iter(1:k);
RE_totcost = 100 * abs(totcost_iter - totcost_ref) / totcost_scale;

record_idx_outer = numel(RE_totcost);
record_iter = outer_iter(record_idx_outer);
record_extrapolated_outer_by_hold = false;
record_RE_totcost = RE_totcost(record_idx_outer);
record_val_violation = violation_iter(record_idx_outer);
own_final_RE_totcost = record_RE_totcost;

F_l_MW = F_l / 1e3;
y_out_MW = output_pesc / 1e3;
PB_MW = abs(sum(y_out_MW));
LINE_pos_MW = max(pai_il' * y_out_MW - F_l_MW, 0);
LINE_neg_MW = max(-pai_il' * y_out_MW - F_l_MW, 0);
CV_MW = norm([PB_MW; LINE_pos_MW; LINE_neg_MW], 2);
own_final_val_violation = CV_MW;
line_MW = max([0; LINE_pos_MW; LINE_neg_MW]);

outer.k = outer_iter;
outer.iter = outer_iter;
outer.val_violation = violation_iter;
outer.totcost = totcost_iter;
outer.totcost_ref = totcost_ref;
outer.RE_totcost = RE_totcost;
outer.final_RE_totcost = record_RE_totcost;
outer.final_val_violation = CV_MW;
outer.record_iter = record_iter;
outer.record_iter_used = record_idx_outer;
outer.record_source = 'own_final_iter';
outer.record_extrapolated_by_hold = record_extrapolated_outer_by_hold;
outer.record_RE_totcost = record_RE_totcost;
outer.record_val_violation = record_val_violation;
outer.own_final_RE_totcost = own_final_RE_totcost;
outer.own_final_val_violation = own_final_val_violation;
outer.CV_MW = CV_MW;
outer.PB_MW = PB_MW;
outer.line_MW = line_MW;
outer.final_gap_mean_MW = mean(outer.final_gap_MW,'omitnan');
outer.final_gap_max_MW = max(outer.final_gap_MW);
outer.final_gap_p95_MW = prctile(outer.final_gap_MW,95);
outer.outer_avg_cnt = outer_avg_cnt;
outer.sync_cnt = sync_cnt;
outer.outer_price_step = outer_price_step(1:k);
outer.sync_price_step = sync_price_step(1:sync_cnt);
outer.outer_wait_steps = outer_wait_steps(1:k);
outer.total_outer_wait_steps = sum(outer.outer_wait_steps);
outer.min_outer_iter = min_outer_iter;
outer.stable_outer_window = stable_outer_window;
outer.max_outer_iter = max_outer_iter;
outer.errTol_UESM = errTol_UESM;
outer.errTol_output_CV_MW = errTol_output_CV_MW;
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
outer.pi = pi(:);
outer.pi_0 = pi_0(:);
outer.pes = output_pes;
outer.pg = output_pg;
outer.pb = output_pb;
outer.ps = output_ps;
outer.pesc = output_pesc / 1e3;

save('inner_k_data_2.mat','inner_k');
save('outer_data_2.mat','outer');
toc;
