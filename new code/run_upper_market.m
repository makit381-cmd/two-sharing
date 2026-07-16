function run_upper_market(scenario)
% Unified outer loop. Scenario 0/2 waits for a full community response cycle;
% scenario 1/3 updates from held asynchronous community outputs.
root_dir = fileparts(mfilename('fullpath'));
data = load(fullfile(root_dir, 'param.mat'));
is_synchronous = any(scenario == [0,2]);
solver = str2func(sprintf('localMarketS%d', scenario));

a = data.a; b = data.b; c = data.c; D = data.D;
start_idx = data.start_idx; end_idx = data.end_idx;
pg_max = data.pg_max; pg_min = data.pg_min; beta_qre = data.beta_qre;
pi_max = data.pi_max; pi_min = data.pi_min; F_l = data.F_l;
PTDF = data.PTDF; alpha_PB = data.alpha_PB; alpha_l = data.alpha_l;
init_rho = data.init_rho; errTol_LESMs = data.errTol_LESMs;
errTol_UESM = data.errTol_UESM; a_extend = data.a_extend;
num_LESMs = data.num_LESMs; num_prosumers = data.num_prosumers;

pg = data.init_pg; pes = data.init_pes; pb = data.init_pb; ps = data.init_ps;
pesc = data.init_pesc; pi = data.init_pi; pi_0 = data.init_pi_0;
pi_PB = data.init_pi_PB; pi_l_pos = data.init_pi_l_pos; pi_l_neg = data.init_pi_l_neg;

min_outer_iter = 100;
max_outer_iter = 1500;
stable_outer_window = 5;
errTol_output_CV_MW = 2e-1;
outer_avg_start_iter = 1;
outer_rng_seed = 999;
k0 = min(10, ceil(0.5 * sqrt(end_idx - start_idx + 1)));
rng(outer_rng_seed, 'twister');
delay = arrayfun(@(x) randi([0,x]), k0);
count = zeros(num_LESMs,1);
response_mask = false(num_LESMs,1);

sum_pes = zeros(num_prosumers,1);
sum_pg = zeros(num_prosumers,1);
sum_pb = zeros(num_prosumers,1);
sum_ps = zeros(num_prosumers,1);
sum_pesc = zeros(num_LESMs,1);
outer_avg_cnt = 0;
price_steps = nan(max_outer_iter,1);
sync_steps = nan(max_outer_iter,1);
sync_cnt = 0;
totcost = nan(max_outer_iter,1);
violation = nan(max_outer_iter,1);
inner_k = struct('h',[],'val_violation',[],'val_obj',[],'obj_ref',[],'RE_obj',[]);
recorded_inner = false;
total_inner_iter = 0;
num_local_call = 0;
max_inner_iter = 0;
total_local_time = 0;
state.success = 1;

for k = 1:max_outer_iter
    for i = 1:num_LESMs
        idx = start_idx(i):end_idx(i);
        local_ready = count(i) >= delay(i);
        if local_ready
            record_trace = ~recorded_inner && i == 71 && k >= 100;
            local_tic = tic;
            res = solver(c(idx),b(idx),a(i),D(idx),init_rho(i), ...
                pg_max(idx),pg_min(idx),pi(i),pes(idx),pg(idx),pb(idx),ps(idx), ...
                pesc(i),errTol_LESMs(i),pi_0(i),pi_max,pi_min,record_trace,beta_qre(idx));
            total_local_time = total_local_time + toc(local_tic);
            total_inner_iter = total_inner_iter + res.iter;
            num_local_call = num_local_call + 1;
            max_inner_iter = max(max_inner_iter,res.iter);
            pi(i) = res.pi;
            pes(idx) = res.pes;
            pg(idx) = res.pg;
            pb(idx) = res.pb;
            ps(idx) = res.ps;
            pesc(i) = res.pesc;
            count(i) = 0;
            delay(i) = randi([0,k0(i)]);
            response_mask(i) = true;
            if record_trace
                ref = exactLowerObj(c(idx),b(idx),a(i),D(idx),pg_max(idx),pg_min(idx),pi_0(i),pi_max,pi_min);
                inner_k.h = res.trace_h;
                inner_k.val_violation = res.trace_violation;
                inner_k.val_obj = res.trace_obj;
                inner_k.obj_ref = ref;
                inner_k.RE_obj = 100 * abs(res.trace_obj - ref) / max(abs(ref),eps);
                recorded_inner = true;
            end
        else
            count(i) = count(i) + 1;
        end
    end

    % Every scenario reports a time-slot ergodic average from its first outer slot.
    if k >= outer_avg_start_iter
        outer_avg_cnt = outer_avg_cnt + 1;
        sum_pes = sum_pes + pes;
        sum_pg = sum_pg + pg;
        sum_pb = sum_pb + pb;
        sum_ps = sum_ps + ps;
        sum_pesc = sum_pesc + pesc;
    end
    bar_pes = sum_pes / outer_avg_cnt;
    bar_pg = sum_pg / outer_avg_cnt;
    bar_pb = sum_pb / outer_avg_cnt;
    bar_ps = sum_ps / outer_avg_cnt;
    bar_pesc = sum_pesc / outer_avg_cnt;

    should_update_price = ~is_synchronous || all(response_mask);
    if should_update_price
        old_pi_0 = pi_0;
        pi_PB = pi_PB - alpha_PB * sum(pesc);
        pi_l_pos = min(pi_l_pos - alpha_l * (PTDF' * pesc - F_l),0);
        pi_l_neg = min(pi_l_neg - alpha_l * (-PTDF' * pesc - F_l),0);
        pi_0 = pi_PB + PTDF * (pi_l_pos - pi_l_neg);
        step = norm(pi_0 - old_pi_0,inf);
        price_steps(k) = step;
        if is_synchronous
            sync_cnt = sync_cnt + 1;
            sync_steps(sync_cnt) = step;
            response_mask(:) = false;
        end
    end

    totcost(k) = sum(0.5 .* c .* bar_pg.^2 + b .* bar_pg + pi_max .* bar_pb ...
        - pi_min .* bar_ps + 0.5 .* a_extend .* bar_pes.^2) + sum(0.5 .* a .* bar_pesc.^2);
    y_MW = bar_pesc / 1e3;
    line_pos = max(PTDF' * y_MW - F_l / 1e3,0);
    line_neg = max(-PTDF' * y_MW - F_l / 1e3,0);
    violation(k) = norm([abs(sum(y_MW)); line_pos; line_neg],2);

    if is_synchronous
        ready_to_stop = sync_cnt >= stable_outer_window && ...
            all(sync_steps(sync_cnt-stable_outer_window+1:sync_cnt) <= errTol_UESM);
    else
        ready_to_stop = k >= stable_outer_window && ...
            all(price_steps(k-stable_outer_window+1:k) <= errTol_UESM);
    end
    if k >= min_outer_iter && ready_to_stop && violation(k) <= errTol_output_CV_MW
        break;
    end
    if k == max_outer_iter
        state.success = 0;
    end
end

totcost_ref = exactUpperObj(a,c,b,pi_max,pi_min,D,F_l,PTDF,start_idx,end_idx,pg_max,pg_min,a_extend);
RE_totcost = 100 * abs(totcost(1:k) - totcost_ref) / max(abs(totcost_ref),eps);
final_gap = zeros(num_LESMs,1);
for i = 1:num_LESMs
    final_gap(i) = abs(bar_pesc(i) - sum(bar_pes(start_idx(i):end_idx(i)))) / 1e3;
end

outer = struct();
outer.k = (1:k)';
outer.iter = outer.k;
outer.totcost = totcost(1:k);
outer.totcost_ref = totcost_ref;
outer.RE_totcost = RE_totcost;
outer.val_violation = violation(1:k);
outer.final_RE_totcost = RE_totcost(end);
outer.final_val_violation = violation(k);
outer.CV_MW = violation(k);
outer.PB_MW = abs(sum(bar_pesc / 1e3));
outer.line_MW = max([0; line_pos; line_neg]);
outer.final_gap_MW = final_gap;
outer.final_gap_mean_MW = mean(final_gap);
outer.final_gap_max_MW = max(final_gap);
outer.final_gap_p95_MW = prctile(final_gap,95);
outer.outer_avg_start_iter = outer_avg_start_iter;
outer.outer_avg_cnt = outer_avg_cnt;
outer.min_outer_iter = min_outer_iter;
outer.max_outer_iter = max_outer_iter;
outer.errTol_UESM = errTol_UESM;
outer.errTol_output_CV_MW = errTol_output_CV_MW;
outer.outer_price_step = price_steps(1:k);
outer.sync_cnt = sync_cnt;
outer.random_seed = outer_rng_seed;
outer.state = state;
outer.initial_source = data.formal_metadata.initial_source;
outer.parameter_b_range = data.formal_metadata.b_range;
outer.parameter_c_range = data.formal_metadata.c_range;
outer.num_local_call = num_local_call;
outer.total_inner_iter = total_inner_iter;
outer.max_inner_iter = max_inner_iter;
outer.total_local_time = total_local_time;
outer.avg_inner_iter = total_inner_iter / max(num_local_call,1);
outer.pi = pi;
outer.pi_0 = pi_0;
outer.pes = bar_pes;
outer.pg = bar_pg;
outer.pb = bar_pb;
outer.ps = bar_ps;
outer.pesc = bar_pesc / 1e3;

save(fullfile(root_dir, sprintf('inner_k_data_%d.mat',scenario)), 'inner_k');
save(fullfile(root_dir, sprintf('outer_data_%d.mat',scenario)), 'outer');
fprintf('S%d finished at k=%d, RE=%.6g%%, CV=%.6g MW.\n',scenario,k,outer.final_RE_totcost,outer.CV_MW);
end
