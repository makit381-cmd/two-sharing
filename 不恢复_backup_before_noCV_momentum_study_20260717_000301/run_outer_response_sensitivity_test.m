function summary = run_outer_response_sensitivity_test()
%RUN_OUTER_RESPONSE_SENSITIVITY_TEST  Fixed-state finite differences of raw Agg.
% The local analytic/QRE response is used; no full outer solver is called.

P = load('params.mat');
comm_list = [1 3 49 71 123];
delta_list = [1e-6 1e-5 1e-4 1e-3];
num_comm = numel(comm_list);
num_delta = numel(delta_list);

formal_parameters = struct('qre_noise_enabled',P.qre_noise_enabled, ...
    'qre_audit_enabled',P.qre_audit_enabled,'inner_avg_policy',P.inner_avg_policy, ...
    'inner_avg_policy_version',P.inner_avg_policy_version,'max_inner_iter',P.max_inner_iter, ...
    'inner_cv_stop_enabled',P.inner_cv_stop_enabled,'h0_max',P.h0_max);
test_overrides = struct('qre_noise_enabled',false,'qre_audit_enabled',false, ...
    'fixed_warm_start',true,'fixed_local_call_count',1, ...
    'label','finite-difference diagnostic only; not formal profile');

rng_config = P.rng_seeds;
avg_config = struct('policy',P.inner_avg_policy, ...
    'policy_version',P.inner_avg_policy_version, ...
    'start_cv_factor',P.inner_avg_start_cv_factor, ...
    'start_price_factor',P.inner_avg_start_price_factor, ...
    'start_stable_window',P.inner_avg_start_stable_window, ...
    'min_price_updates',P.inner_avg_min_price_updates, ...
    'min_samples',P.inner_avg_min_samples, ...
    'formal_cv_stable_window',P.inner_formal_cv_stable_window);

records = repmat(struct('community_id',NaN,'delta',NaN,'pi0',NaN, ...
    'pesc_raw',NaN,'pes_sum',NaN,'inner_cv_kW',NaN,'lambda_final',NaN, ...
    'local_iterations',NaN,'price_updates',NaN,'wall_iter',NaN,'runtime_s',NaN, ...
    'stop_reason',''),num_comm*num_delta*2,1);
row = 0;
for ci = 1:num_comm
    i = comm_list(ci);
    for di = 1:num_delta
        delta = delta_list(di);
        pi0_values = [P.init_pi_0(i)-delta, P.init_pi_0(i)+delta];
        for si = 1:2
            row = row + 1;
            t = tic;
            r = localMarket_S3(P.c(P.start_idx(i):P.end_idx(i)), ...
                P.b(P.start_idx(i):P.end_idx(i)),P.a(i), ...
                P.D(P.start_idx(i):P.end_idx(i)),P.init_rho(i), ...
                P.pg_max(P.start_idx(i):P.end_idx(i)),P.pg_min(P.start_idx(i):P.end_idx(i)), ...
                P.init_pi(i),P.init_pes(P.start_idx(i):P.end_idx(i)), ...
                P.init_pg(P.start_idx(i):P.end_idx(i)),P.init_pb(P.start_idx(i):P.end_idx(i)), ...
                P.init_ps(P.start_idx(i):P.end_idx(i)),P.init_pesc(i),P.errTol_LESMs(i), ...
                P.inner_cv_tol_kW(i),pi0_values(si),P.pi_max,P.pi_min,false, ...
                P.beta_qre(P.start_idx(i):P.end_idx(i)),P.qre_epsilon,P.agg_epsilon_i(i), ...
                P.qre_z_cap,P.qre_backoff_factor,P.qre_max_backoffs,P.inner_momentum, ...
                P.max_inner_iter,P.h0_max,P.min_inner_iter,P.inner_avg_start_iter, ...
                P.inner_cv_stop_enabled,P.qre_certificate_enabled,P.agg_certificate_enabled, ...
                P.agg_cert_tol,P.stable_inner_window,false,0,false,P.qre_audit_seed, ...
                P.agg_gap_diagnostic_enabled,false,P.diagnostic_record_every, ...
                P.exact_sync_diagnostic_every,P.rolling_window,false,inf,1,i,rng_config, ...
                avg_config,1);
            records(row).community_id = i;
            records(row).delta = delta;
            records(row).pi0 = pi0_values(si);
            records(row).pesc_raw = r.pesc_agg_raw;
            records(row).pes_sum = r.pes_sum;
            records(row).inner_cv_kW = r.inner_cv_kW;
            records(row).lambda_final = r.pi;
            records(row).local_iterations = r.iter;
            records(row).price_updates = r.price_updates;
            records(row).wall_iter = r.wall_iter;
            records(row).runtime_s = toc(t);
            records(row).stop_reason = r.stop_reason;
        end
    end
end

derivative = repmat(struct('community_id',NaN,'delta',NaN,'d_pesc_d_pi0',NaN, ...
    'one_over_a',NaN,'one_over_c_out',NaN,'ratio_to_one_over_a',NaN, ...
    'ratio_to_one_over_c_out',NaN,'plus_cv_kW',NaN,'minus_cv_kW',NaN), ...
    num_comm*num_delta,1);
row = 0;
for ci = 1:num_comm
    i = comm_list(ci);
    for di = 1:num_delta
        row = row + 1;
        idx0 = ((ci-1)*num_delta + di)*2 - 1;
        rminus = records(idx0);
        rplus = records(idx0+1);
        d = (rplus.pesc_raw-rminus.pesc_raw)/(2*delta_list(di));
        derivative(row).community_id = i;
        derivative(row).delta = delta_list(di);
        derivative(row).d_pesc_d_pi0 = d;
        derivative(row).one_over_a = 1/P.a(i);
        derivative(row).one_over_c_out = 1/P.c_out(i);
        derivative(row).ratio_to_one_over_a = d/(1/P.a(i));
        derivative(row).ratio_to_one_over_c_out = d/(1/P.c_out(i));
        derivative(row).plus_cv_kW = rplus.inner_cv_kW;
        derivative(row).minus_cv_kW = rminus.inner_cv_kW;
    end
end

finite_d = [derivative.d_pesc_d_pi0]';
summary = struct();
summary.community_list = comm_list;
summary.delta_list = delta_list;
summary.formal_parameters = formal_parameters;
summary.test_overrides = test_overrides;
summary.records = records;
summary.derivative = derivative;
summary.derivative_min = min(finite_d);
summary.derivative_max = max(finite_d);
summary.theory_one_over_a = 1./P.a(comm_list(:));
summary.theory_one_over_c_out = 1./P.c_out(comm_list(:));
summary.success_gate = all(isfinite(finite_d)) && all([records.inner_cv_kW] >= 0) && ...
    all(isfinite([records.pesc_raw]));
save('outer_response_sensitivity_test.mat','summary','-v7.3');
fprintf('Outer response sensitivity test: %s; derivative range=[%.6e, %.6e]\n', ...
    ternary(summary.success_gate,'PASS','FAIL'),summary.derivative_min,summary.derivative_max);

    function out = ternary(condition,yes_value,no_value)
        if condition, out = yes_value; else, out = no_value; end
    end
end
