function summary = run_outer_units_audit()
%RUN_OUTER_UNITS_AUDIT  Audit the formal outer scales and finite-line mask.

P = load('params.mat');
finite_line_mask = isfinite(P.F_l);
summary = struct();
summary.power_unit = 'kW';
summary.line_limit_unit = 'kW';
summary.num_communities = P.num_LESMs;
summary.num_finite_lines = nnz(finite_line_mask);
summary.finite_line_indices = find(finite_line_mask);
summary.min_a = min(P.a);
summary.max_a = max(P.a);
summary.min_c_out = min(P.c_out);
summary.max_c_out = max(P.c_out);
summary.LD_out = P.LD_out;
summary.alpha_out_limit = P.alpha_out_limit;
summary.alpha_PB = P.alpha_PB;
summary.alpha_l = P.alpha_l;
summary.alpha_ratio = P.alpha_PB / P.alpha_out_limit;
summary.min_pi_0 = min(P.init_pi_0);
summary.max_pi_0 = max(P.init_pi_0);
summary.min_pesc_init = min(P.init_pesc);
summary.max_pesc_init = max(P.init_pesc);
summary.sum_abs_pesc_init = sum(abs(P.init_pesc));
summary.system_balance_scale_kW = P.system_balance_scale_kW;
summary.outer_pb_cv_tol_kW = P.outer_pb_cv_tol_kW;
summary.outer_line_cv_tol_kW = P.outer_line_cv_tol_kW(finite_line_mask);
summary.capacity_min_kW = min(P.F_l(finite_line_mask));
summary.capacity_max_kW = max(P.F_l(finite_line_mask));
summary.formal_parameters = struct('run_profile',P.run_profile, ...
    'qre_noise_enabled',P.qre_noise_enabled, ...
    'qre_certificate_enabled',P.qre_certificate_enabled, ...
    'outer_momentum',P.outer_momentum, ...
    'outer_avg_start_price_update',P.outer_avg_start_price_update);
summary.success_gate = strcmp(summary.power_unit,'kW') && ...
    strcmp(summary.line_limit_unit,'kW') && summary.num_finite_lines == 7 && ...
    summary.alpha_PB > 0 && summary.alpha_l > 0 && ...
    abs(summary.alpha_PB-summary.alpha_l) <= eps(max(1,abs(summary.alpha_PB)));
save('outer_units_audit.mat','summary');
fprintf('Outer units audit: %s; finite lines=%d; LD_out=%.8e; alpha=%.8e\n', ...
    ternary(summary.success_gate,'PASS','FAIL'),summary.num_finite_lines, ...
    summary.LD_out,summary.alpha_PB);

    function out = ternary(condition,yes_value,no_value)
        if condition, out = yes_value; else, out = no_value; end
    end
end
