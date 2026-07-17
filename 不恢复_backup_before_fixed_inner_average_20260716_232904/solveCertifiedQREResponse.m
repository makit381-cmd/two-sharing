function sol = solveCertifiedQREResponse( ...
    c_j,b_j,a_i,D_j,pg_min_j,pg_max_j,lambda,pi_max,pi_min,beta_j, ...
    qre_epsilon,qre_z_cap,qre_backoff_factor,qre_max_backoffs,qre_stream)
%SOLVECERTIFIEDQRERESPONSE Generate and certify one QRE response.
% One bounded random direction is reused for all geometric backoffs.

if nargin < 15 || isempty(qre_stream)
    qre_stream = RandStream('mt19937ar','Seed',31002);
end

validateattributes(beta_j, {'numeric'}, {'scalar','real','finite','positive'});
validateattributes(qre_epsilon, {'numeric'}, {'scalar','real','finite','nonnegative'});
validateattributes(qre_z_cap, {'numeric'}, {'scalar','real','finite','positive'});
validateattributes(qre_backoff_factor, {'numeric'}, {'scalar','real','finite','positive','<',1});
validateattributes(qre_max_backoffs, {'numeric'}, {'scalar','real','finite','integer','nonnegative'});

exact = solveExactEUResponse( ...
    c_j,b_j,a_i,D_j,pg_min_j,pg_max_j,lambda,pi_max,pi_min);

z_pg_raw = max(min(randn(qre_stream,1), qre_z_cap), -qre_z_cap);
z_pes_raw = max(min(randn(qre_stream,1), qre_z_cap), -qre_z_cap);
cert_tol = 1e-10 * max(1, abs(exact.obj));

gamma = 1;
accepted = false;
backoff_count = 0;
qre = exact;
qre_gap = 0;

for attempt = 0:qre_max_backoffs
    qre = generate_qre(gamma, z_pg_raw, z_pes_raw);
    gap_raw = qre.obj - exact.obj;
    if gap_raw < -cert_tol
        error('solveCertifiedQREResponse:objectiveInconsistency', ...
            'QRE objective is below the exact objective beyond tolerance.');
    end
    qre_gap = max(gap_raw, 0);
    if qre_gap <= qre_epsilon + cert_tol
        accepted = true;
        backoff_count = attempt;
        break;
    end
    if attempt < qre_max_backoffs
        gamma = qre_backoff_factor * gamma;
    end
end

fallback_to_exact = ~accepted;
if fallback_to_exact
    qre = exact;
    gamma = 0;
    qre_gap = 0;
    backoff_count = qre_max_backoffs;
end

sol = qre;
sol.exact_obj = exact.obj;
sol.epsilon_gap = qre_gap;
sol.gamma = gamma;
sol.backoff_count = backoff_count;
sol.fallback_to_exact = fallback_to_exact;
sol.branch_qre = qre.branch;
sol.branch_exact = exact.branch;
sol.branch_switched = qre.branch ~= exact.branch;
sol.z_pg_raw = z_pg_raw;
sol.z_pes_raw = z_pes_raw;
sol.certificate_tolerance = cert_tol;
sol.certificate_pass = qre_gap <= qre_epsilon + cert_tol;
sol.response_source = ternary(fallback_to_exact, 'exact_fallback', 'qre');

    function out = generate_qre(gamma_now, z_pg, z_pes)
        tol_F = 1e-10;
        pg_buyer = min(max((pi_max - b_j) / c_j + ...
            gamma_now * sqrt(1 / (beta_j * c_j)) * z_pg, ...
            pg_min_j), pg_max_j);
        pes_buyer = (lambda - pi_max) / a_i + ...
            gamma_now * sqrt(1 / (beta_j * a_i)) * z_pes;
        F_buyer = D_j + pes_buyer - pg_buyer;

        pg_seller = min(max((pi_min - b_j) / c_j + ...
            gamma_now * sqrt(1 / (beta_j * c_j)) * z_pg, ...
            pg_min_j), pg_max_j);
        pes_seller = (lambda - pi_min) / a_i + ...
            gamma_now * sqrt(1 / (beta_j * a_i)) * z_pes;
        F_seller = D_j + pes_seller - pg_seller;

        if F_buyer > tol_F
            out.pg = pg_buyer;
            out.pes = pes_buyer;
            out.pb = F_buyer;
            out.ps = 0;
            out.branch = 1;
        elseif F_seller < -tol_F
            out.pg = pg_seller;
            out.pes = pes_seller;
            out.pb = 0;
            out.ps = -F_seller;
            out.branch = 2;
        else
            z_middle = (sqrt(c_j) * z_pg + sqrt(a_i) * z_pes) / ...
                sqrt(c_j + a_i);
            pes_middle = (lambda - b_j - c_j * D_j) / (c_j + a_i) + ...
                gamma_now * sqrt(1 / (beta_j * (c_j + a_i))) * z_middle;
            out.pes = min(max(pes_middle, pg_min_j - D_j), pg_max_j - D_j);
            out.pg = D_j + out.pes;
            out.pb = 0;
            out.ps = 0;
            out.branch = 3;
        end
        out.obj = euLocalLagrangian( ...
            out.pg, out.pes, out.pb, out.ps, ...
            c_j, b_j, a_i, lambda, pi_max, pi_min);
        out.balance_residual = D_j + out.pes + out.ps - out.pg - out.pb;
    end

    function value = ternary(condition, true_value, false_value)
        if condition
            value = true_value;
        else
            value = false_value;
        end
    end
end
