function result = solveExactCommunityLambda( ...
    c,b,a,D,pg_min,pg_max,pi_0,pi_max,pi_min,lambda_initial)
%SOLVEEXACTCOMMUNITYLAMBDA Solve the exact fixed-pi0 community balance root.
% This diagnostic solves only the scalar exact community residual. It does
% not alter any held state or any formal market update.

N = numel(c);
lambda_initial = double(lambda_initial);
root_tol_lambda = 1e-12;
root_tol_residual_kW = 1e-8;
root_max_iter = 200;
expand_factor = 2;
max_expand_iter = 50;

lambda_left = lambda_initial - 0.1;
lambda_right = lambda_initial + 0.1;
res_left = exactCommunityResidual(lambda_left);
res_right = exactCommunityResidual(lambda_right);
expand_iter = 0;

while res_left * res_right > 0 && expand_iter < max_expand_iter
    half_width = 0.1 * expand_factor^(expand_iter + 1);
    lambda_left = lambda_initial - half_width;
    lambda_right = lambda_initial + half_width;
    res_left = exactCommunityResidual(lambda_left);
    res_right = exactCommunityResidual(lambda_right);
    expand_iter = expand_iter + 1;
end

result = struct();
result.lambda_exact = NaN;
result.residual_exact = NaN;
result.iter = 0;
result.success = false;
result.lower_bound = lambda_left;
result.upper_bound = lambda_right;
result.bracket_expansions = expand_iter;

if res_left * res_right > 0
    result.failure_reason = 'root_bracket_not_found';
    return;
end

for root_iter = 1:root_max_iter
    lambda_mid = 0.5 * (lambda_left + lambda_right);
    res_mid = exactCommunityResidual(lambda_mid);
    result.iter = root_iter;
    if abs(res_mid) <= root_tol_residual_kW || ...
            abs(lambda_right - lambda_left) <= root_tol_lambda
        result.lambda_exact = lambda_mid;
        result.residual_exact = res_mid;
        result.success = true;
        result.lower_bound = lambda_left;
        result.upper_bound = lambda_right;
        result.failure_reason = '';
        return;
    end
    if res_left * res_mid <= 0
        lambda_right = lambda_mid;
        res_right = res_mid;
    else
        lambda_left = lambda_mid;
        res_left = res_mid;
    end
end

result.lambda_exact = 0.5 * (lambda_left + lambda_right);
result.residual_exact = exactCommunityResidual(result.lambda_exact);
result.success = abs(result.residual_exact) <= root_tol_residual_kW;
result.lower_bound = lambda_left;
result.upper_bound = lambda_right;
if result.success
    result.failure_reason = '';
else
    result.failure_reason = 'root_max_iter';
end

    function residual = exactCommunityResidual(lambda_value)
        pes_exact = zeros(N,1);
        for j = 1:N
            sol = solveExactEUResponse( ...
                c(j),b(j),a,D(j),pg_min(j),pg_max(j), ...
                lambda_value,pi_max,pi_min);
            pes_exact(j) = sol.pes;
        end
        pesc_exact = (pi_0 - lambda_value) / a;
        residual = pesc_exact - sum(pes_exact);
    end
end
