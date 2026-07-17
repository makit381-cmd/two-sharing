function [obj_star, meta] = exactLowerObjFast(c,b,a,D,pg_max,pg_min,pi_0,pi_max,pi_min)
% 精确社区参考：以一维 lambda 二分求解当前闭式分支模型，避免每次社区调用 Gurobi。
% 与 exactLowerObj 的目标和约束一致，返回固定 pi_0 下的精确目标及求解状态。

N = numel(c);
lambda_lo = -1;
lambda_hi = 1;
[~, pes_lo] = eu_response(lambda_lo);
[~, pes_hi] = eu_response(lambda_hi);
r_lo = (pi_0 - lambda_lo) / a - sum(pes_lo);
r_hi = (pi_0 - lambda_hi) / a - sum(pes_hi);

expand_iter = 0;
while ~(r_lo >= 0 && r_hi <= 0) && expand_iter < 80
    lambda_lo = 2 * lambda_lo - 1;
    lambda_hi = 2 * lambda_hi + 1;
    [~, pes_lo] = eu_response(lambda_lo);
    [~, pes_hi] = eu_response(lambda_hi);
    r_lo = (pi_0 - lambda_lo) / a - sum(pes_lo);
    r_hi = (pi_0 - lambda_hi) / a - sum(pes_hi);
    expand_iter = expand_iter + 1;
end
if ~(r_lo >= 0 && r_hi <= 0)
    error('exactLowerObjFast failed to bracket the unique inner dual root.');
end

for root_iter = 1:80
    lambda = 0.5 * (lambda_lo + lambda_hi);
    [~, pes, ~, ~] = eu_response(lambda);
    residual = (pi_0 - lambda) / a - sum(pes);
    if residual >= 0
        lambda_lo = lambda;
    else
        lambda_hi = lambda;
    end
end

lambda = 0.5 * (lambda_lo + lambda_hi);
[pg, pes, pb, ps] = eu_response(lambda);
y = sum(pes);
obj_star = sum(0.5 .* c(:) .* pg.^2 + b(:) .* pg + ...
    pi_max .* pb - pi_min .* ps + 0.5 .* a .* pes.^2) + ...
    0.5 * a * y^2 - pi_0 * y;

meta = struct('problem',0, 'info','analytic_bisection', 'solver','analytic_bisection', ...
    'lambda',lambda, 'residual', (pi_0-lambda)/a-sum(pes), ...
    'root_iterations',root_iter, 'bracket_expansions',expand_iter, ...
    'objective_source','exactLowerObjFast');

    function [pg_v, pes_v, pb_v, ps_v] = eu_response(lambda_v)
        pg_v = zeros(N,1);
        pes_v = zeros(N,1);
        pb_v = zeros(N,1);
        ps_v = zeros(N,1);
        tol_F = 1e-10;
        for j = 1:N
            pg_buyer = min(max((pi_max-b(j))/c(j), pg_min(j)), pg_max(j));
            pes_buyer = (lambda_v-pi_max)/a;
            F_buyer = D(j) + pes_buyer - pg_buyer;
            pg_seller = min(max((pi_min-b(j))/c(j), pg_min(j)), pg_max(j));
            pes_seller = (lambda_v-pi_min)/a;
            F_seller = D(j) + pes_seller - pg_seller;
            if F_buyer > tol_F
                pg_v(j) = pg_buyer; pes_v(j) = pes_buyer; pb_v(j) = F_buyer;
            elseif F_seller < -tol_F
                pg_v(j) = pg_seller; pes_v(j) = pes_seller; ps_v(j) = -F_seller;
            else
                pes_v(j) = min(max((lambda_v-b(j)-c(j)*D(j))/(c(j)+a), ...
                    pg_min(j)-D(j)), pg_max(j)-D(j));
                pg_v(j) = D(j) + pes_v(j);
            end
        end
    end
end
