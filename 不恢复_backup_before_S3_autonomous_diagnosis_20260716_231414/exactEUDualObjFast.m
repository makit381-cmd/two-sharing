function obj_star = exactEUDualObjFast(c,b,a,D,pg_max,pg_min,pi_max,pi_min,lambda)
% 精确 EU 对偶子问题目标。每个 lambda 均按 exactLowerObj 的同一约束和目标求解，
% 但使用当前分段二次模型的闭式响应，避免记录轨迹时重复调用 Gurobi。

lambda = lambda(:);
N = numel(c);
obj_star = zeros(size(lambda));
tol_F = 1e-10;

for q = 1:numel(lambda)
    lambda_q = lambda(q);
    pg = zeros(N,1);
    pes = zeros(N,1);
    pb = zeros(N,1);
    ps = zeros(N,1);

    for j = 1:N
        pg_buyer = min(max((pi_max-b(j))/c(j), pg_min(j)), pg_max(j));
        pes_buyer = (lambda_q-pi_max)/a;
        F_buyer = D(j) + pes_buyer - pg_buyer;

        pg_seller = min(max((pi_min-b(j))/c(j), pg_min(j)), pg_max(j));
        pes_seller = (lambda_q-pi_min)/a;
        F_seller = D(j) + pes_seller - pg_seller;

        if F_buyer > tol_F
            pg(j) = pg_buyer;
            pes(j) = pes_buyer;
            pb(j) = F_buyer;
        elseif F_seller < -tol_F
            pg(j) = pg_seller;
            pes(j) = pes_seller;
            ps(j) = -F_seller;
        else
            pes(j) = min(max((lambda_q-b(j)-c(j)*D(j))/(c(j)+a), ...
                pg_min(j)-D(j)), pg_max(j)-D(j));
            pg(j) = D(j) + pes(j);
        end
    end

    obj_star(q) = sum(0.5 .* c(:) .* pg.^2 + b(:) .* pg + ...
        pi_max .* pb - pi_min .* ps + 0.5 .* a .* pes.^2 - lambda_q .* pes);
end
end
