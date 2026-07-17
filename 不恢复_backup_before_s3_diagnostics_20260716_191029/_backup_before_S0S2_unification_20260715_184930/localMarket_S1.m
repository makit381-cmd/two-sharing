function res = localMarket_S1(c,b,a,D,init_rho,pg_max,pg_min,init_pi, ...
    init_pes,init_pg,init_pb,init_ps,init_pesc,errTol, ...
    pi_0,pi_max,pi_min,record_trace,beta_qre,~,qre_z_cap,inner_momentum,max_inner_iter)

% S1: 异步保持 + 即时输出
history_limit = max_inner_iter + 2;
A = min(history_limit, 128);
min_inner_iter = 100;
stable_inner_window = 5;
lambda_step = nan(A,1);
cnt = 0;
N = length(b);

theta = 1.8;
rho_max = theta / (N + 1);
rho_min = 0.05 * rho_max;
rho = min(max(init_rho, rho_min), rho_max);

%% 异步保持设置
h0 = 6;
delay = randi([0,h0],N,1);
count = zeros(N,1);

h = 1;

val_lambda = zeros(A,1);
val_lambda(1) = init_pi;
val_lambda_au = zeros(A,1);
val_lambda_au(1) = init_pi;

val_pes = zeros(N,A);
val_pes(:,1) = init_pes;

val_pg = zeros(N,A);
val_pg(:,1) = init_pg;

val_pb = zeros(N,A);
val_pb(:,1) = init_pb;

val_ps = zeros(N,A);
val_ps(:,1) = init_ps;

val_pesc = zeros(A,1);
val_pesc(1) = init_pesc;

val_violation = zeros(A,1);
val_violation(1) = abs((val_pesc(1) - sum(init_pes)) / 1e3);

val_obj = nan(A,1);
val_obj(1) = originalLowerObj(val_pg(:,1), val_pb(:,1), ...
    val_ps(:,1), val_pes(:,1), val_pesc(1));

fast_pg  = init_pg;
fast_pes = init_pes;
fast_pb  = init_pb;
fast_ps  = init_ps;
tol_F = 1e-8;

res.success = 1;

while true

    if h + 1 > A
        nextA = growHistoryCapacity(A, h + 1, history_limit);
        lambda_step(nextA,1) = NaN;
        val_lambda(nextA,1) = 0;
        val_lambda_au(nextA,1) = 0;
        val_pes(:,nextA) = 0;
        val_pg(:,nextA) = 0;
        val_pb(:,nextA) = 0;
        val_ps(:,nextA) = 0;
        val_pesc(nextA,1) = 0;
        val_violation(nextA,1) = 0;
        val_obj(nextA,1) = NaN;
        A = nextA;
    end

    alpha = rho * a;

    for j = 1:N
        if count(j) < delay(j)
            count(j) = count(j) + 1;
            fast_pes(j) = val_pes(j,h);
            fast_pg(j)  = val_pg(j,h);
            fast_pb(j)  = val_pb(j,h);
            fast_ps(j)  = val_ps(j,h);
        else
            % 与 S3 对齐：内联未认证 QRE 解析响应，保留随机数顺序与分支。
            z_pg = max(min(randn(1), qre_z_cap), -qre_z_cap);
            z_pes = max(min(randn(1), qre_z_cap), -qre_z_cap);
            c_j = c(j); b_j = b(j); D_j = D(j);
            pg_min_j = pg_min(j); pg_max_j = pg_max(j); beta_j = beta_qre(j);
            z_middle = (sqrt(c_j) * z_pg + sqrt(a) * z_pes) / sqrt(c_j + a);
            pg_buyer = min(max((pi_max - b_j) / c_j + sqrt(1 / (beta_j * c_j)) * z_pg, pg_min_j), pg_max_j);
            pes_buyer = (val_lambda(h) - pi_max) / a + sqrt(1 / (beta_j * a)) * z_pes;
            F_buyer = D_j + pes_buyer - pg_buyer;
            pg_seller = min(max((pi_min - b_j) / c_j + sqrt(1 / (beta_j * c_j)) * z_pg, pg_min_j), pg_max_j);
            pes_seller = (val_lambda(h) - pi_min) / a + sqrt(1 / (beta_j * a)) * z_pes;
            F_seller = D_j + pes_seller - pg_seller;
            if F_buyer > tol_F
                fast_pg(j) = pg_buyer; fast_pes(j) = pes_buyer; fast_pb(j) = F_buyer; fast_ps(j) = 0;
            elseif F_seller < -tol_F
                fast_pg(j) = pg_seller; fast_pes(j) = pes_seller; fast_pb(j) = 0; fast_ps(j) = -F_seller;
            else
                pes_unclipped = (val_lambda(h) - b_j - c_j * D_j) / (c_j + a) + sqrt(1 / (beta_j * (c_j + a))) * z_middle;
                fast_pes(j) = min(max(pes_unclipped, pg_min_j - D_j), pg_max_j - D_j);
                fast_pg(j) = D_j + fast_pes(j); fast_pb(j) = 0; fast_ps(j) = 0;
            end
            count(j) = 0;
            delay(j) = randi([0,h0]);
        end
    end

    val_pes(:,h+1) = fast_pes;
    val_pg(:,h+1)  = fast_pg;
    val_pb(:,h+1)  = fast_pb;
    val_ps(:,h+1)  = fast_ps;

    val_pesc(h+1) = (pi_0 - val_lambda(h)) / a;

    val_lambda_au(h+1) = val_lambda(h) + ...
        alpha * (val_pesc(h+1) - sum(val_pes(:,h+1)));
    val_lambda(h+1) = val_lambda_au(h+1) + inner_momentum * ...
        (val_lambda_au(h+1) - val_lambda_au(h));

    if record_trace
        val_violation(h+1) = abs((val_pesc(h+1) - sum(fast_pes)) / 1e3);
        val_obj(h+1) = originalLowerObj(fast_pg, fast_pb, ...
            fast_ps, fast_pes, sum(fast_pes));
    end

    lambda_step(h) = abs(val_lambda(h+1) - val_lambda(h));

    if h >= min_inner_iter
        if all(lambda_step(h-stable_inner_window+1:h) <= errTol) && ...
                abs((val_pesc(h+1) - sum(fast_pes)) / 1e3) <= 0.01
            break;
        end
    end

    if h >= max_inner_iter
        res.success = 0;
        break;
    end

    if h > 1
        if (val_lambda(h+1) - val_lambda(h)) * ...
                (val_lambda(h) - val_lambda(h-1)) < 0
            rho = min(rho / 2,rho_min);
            cnt = 0;
        elseif h > 10 && ...
                (all(diff(val_lambda(h-5:h+1)) > 0) || ...
                 all(diff(val_lambda(h-5:h+1)) < 0))
            rho = min(rho * 1.2,rho_max);
        elseif h > 100
            if (val_lambda(h+1) - val_lambda(h) > errTol && ...
                    val_lambda(h-1) - val_lambda(h) > errTol) || ...
               (val_lambda(h) - val_lambda(h+1) > errTol && ...
                    val_lambda(h) - val_lambda(h-1) > errTol)
                cnt = cnt + 1;
                if cnt >= 5
                    rho = max(rho / 2,rho_min);
                    cnt = 0;
                end
            end
        end
    end

    h = h + 1;

end

if record_trace
    trace_end = h + 1;
    res.trace_h = (0:h)';
    res.trace_violation = abs(val_violation(1:trace_end));
    res.trace_obj = val_obj(1:trace_end);
    res.trace_lambda = val_lambda(1:trace_end);
else
    res.trace_h = [];
    res.trace_violation = [];
    res.trace_obj = [];
    res.trace_lambda = [];
end

res.pi = val_lambda(h+1);

res.pes = fast_pes;
res.pg  = fast_pg;
res.pb  = fast_pb;
res.ps  = fast_ps;

res.pesc = sum(fast_pes);
res.xj2 = sum(res.pes .^ 2);
res.iter = h;
res.price_updates = h;

    function obj = originalLowerObj(pg_v, pb_v, ps_v, pes_v, y_v)

        obj = sum(0.5 .* c(:) .* pg_v(:).^2 ...
            + b(:) .* pg_v(:) ...
            + pi_max .* pb_v(:) ...
            - pi_min .* ps_v(:) ...
            + 0.5 .* a .* pes_v(:).^2) ...
            + 0.5 * a * y_v^2 ...
            - pi_0 * y_v;

    end

end
