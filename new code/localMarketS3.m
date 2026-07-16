function res = localMarketS3(c,b,a,D,init_rho,pg_max,pg_min,init_pi, ...
    init_pes,init_pg,init_pb,init_ps,init_pesc,errTol, ...
    pi_0,pi_max,pi_min,record_trace,beta_qre)

max_inner_iter = 1500;
A = max_inner_iter + 2;
min_inner_iter = 220;
stable_inner_window = 5;
errTol_output_CV_MW = 5e-3;
lambda_step = nan(A,1);
cnt = 0;
N = length(b);

theta = 1.8;
rho_max = theta / (N + 1);
rho_min = 0.05 * rho_max;
rho = min(max(init_rho, rho_min), rho_max);
%% 加异步
h0 = 6;
delay = randi([0,h0],N,1);
count = zeros(N,1);

h = 1;

val_lambda = zeros(A,1);
val_lambda(1) = init_pi;

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
val_violation(1) = abs((val_pesc(1) - sum(init_pes)) / 1e3);%约束违反单位用MW

val_obj = nan(A,1);
val_obj(1) = originalLowerObj(val_pg(:,1), val_pb(:,1), ...
    val_ps(:,1), val_pes(:,1), val_pesc(1));

fast_pg  = init_pg;
fast_pes = init_pes;
fast_pb  = init_pb;
fast_ps  = init_ps;

sum_pes = zeros(N,1);
sum_pg  = zeros(N,1);
sum_pb  = zeros(N,1);
sum_ps  = zeros(N,1);
avg_cnt = 0;

res.success = 1;

while true

    alpha = rho * a;%统一单位

    for j = 1:N
        if count(j) < delay(j)
            count(j) = count(j) + 1;
            fast_pes(j) = val_pes(j,h);
            fast_pg(j)  = val_pg(j,h);
            fast_pb(j)  = val_pb(j,h);
            fast_ps(j)  = val_ps(j,h);
        else
            [fast_pg(j), fast_pes(j), fast_pb(j), fast_ps(j)] = ...
                solveBranch(j, val_lambda(h));
            count(j) = 0;
            delay(j) = randi([0,h0]);
        end

    end

    val_pes(:,h+1) = fast_pes;
    val_pg(:,h+1)  = fast_pg;
    val_pb(:,h+1)  = fast_pb;
    val_ps(:,h+1)  = fast_ps;

    val_pesc(h+1) = (pi_0 - val_lambda(h)) / a;

    avg_cnt = avg_cnt + 1;

    sum_pes = sum_pes + fast_pes;
    sum_pg  = sum_pg  + fast_pg;
    sum_pb  = sum_pb  + fast_pb;
    sum_ps  = sum_ps  + fast_ps;

    bar_pes = sum_pes / avg_cnt;
    bar_pg  = sum_pg  / avg_cnt;
    bar_pb  = sum_pb  / avg_cnt;
    bar_ps  = sum_ps  / avg_cnt;
    bar_pesc = sum(bar_pes);

    val_lambda(h+1) = val_lambda(h) + ...
        alpha * (val_pesc(h+1) - sum(val_pes(:,h+1)));

    lower_output_CV_MW = abs((val_pesc(h+1) - sum(fast_pes)) / 1e3);

    if record_trace
        val_violation(h+1) = lower_output_CV_MW;%仍然是约束违反单位取MW
        val_obj(h+1) = originalLowerObj(bar_pg, bar_pb, ...
            bar_ps, bar_pes, sum(bar_pes));
    end

    lambda_step(h) = abs(val_lambda(h+1) - val_lambda(h));

    if h >= min_inner_iter%内层最小迭代次数
        if all(lambda_step(h-stable_inner_window+1:h) <= errTol) && ...
                lower_output_CV_MW <= errTol_output_CV_MW
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
            rho = max(rho / 2, rho_min);
            cnt = 0;
        elseif h > 10 && ...
                (all(diff(val_lambda(h-5:h+1)) > 0) || ...
                 all(diff(val_lambda(h-5:h+1)) < 0))
            rho = min(rho * 1.2, rho_max);
        elseif h > 100
            if (val_lambda(h+1) - val_lambda(h) > errTol && ...
                    val_lambda(h-1) - val_lambda(h) > errTol) || ...
               (val_lambda(h) - val_lambda(h+1) > errTol && ...
                    val_lambda(h) - val_lambda(h-1) > errTol)
                cnt = cnt + 1;
                if cnt >= 5
                    rho = max(rho / 2, rho_min);
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

res.pes = bar_pes;
res.pg  = bar_pg;
res.pb  = bar_pb;
res.ps  = bar_ps;

res.pesc = bar_pesc;
res.xj2 = sum(res.pes .^ 2);
res.iter = h;

    function [pg, pes, pb, ps] = solveBranch(j, lambda)

        tol_F = 1e-8;
        z_pg = randn(1);
        z_pes = randn(1);
        z_middle = (sqrt(c(j)) * z_pg + sqrt(a) * z_pes) / sqrt(c(j) + a);
        pg_buyer = max(min((pi_max - b(j)) / c(j) + sqrt(1 / (beta_qre(j) * c(j))) * z_pg, pg_max(j)), pg_min(j)) ...
            ;
        pes_buyer = (lambda - pi_max) / a ...
            + sqrt(1 / (beta_qre(j) * a)) * z_pes;

        F_buyer = D(j) + pes_buyer - pg_buyer;

        pg_seller = max(min((pi_min - b(j)) / c(j) + sqrt(1 / (beta_qre(j) * c(j))) * z_pg, pg_max(j)), pg_min(j)) ...
           ;
        pes_seller = (lambda - pi_min) / a ...
            + sqrt(1 / (beta_qre(j) * a)) * z_pes;

        F_seller = D(j) + pes_seller - pg_seller;

        if F_buyer > tol_F

            pg = pg_buyer;
            pes = pes_buyer;

            pb = D(j) + pes - pg;
            ps = 0;

        elseif F_seller < -tol_F

            pg = pg_seller;
            pes = pes_seller;

            pb = 0;
            ps = pg - D(j) - pes;

        else

            pes = (lambda - b(j) - c(j) * D(j)) / (c(j) + a) ...
                + sqrt(1 / (beta_qre(j) * (c(j) + a))) * z_middle;

            pes = max(min(pes, pg_max(j) - D(j)), pg_min(j) - D(j));

            pg = D(j) + pes;
            pb = 0;
            ps = 0;

        end

    end

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
