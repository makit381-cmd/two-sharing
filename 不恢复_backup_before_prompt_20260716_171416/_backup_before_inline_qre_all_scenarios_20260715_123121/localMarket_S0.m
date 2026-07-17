function res = localMarket_S0(c,b,a,D,init_rho,pg_max,pg_min,init_pi, ...
    init_pes,init_pg,init_pb,init_ps,init_pesc,errTol, ...
    pi_0,pi_max,pi_min,record_trace,beta_qre,qre_error_rel,qre_z_cap,inner_momentum,max_inner_iter)

% S0: 同步等待处理内部延迟 + 即时输出
h0 = 6;
max_inner_wall_iter = (h0 + 1) * max_inner_iter;
history_limit = max_inner_wall_iter + 2;
A = min(history_limit, 128);
min_inner_iter = 100;
stable_inner_window = 5;

% 记录真实同步价格更新
sync_lambda = nan(A,1);
sync_lambda(1) = init_pi;
sync_lambda_step = nan(A,1);
sync_cnt = 0;

cnt = 0;
N = length(b);

theta = 1.8;
rho_max = theta / (N + 1);
rho_min = 0.05 * rho_max;
rho = min(max(init_rho, rho_min), rho_max);

delay = randi([0,h0],N,1);
count = zeros(N,1);
mask = false(N,1);

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

res.success = 1;

while true

    if h + 1 > A
        nextA = growHistoryCapacity(A, h + 1, history_limit);
        sync_lambda(nextA,1) = NaN;
        sync_lambda_step(nextA,1) = NaN;
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
        if ~mask(j)
            if count(j) == delay(j)
                [fast_pg(j), fast_pes(j), fast_pb(j), fast_ps(j)] = ...
                    solveBranch(j, val_lambda(h));

                mask(j) = true;
                count(j) = 0;
                delay(j) = randi([0,h0]);
            else
                count(j) = count(j) + 1;
            end
        end
    end

    val_pes(:,h+1) = fast_pes;
    val_pg(:,h+1)  = fast_pg;
    val_pb(:,h+1)  = fast_pb;
    val_ps(:,h+1)  = fast_ps;

    val_pesc(h+1) = (pi_0 - val_lambda(h)) / a;

    val_lambda(h+1) = val_lambda(h);
    val_lambda_au(h+1) = val_lambda_au(h);

    if all(mask)

        val_lambda_au(h+1) = val_lambda(h) + ...
            alpha * (val_pesc(h+1) - sum(fast_pes));
        val_lambda(h+1) = val_lambda_au(h+1) + inner_momentum * ...
            (val_lambda_au(h+1) - val_lambda_au(h));

        sync_cnt = sync_cnt + 1;

        % 保存本次真实同步更新后的价格
        sync_lambda(sync_cnt+1) = val_lambda(h+1);

        % 相邻两次真实同步价格更新量
        sync_lambda_step(sync_cnt) = abs( ...
            sync_lambda(sync_cnt+1) - sync_lambda(sync_cnt));

        mask(:) = false;
    end

    if record_trace
        val_violation(h+1) = abs( ...
            (val_pesc(h+1) - sum(fast_pes)) / 1e3);

        val_obj(h+1) = originalLowerObj( ...
            fast_pg, fast_pb, fast_ps, fast_pes, sum(fast_pes));
    end

    if sync_cnt >= min_inner_iter && sync_cnt >= stable_inner_window
        recent_step = sync_lambda_step( ...
            sync_cnt-stable_inner_window+1:sync_cnt);

        if all(recent_step <= errTol) && ...
                abs((val_pesc(h+1) - sum(fast_pes)) / 1e3) <= 1e-2
            break;
        end
    end

    if sync_cnt >= max_inner_iter
        res.success = 0;
        break;
    end

    if sync_cnt > 1
        if (val_lambda(h+1) - val_lambda(h)) * ...
                (val_lambda(h) - val_lambda(max(h-1,1))) < 0

            rho = max(rho / 2, rho_min);
            cnt = 0;

        elseif sync_cnt > 10

            recent_lambda = val_lambda(max(1,h-5):h+1);

            if all(diff(recent_lambda) >= 0) || ...
                    all(diff(recent_lambda) <= 0)
                rho = min(rho * 1.2, rho_max);
            end

        elseif sync_cnt > 100

            if (val_lambda(h+1) - val_lambda(h) > errTol && ...
                    val_lambda(max(h-1,1)) - val_lambda(h) > errTol) || ...
               (val_lambda(h) - val_lambda(h+1) > errTol && ...
                    val_lambda(h) - val_lambda(max(h-1,1)) > errTol)

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

res.pes = fast_pes;
res.pg  = fast_pg;
res.pb  = fast_pb;
res.ps  = fast_ps;
res.pesc = sum(fast_pes);
res.xj2 = sum(res.pes .^ 2);
res.iter = h;
res.price_updates = sync_cnt;

    function [pg, pes, pb, ps] = solveBranch(j, lambda)
        [pg, pes, pb, ps] = boundedQREBranch(c(j), b(j), a, D(j), ...
            pg_min(j), pg_max(j), lambda, pi_max, pi_min, beta_qre(j), ...
            qre_error_rel, qre_z_cap);
    end

    function obj = originalLowerObj(pg_v, pb_v, ps_v, pes_v, y_v)

        obj = sum( ...
            0.5 .* c(:) .* pg_v(:).^2 ...
            + b(:) .* pg_v(:) ...
            + pi_max .* pb_v(:) ...
            - pi_min .* ps_v(:) ...
            + 0.5 .* a .* pes_v(:).^2) ...
            + 0.5 * a * y_v^2 ...
            - pi_0 * y_v;

    end

end
