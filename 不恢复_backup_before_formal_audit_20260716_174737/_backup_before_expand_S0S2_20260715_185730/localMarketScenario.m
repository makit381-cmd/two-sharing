function res = localMarketScenario(c,b,a,D,init_rho,pg_max,pg_min,init_pi, ...
    init_pes,init_pg,init_pb,init_ps,init_pesc,errTol,inner_cv_tol_kW, ...
    pi_0,pi_max,pi_min,record_trace,beta_qre,~,~,qre_z_cap,~, ...
    inner_momentum,max_inner_iter,h0,min_inner_iter,inner_avg_start_iter,inner_cv_stop_enabled, ...
    qre_certificate_enabled,agg_certificate_enabled,progress_print_enabled, ...
    progress_print_inner_every,outer_iter,community_id,inner_async_hold,use_ergodic_average)
% 三个非 S3 情形的统一内层实现。
% 价格更新始终使用瞬时 raw Agg--EU 残差；上传、RE/CV 诊断按情形选择瞬时或遍历平均。

if qre_certificate_enabled || agg_certificate_enabled
    error('当前 profile 为 unchecked QRE，不能打开 epsilon 证书开关。');
end
N = numel(b);
alpha = init_rho * a;
lambda = init_pi;
lambda_aux_prev = init_pi;
fast_pg = init_pg; fast_pes = init_pes; fast_pb = init_pb; fast_ps = init_ps;

delay = randi([0,h0], N, 1);
count = zeros(N,1);
ready = false(N,1);

sum_pes = zeros(N,1); sum_pg = zeros(N,1); sum_pb = zeros(N,1); sum_ps = zeros(N,1);
sum_pesc = 0; avg_cnt = 0;
bar_pes = init_pes; bar_pg = init_pg; bar_pb = init_pb; bar_ps = init_ps; bar_pesc = init_pesc;

if record_trace
    trace_h = 0; trace_violation = abs(init_pesc - sum(init_pes));
    trace_obj = communityPrimalObj(init_pg,init_pb,init_ps,init_pes,init_pesc);
    trace_lambda = init_pi; trace_lambda_used = init_pi;
end

price_updates = 0;
wall_iter = 0;
res.success = 1;
res.stop_reason = 'running';
if progress_print_enabled
    fprintf('[S%d inner start] k=%d, community=%d, EU=%d\n', ...
        1 + use_ergodic_average + (~inner_async_hold), outer_iter, community_id, N);
end

while true
    wall_iter = wall_iter + 1;
    lambda_used = lambda;
    if inner_async_hold
        for j = 1:N
            if count(j) < delay(j)
                count(j) = count(j) + 1;
            else
                updateEU(j, lambda_used);
                count(j) = 0;
                delay(j) = randi([0,h0]);
            end
        end
        complete_round = true;
    else
        % 同步等待：一轮内所有 EU 都以相同 lambda 求解；全体完成后才更新 lambda。
        for j = 1:N
            if ~ready(j)
                if count(j) < delay(j)
                    count(j) = count(j) + 1;
                else
                    updateEU(j, lambda_used);
                    ready(j) = true;
                end
            end
        end
        complete_round = all(ready);
    end
    if ~complete_round
        continue;
    end

    price_updates = price_updates + 1;
    pesc_now = (pi_0 - lambda_used) / a; % Algorithm 2 的原始 Agg 输出，绝不重写为 sum(pes)。
    if use_ergodic_average && price_updates >= inner_avg_start_iter
        avg_cnt = avg_cnt + 1;
        sum_pes = sum_pes + fast_pes; sum_pg = sum_pg + fast_pg;
        sum_pb = sum_pb + fast_pb; sum_ps = sum_ps + fast_ps; sum_pesc = sum_pesc + pesc_now;
        bar_pes = sum_pes / avg_cnt; bar_pg = sum_pg / avg_cnt;
        bar_pb = sum_pb / avg_cnt; bar_ps = sum_ps / avg_cnt; bar_pesc = sum_pesc / avg_cnt;
    else
        bar_pes = fast_pes; bar_pg = fast_pg; bar_pb = fast_pb; bar_ps = fast_ps; bar_pesc = pesc_now;
    end

    lambda_aux = lambda_used + alpha * (pesc_now - sum(fast_pes));
    lambda_next = lambda_aux + inner_momentum * (lambda_aux - lambda_aux_prev);
    lambda_step = abs(lambda_next - lambda_used);
    cv_now_kW = abs(bar_pesc - sum(bar_pes));

    if record_trace
        trace_h(end+1,1) = price_updates;
        trace_violation(end+1,1) = cv_now_kW;
        trace_obj(end+1,1) = communityPrimalObj(bar_pg,bar_pb,bar_ps,bar_pes,bar_pesc);
        trace_lambda(end+1,1) = lambda_next;
        trace_lambda_used(end+1,1) = lambda_used;
    end
    if progress_print_enabled && (price_updates == 1 || mod(price_updates,progress_print_inner_every) == 0)
        fprintf('[inner] k=%d, community=%d, h=%d, CV=%.3f kW, lambdaStep=%.3e\n', ...
            outer_iter,community_id,price_updates,cv_now_kW,lambda_step);
    end

    lambda_aux_prev = lambda_aux;
    lambda = lambda_next;
    if ~inner_async_hold
        ready(:) = false; count(:) = 0;
        delay = randi([0,h0], N, 1);
    end

    cv_pass = ~logical(inner_cv_stop_enabled) || (cv_now_kW <= inner_cv_tol_kW);
    stop_start = min_inner_iter;
    if use_ergodic_average
        stop_start = max(stop_start, inner_avg_start_iter);
    end
    if price_updates >= stop_start && lambda_step <= errTol && cv_pass
        if inner_cv_stop_enabled
            res.stop_reason = 'price_step_and_selected_cv';
        else
            res.stop_reason = 'price_step_only';
        end
        break;
    end
    if isfinite(max_inner_iter) && price_updates >= max_inner_iter
        res.success = 0;
        res.stop_reason = 'max_inner_iter';
        break;
    end
end

res.pi = lambda; res.pes = bar_pes; res.pg = bar_pg; res.pb = bar_pb; res.ps = bar_ps;
res.pesc = bar_pesc; res.pesc_agg_raw = bar_pesc; res.pes_sum = sum(bar_pes);
res.inner_cv_kW = abs(res.pesc - res.pes_sum);
res.inner_avg_cnt = avg_cnt; res.inner_avg_start_iter = inner_avg_start_iter;
res.inner_cv_stop_enabled = logical(inner_cv_stop_enabled);
res.iter = wall_iter; res.price_updates = price_updates;
if record_trace
    res.trace_h=trace_h; res.trace_violation=trace_violation; res.trace_obj=trace_obj;
    res.trace_lambda=trace_lambda; res.trace_lambda_used=trace_lambda_used;
else
    res.trace_h=[]; res.trace_violation=[]; res.trace_obj=[]; res.trace_lambda=[]; res.trace_lambda_used=[];
end

    function updateEU(j,lambda_used_now)
        z_pg = max(min(randn(1),qre_z_cap),-qre_z_cap);
        z_pes = max(min(randn(1),qre_z_cap),-qre_z_cap);
        cj=c(j); bj=b(j); Dj=D(j); amin=a; betaj=beta_qre(j);
        z_mid=(sqrt(cj)*z_pg+sqrt(amin)*z_pes)/sqrt(cj+amin);
        pg_b=min(max((pi_max-bj)/cj+sqrt(1/(betaj*cj))*z_pg,pg_min(j)),pg_max(j));
        pes_b=(lambda_used_now-pi_max)/amin+sqrt(1/(betaj*amin))*z_pes;
        F_b=Dj+pes_b-pg_b;
        pg_s=min(max((pi_min-bj)/cj+sqrt(1/(betaj*cj))*z_pg,pg_min(j)),pg_max(j));
        pes_s=(lambda_used_now-pi_min)/amin+sqrt(1/(betaj*amin))*z_pes;
        F_s=Dj+pes_s-pg_s;
        if F_b > 1e-8
            fast_pg(j)=pg_b; fast_pes(j)=pes_b; fast_pb(j)=F_b; fast_ps(j)=0;
        elseif F_s < -1e-8
            fast_pg(j)=pg_s; fast_pes(j)=pes_s; fast_pb(j)=0; fast_ps(j)=-F_s;
        else
            pes0=(lambda_used_now-bj-cj*Dj)/(cj+amin)+sqrt(1/(betaj*(cj+amin)))*z_mid;
            fast_pes(j)=min(max(pes0,pg_min(j)-Dj),pg_max(j)-Dj);
            fast_pg(j)=Dj+fast_pes(j); fast_pb(j)=0; fast_ps(j)=0;
        end
    end
    function obj = communityPrimalObj(pg_v,pb_v,ps_v,pes_v,pesc_v)
        obj=sum(0.5*c(:).*pg_v(:).^2+b(:).*pg_v(:)+pi_max*pb_v(:)-pi_min*ps_v(:)+0.5*a*pes_v(:).^2) ...
            +0.5*a*pesc_v.^2-pi_0*pesc_v;
    end
end
