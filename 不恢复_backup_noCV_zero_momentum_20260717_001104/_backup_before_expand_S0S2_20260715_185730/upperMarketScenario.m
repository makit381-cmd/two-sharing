function upperMarketScenario(scenario_id)
% S0--S2 公共外层实现：统一使用 S3 的 raw Agg、RE、总 CV 和存档口径。
% S0=同步等待/不平均，S1=异步保持/不平均，S2=同步等待/遍历平均。

if ~ismember(scenario_id,[0 1 2])
    error('scenario_id 必须为 0、1 或 2。');
end
P=load('params.mat');
params_snapshot=P;
switch scenario_id
    case 0
        scenario_name='S0'; outer_async_hold=false; use_ergodic_average=false;
    case 1
        scenario_name='S1'; outer_async_hold=true;  use_ergodic_average=false;
    case 2
        scenario_name='S2'; outer_async_hold=false; use_ergodic_average=true;
end
if P.qre_certificate_enabled || P.agg_certificate_enabled
    error('%s 目前仅迁移 unchecked-QRE profile；请先关闭证书开关。',scenario_name);
end

rng(P.rng_seeds.outer_schedule,'twister');
finite_mask=isfinite(P.F_l); finite_idx=find(finite_mask); nfinite=numel(finite_idx);
pi=P.init_pi; pi0=P.init_pi_0; pes=P.init_pes; pg=P.init_pg; pb=P.init_pb; ps=P.init_ps; pesc=P.init_pesc;
piPB=P.init_pi_PB; piLpos=P.init_pi_l_pos; piLneg=P.init_pi_l_neg;
auxPB=piPB; auxLpos=piLpos; auxLneg=piLneg;
count=zeros(P.num_LESMs,1); delay=zeros(P.num_LESMs,1);
if outer_async_hold
    for i=1:P.num_LESMs, delay(i)=randi([0,P.k0_max(i)]); end
end

cap=128;
re_hist=nan(cap,1); cv_hist=nan(cap,1); pb_hist=nan(cap,1); maxline_hist=nan(cap,1);
comp_hist=nan(cap,1); price_hist=nan(cap,1); line_inst=nan(nfinite,cap); line_metric=nan(nfinite,cap);
sum_pes=zeros(P.num_prosumers,1); sum_pg=sum_pes; sum_pb=sum_pes; sum_ps=sum_pes; sum_pesc=zeros(P.num_LESMs,1); avg_cnt=0;
inner_k=struct('h',[],'val_violation',[],'val_obj',[],'obj_ref',[],'RE_obj',[]);
over_trace=false; num_local=0; total_inner_wall=0; total_inner_updates=0; max_wall=0; max_updates=0; total_time=0;
ref=exactUpperObj(P.a,P.c,P.b,P.pi_max,P.pi_min,P.D,P.F_l,P.PTDF,P.start_idx,P.end_idx,P.pg_max,P.pg_min,P.a_extend);
refscale=max(abs(ref),eps);
k=0; stop_reason='running';
while true
    k=k+1;
    if k>cap
        newcap=ceil(1.5*cap);
        re_hist(newcap,1)=NaN; cv_hist(newcap,1)=NaN; pb_hist(newcap,1)=NaN; maxline_hist(newcap,1)=NaN;
        comp_hist(newcap,1)=NaN; price_hist(newcap,1)=NaN; line_inst(:,newcap)=NaN; line_metric(:,newcap)=NaN; cap=newcap;
    end
    if P.progress_print_enabled && (k==1 || mod(k,P.progress_print_outer_every)==0)
        fprintf('[%s outer start] k=%d, local calls=%d\n',scenario_name,k,num_local);
    end
    for i=1:P.num_LESMs
        if outer_async_hold && count(i)<delay(i)
            count(i)=count(i)+1;
            continue;
        end
        idx=P.start_idx(i):P.end_idx(i);
        record_this=ismember(k,P.record_k) && i==P.record_inner_comm && ~over_trace;
        local_tic=tic;
        r=localMarketScenario(P.c(idx),P.b(idx),P.a(i),P.D(idx),P.init_rho(i),P.pg_max(idx),P.pg_min(idx),pi(i), ...
            pes(idx),pg(idx),pb(idx),ps(idx),pesc(i),P.errTol_LESMs(i),P.inner_cv_tol_kW(i),pi0(i),P.pi_max,P.pi_min,record_this, ...
            P.beta_qre(idx),P.qre_epsilon,P.agg_epsilon,P.qre_z_cap,P.qre_max_backoffs,P.inner_momentum,P.max_inner_iter,P.h0_max, ...
            P.min_inner_iter,P.inner_avg_start_iter,P.inner_cv_stop_enabled,P.qre_certificate_enabled,P.agg_certificate_enabled, ...
            P.progress_print_enabled,P.progress_print_inner_every,k,i,outer_async_hold,use_ergodic_average);
        total_time=total_time+toc(local_tic); num_local=num_local+1; total_inner_wall=total_inner_wall+r.iter; total_inner_updates=total_inner_updates+r.price_updates;
        max_wall=max(max_wall,r.iter); max_updates=max(max_updates,r.price_updates);
        pi(i)=r.pi; pes(idx)=r.pes; pg(idx)=r.pg; pb(idx)=r.pb; ps(idx)=r.ps; pesc(i)=r.pesc;
        if outer_async_hold
            count(i)=0; delay(i)=randi([0,P.k0_max(i)]);
        end
        if record_this
            [obj_ref,meta]=exactLowerObjFast(P.c(idx),P.b(idx),P.a(i),P.D(idx),P.pg_max(idx),P.pg_min(idx),pi0(i),P.pi_max,P.pi_min);
            inner_k.h=r.trace_h(:); inner_k.val_violation=r.trace_violation(:); inner_k.val_obj=r.trace_obj(:);
            inner_k.obj_ref=repmat(obj_ref,size(inner_k.val_obj)); inner_k.RE_obj=100*abs(inner_k.val_obj-obj_ref)/max(abs(obj_ref),eps);
            inner_k.CV_inner_selected_kW=inner_k.val_violation; inner_k.obj_ref_meta=meta; inner_k.record_outer_iter=k; inner_k.record_comm=i;
            over_trace=true;
        end
    end

    if use_ergodic_average && k>=P.outer_avg_start_iter
        avg_cnt=avg_cnt+1; sum_pes=sum_pes+pes; sum_pg=sum_pg+pg; sum_pb=sum_pb+pb; sum_ps=sum_ps+ps; sum_pesc=sum_pesc+pesc;
        m_pes=sum_pes/avg_cnt; m_pg=sum_pg/avg_cnt; m_pb=sum_pb/avg_cnt; m_ps=sum_ps/avg_cnt; m_pesc=sum_pesc/avg_cnt;
    else
        m_pes=pes; m_pg=pg; m_pb=pb; m_ps=ps; m_pesc=pesc;
    end
    cost=sum(0.5*P.c.*m_pg.^2+P.b.*m_pg+P.pi_max*m_pb-P.pi_min*m_ps+0.5*P.a_extend.*m_pes.^2)+sum(0.5*P.a.*m_pesc.^2);
    pbv=abs(sum(m_pesc)); lp=max(P.PTDF'*m_pesc-P.F_l,0); ln=max(-P.PTDF'*m_pesc-P.F_l,0); lv=max(lp,ln);
    re_hist(k)=100*abs(cost-ref)/refscale; cv_hist(k)=norm([pbv;lp;ln],2); pb_hist(k)=pbv; maxline_hist(k)=max(lv(finite_mask));
    line_metric(:,k)=lv(finite_mask)./P.F_l(finite_mask);
    ilp=max(P.PTDF'*pesc-P.F_l,0); iln=max(-P.PTDF'*pesc-P.F_l,0); line_inst(:,k)=max(ilp,iln); line_inst(:,k)=line_inst(:,k)./P.F_l(finite_mask);
    comp_hist(k)=max([pbv/P.system_balance_scale_kW;line_metric(:,k)]);

    % 外层梯度始终使用本轮瞬时 raw Agg 上传；平均值仅用于 RE/CV 输出诊断。
    nextAuxPB=piPB-P.alpha_PB*sum(pesc); piPB=nextAuxPB+P.outer_momentum*(nextAuxPB-auxPB); auxPB=nextAuxPB;
    nextAuxPos=min(piLpos-P.alpha_l*(P.PTDF'*pesc-P.F_l),0); nextAuxNeg=min(piLneg-P.alpha_l*(-P.PTDF'*pesc-P.F_l),0);
    piLpos=nextAuxPos+P.outer_momentum*(nextAuxPos-auxLpos); piLneg=nextAuxNeg+P.outer_momentum*(nextAuxNeg-auxLneg); auxLpos=nextAuxPos; auxLneg=nextAuxNeg;
    pi0_next=piPB+P.PTDF*(piLpos-piLneg); price_hist(k)=norm(pi0_next-pi0,inf); pi0=pi0_next;

    start_stop=P.min_outer_iter;
    if use_ergodic_average, start_stop=max(start_stop,P.outer_avg_start_iter); end
    cv_ok=~logical(P.outer_cv_stop_enabled) || (pbv<=P.outer_pb_cv_tol_kW && all(lv(finite_mask)<=P.outer_line_cv_tol_kW(finite_mask)));
    if k>=start_stop && price_hist(k)<=P.errTol_UESM && cv_ok
        if P.outer_cv_stop_enabled, stop_reason='price_step_and_selected_cv'; else, stop_reason='price_step_only'; end
        break;
    end
    if isfinite(P.max_outer_iter) && k>=P.max_outer_iter
        stop_reason='max_outer_iter'; break;
    end
end

outer=struct;
outer.k=(1:k)'; outer.iter=outer.k; outer.totcost_ref=ref; outer.RE_totcost=re_hist(1:k); outer.CV_l2_history_kW=cv_hist(1:k);
outer.PB_history_kW=pb_hist(1:k); outer.max_finite_line_cv_history_kW=maxline_hist(1:k); outer.componentwise_CV_history=comp_hist(1:k);
outer.outer_history_semantics=sprintf('%s: %s raw Agg output; outer price update always uses instantaneous raw Agg.',scenario_name,ternary(use_ergodic_average,'suffix-average','instantaneous'));
outer.final_RE_totcost=outer.RE_totcost(end); outer.CV_kW=outer.CV_l2_history_kW(end); outer.PB_kW=pb_hist(k); outer.line_kW=max([0;lv]);
outer.system_balance_scale_kW=P.system_balance_scale_kW; outer.cv_policy=P.cv_policy; outer.finite_line_idx=finite_idx;
outer.line_ratio_instantaneous=line_inst(:,1:k); outer.line_ratio_average=line_metric(:,1:k); outer.outer_price_step=price_hist(1:k);
outer.outer_avg_cnt=avg_cnt; outer.outer_avg_start_iter=P.outer_avg_start_iter; outer.inner_avg_start_iter=P.inner_avg_start_iter;
outer.stop_reason=stop_reason; outer.errTol_UESM=P.errTol_UESM; outer.min_outer_iter=P.min_outer_iter; outer.max_outer_iter=P.max_outer_iter;
outer.pes=m_pes; outer.pg=m_pg; outer.pb=m_pb; outer.ps=m_ps; outer.pesc=m_pesc; outer.pi=pi; outer.pi_0=pi0;
outer.final_gap_kW=zeros(P.num_LESMs,1); for i=1:P.num_LESMs, outer.final_gap_kW(i)=abs(m_pesc(i)-sum(m_pes(P.start_idx(i):P.end_idx(i)))); end
outer.final_gap_mean_kW=mean(outer.final_gap_kW); outer.final_gap_max_kW=max(outer.final_gap_kW);
outer.raw_l2_cv_diagnostic_only=true; outer.random_seed=P.rng_seeds.outer_schedule; outer.num_local_call=num_local;
outer.total_inner_iter=total_inner_wall; outer.total_inner_price_updates=total_inner_updates; outer.max_inner_iter=max_wall; outer.max_inner_price_updates=max_updates;
outer.total_local_time=total_time; outer.avg_inner_iter=total_inner_wall/max(num_local,1); outer.avg_inner_price_updates=total_inner_updates/max(num_local,1);
outer.scenario=scenario_name; outer.outer_async_hold=outer_async_hold; outer.use_ergodic_average=use_ergodic_average;

save(sprintf('inner_k_data_%d.mat',scenario_id),'inner_k','params_snapshot');
save(sprintf('outer_data_%d.mat',scenario_id),'outer','params_snapshot');
fprintf('[%s finished] k=%d, RE=%.6f%%, total CV=%.8f p.u., reason=%s\n',scenario_name,k,outer.final_RE_totcost,outer.CV_kW/P.system_balance_scale_kW,stop_reason);
end

function s=ternary(condition,true_value,false_value)
if condition, s=true_value; else, s=false_value; end
end
