function focus_table = run_S3_A2_fast_focus()
%RUN_S3_A2_FAST_FOCUS Recheck the deterministic screening focus set with fast QRE.

root = fileparts(mfilename('fullpath'));
cd(root);
screen = load('S3_A2_all_communities_screening.mat');
params = load('params.mat');
t = screen.summary_table;
num_comm = params.num_LESMs;
pick = [];
pick = [pick; topIds(t.stop_updates,10,'descend')];
pick = [pick; topIds(t.eu_count,10,'descend')];
pick = [pick; topIds(t.cv_tol,10,'ascend')];
distance_to_cap = abs(params.max_inner_iter - t.stop_updates);
pick = [pick; topIds(distance_to_cap,10,'ascend')];
pick = [pick; [1;3;49;71;123]];
pick = unique(pick(isfinite(pick) & pick >= 1 & pick <= num_comm));

rng_config = struct('inner_delay_seed',params.rng_seeds.inner_delay, ...
    'qre_noise_seed',params.rng_seeds.qre_noise, ...
    'qre_audit_seed',params.rng_seeds.qre_audit);
avg_config = struct('policy','dynamic_stable_start', ...
    'start_cv_factor',5,'start_price_factor',10,'start_stable_window',20, ...
    'min_price_updates',100,'min_samples',200, ...
    'formal_cv_stable_window',5);

m = numel(pick);
community_id = pick(:);
eu_count = zeros(m,1); cv_tol = zeros(m,1); avg_start = nan(m,1);
stop_updates = nan(m,1); avg_samples = nan(m,1); formal_cv = nan(m,1);
held_cv = nan(m,1); lambda_step = nan(m,1); max_age = nan(m,1);
runtime_s = nan(m,1); success = false(m,1); formal_stable = false(m,1);
audit_pass = false(m,1); stop_reason = repmat({''},m,1);
delay_hash = repmat({''},m,1); update_hash = repmat({''},m,1);
qre_hash = repmat({''},m,1); error_message = repmat({''},m,1);

for q = 1:m
    i = community_id(q);
    eu_count(q) = t.eu_count(i); cv_tol(q) = t.cv_tol(i);
    t0 = tic;
    try
        res = callLocal(params,i,rng_config,avg_config);
        runtime_s(q) = toc(t0);
        avg_start(q) = res.inner_avg_start_price_update;
        stop_updates(q) = res.price_updates;
        avg_samples(q) = res.inner_avg_cnt;
        formal_cv(q) = res.final_formal_cv_kW;
        held_cv(q) = res.final_held_cv_kW;
        lambda_step(q) = res.final_lambda_step;
        max_age(q) = res.max_observed_EU_age;
        success(q) = logical(res.success);
        formal_stable(q) = logical(res.formal_cv_stable_pass);
        audit_pass(q) = logical(res.qre_all_certificate_pass);
        stop_reason{q} = res.stop_reason;
        delay_hash{q} = res.delay_sequence_hash;
        update_hash{q} = res.update_sequence_hash;
        qre_hash{q} = res.qre_noise_sequence_hash;
    catch ME
        runtime_s(q) = toc(t0);
        error_message{q} = ME.message;
        stop_reason{q} = ['exception: ',ME.identifier];
    end
    fprintf('[A2 fast focus] %3d/%3d, community=%3d, stop=%7.0f, CV=%10.4g, ok=%d\n', ...
        q,m,i,stop_updates(q),formal_cv(q),success(q));
end

focus_table = table(community_id,eu_count,cv_tol,avg_start,stop_updates, ...
    avg_samples,formal_cv,held_cv,lambda_step,max_age,runtime_s,success, ...
    formal_stable,audit_pass,stop_reason,delay_hash,update_hash,qre_hash,error_message);
meta = struct('focus_ids',community_id,'qre_noise_enabled',true, ...
    'qre_audit_enabled',true,'qre_audit_rate',1, ...
    'seed_derivation_policy','deriveDeterministicSeed(base_seed,community_id,local_call_count,stream_tag)');
save('S3_A2_fast_focus.mat','focus_table','meta');
writetable(focus_table,'S3_A2_fast_focus.csv');
disp(focus_table);
fprintf('fast focus success=%d/%d, max stop=%g, max runtime=%.3fs\n', ...
    sum(success),m,max(stop_updates),max(runtime_s));

    function ids = topIds(value,n,order)
        ids_all = (1:numel(value))';
        [~,ord] = sort(value,order,'MissingPlacement','last');
        ids = ids_all(ord(1:min(n,numel(ord))));
    end

    function res = callLocal(p,community,rng_cfg,avg_cfg)
        s = p.start_idx(community):p.end_idx(community);
        res = localMarket_S3(p.c(s),p.b(s),p.a(community),p.D(s), ...
            p.init_rho(community),p.pg_max(s),p.pg_min(s),p.init_pi(community), ...
            p.init_pes(s),p.init_pg(s),p.init_pb(s),p.init_ps(s), ...
            p.init_pesc(community),p.errTol_LESMs(community), ...
            p.inner_cv_tol_kW(community),p.init_pi_0(community),p.pi_max,p.pi_min, ...
            false,p.beta_qre(s),p.qre_epsilon,p.agg_epsilon_i(community), ...
            p.qre_z_cap,p.qre_backoff_factor,p.qre_max_backoffs,p.inner_momentum, ...
            p.max_inner_iter,p.h0_max,p.min_inner_iter,p.inner_avg_start_iter, ...
            p.inner_cv_stop_enabled,p.qre_certificate_enabled,p.agg_certificate_enabled, ...
            p.agg_cert_tol,p.stable_inner_window,true,1, ...
            p.qre_audit_trace_community_full,p.qre_audit_seed, ...
            p.agg_gap_diagnostic_enabled,true,p.diagnostic_record_every, ...
            p.exact_sync_diagnostic_every,p.rolling_window,false,inf, ...
            1,community,rng_cfg,avg_cfg,1);
    end
end
