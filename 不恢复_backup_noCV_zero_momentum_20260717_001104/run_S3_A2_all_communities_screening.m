function summary_table = run_S3_A2_all_communities_screening(first_id,last_id)
%RUN_S3_A2_ALL_COMMUNITIES_SCREENING Run only independent S3 inner calls.

root = fileparts(mfilename('fullpath'));
cd(root);
params = load('params.mat');
num_comm = params.num_LESMs;
if nargin < 1 || isempty(first_id), first_id = 1; end
if nargin < 2 || isempty(last_id), last_id = num_comm; end
validateattributes(first_id,{'numeric'},{'scalar','integer','>=',1,'<=',num_comm});
validateattributes(last_id,{'numeric'},{'scalar','integer','>=',first_id,'<=',num_comm});
partial_dir = fullfile(root,'S3_A2_all_communities_partial');
if ~exist(partial_dir,'dir'), mkdir(partial_dir); end
rng_config = struct('inner_delay_seed',params.rng_seeds.inner_delay, ...
    'qre_noise_seed',params.rng_seeds.qre_noise, ...
    'qre_audit_seed',params.rng_seeds.qre_audit);
avg_config = struct('policy','dynamic_stable_start', ...
    'start_cv_factor',5,'start_price_factor',10,'start_stable_window',20, ...
    'min_price_updates',100,'min_samples',200, ...
    'formal_cv_stable_window',5);

community_id = (1:num_comm)';
eu_count = params.end_idx(:) - params.start_idx(:) + 1;
load_kW = params.community_load_scale_kW(:);
cv_tol = params.inner_cv_tol_kW(:);
avg_start = nan(num_comm,1);
stop_updates = nan(num_comm,1);
price_update_count = nan(num_comm,1);
avg_samples = nan(num_comm,1);
formal_cv = nan(num_comm,1);
held_cv = nan(num_comm,1);
lambda_step = nan(num_comm,1);
max_age = nan(num_comm,1);
runtime_s = nan(num_comm,1);
success = false(num_comm,1);
formal_stable = false(num_comm,1);
no_nan_inf = false(num_comm,1);
audit_pass = false(num_comm,1);
stop_reason = repmat({''},num_comm,1);
delay_hash = repmat({''},num_comm,1);
update_hash = repmat({''},num_comm,1);
qre_hash = repmat({''},num_comm,1);
error_message = repmat({''},num_comm,1);

for i = first_id:last_id
    partial_file = fullfile(partial_dir,sprintf('community_%03d.mat',i));
    if exist(partial_file,'file')
        cached = load(partial_file,'one');
        one = cached.one;
        avg_start(i) = one.avg_start;
        stop_updates(i) = one.stop_updates;
        price_update_count(i) = one.price_update_count;
        avg_samples(i) = one.avg_samples;
        formal_cv(i) = one.formal_cv;
        held_cv(i) = one.held_cv;
        lambda_step(i) = one.lambda_step;
        max_age(i) = one.max_age;
        runtime_s(i) = one.runtime_s;
        success(i) = one.success;
        formal_stable(i) = one.formal_stable;
        no_nan_inf(i) = one.no_nan_inf;
        audit_pass(i) = one.audit_pass;
        stop_reason{i} = one.stop_reason;
        delay_hash{i} = one.delay_hash;
        update_hash{i} = one.update_hash;
        qre_hash{i} = one.qre_hash;
        error_message{i} = one.error_message;
        fprintf('[A2 deterministic cached] community=%3d, stop=%7.0f, CV=%10.4g, ok=%d\n', ...
            i,stop_updates(i),formal_cv(i),success(i));
        continue;
    end
    t0 = tic;
    try
        res = callLocal(params,i,rng_config,avg_config);
        runtime_s(i) = toc(t0);
        avg_start(i) = res.inner_avg_start_price_update;
        stop_updates(i) = res.price_updates;
        price_update_count(i) = res.price_update_count;
        avg_samples(i) = res.inner_avg_cnt;
        formal_cv(i) = res.final_formal_cv_kW;
        held_cv(i) = res.final_held_cv_kW;
        lambda_step(i) = res.final_lambda_step;
        max_age(i) = res.max_observed_EU_age;
        success(i) = logical(res.success);
        formal_stable(i) = logical(res.formal_cv_stable_pass);
        no_nan_inf(i) = all(isfinite([formal_cv(i),held_cv(i),lambda_step(i),max_age(i)]));
        audit_pass(i) = logical(res.qre_all_certificate_pass);
        stop_reason{i} = res.stop_reason;
        delay_hash{i} = res.delay_sequence_hash;
        update_hash{i} = res.update_sequence_hash;
        qre_hash{i} = res.qre_noise_sequence_hash;
    catch ME
        runtime_s(i) = toc(t0);
        error_message{i} = ME.message;
        stop_reason{i} = ['exception: ',ME.identifier];
    end
    one = struct('avg_start',avg_start(i),'stop_updates',stop_updates(i), ...
        'price_update_count',price_update_count(i),'avg_samples',avg_samples(i), ...
        'formal_cv',formal_cv(i),'held_cv',held_cv(i),'lambda_step',lambda_step(i), ...
        'max_age',max_age(i),'runtime_s',runtime_s(i),'success',success(i), ...
        'formal_stable',formal_stable(i),'no_nan_inf',no_nan_inf(i), ...
        'audit_pass',audit_pass(i),'stop_reason',stop_reason{i}, ...
        'delay_hash',delay_hash{i},'update_hash',update_hash{i}, ...
        'qre_hash',qre_hash{i},'error_message',error_message{i});
    save(partial_file,'one');
    fprintf('[A2 deterministic] %3d/%3d, EU=%3d, stop=%7.0f, CV=%10.4g, ok=%d\n', ...
        i,last_id,eu_count(i),stop_updates(i),formal_cv(i),success(i));
end

summary_table = table(community_id,eu_count,load_kW,cv_tol,avg_start, ...
    stop_updates,price_update_count,avg_samples,formal_cv,held_cv,lambda_step, ...
    max_age,runtime_s,success,formal_stable,no_nan_inf,audit_pass,stop_reason, ...
    delay_hash,update_hash,qre_hash,error_message);
meta = struct('profile','formal_audit','qre_noise_enabled',false, ...
    'qre_audit_enabled',true,'qre_audit_rate',1, ...
    'inner_avg_policy','dynamic_stable_start','inner_avg_start_cv_factor',5, ...
    'inner_avg_start_price_factor',10,'inner_avg_start_stable_window',20, ...
    'inner_avg_min_price_updates',100,'inner_avg_min_samples',200, ...
    'inner_formal_cv_stable_window',5,'inner_step_safety',params.inner_step_safety, ...
    'max_inner_iter',params.max_inner_iter,'community_local_call_count',ones(num_comm,1), ...
    'seed_derivation_policy','deriveDeterministicSeed(base_seed,community_id,local_call_count,stream_tag)');
batch_file = fullfile(root,sprintf('S3_A2_all_communities_screening_batch_%03d_%03d.mat',first_id,last_id));
save(batch_file,'summary_table','meta');
writetable(summary_table,fullfile(root,sprintf('S3_A2_all_communities_screening_batch_%03d_%03d.csv',first_id,last_id)));
if first_id == 1 && last_id == num_comm && all(isfinite(stop_updates))
    save('S3_A2_all_communities_screening.mat','summary_table','meta');
    writetable(summary_table,'S3_A2_all_communities_screening.csv');
    fprintf('deterministic success=%d/%d, max stop=%g, max runtime=%.3fs\n', ...
        sum(success),num_comm,max(stop_updates),max(runtime_s));
else
    fprintf('batch complete: communities %d-%d; rerun remaining batches, then call without arguments.\n', ...
        first_id,last_id);
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
            p.agg_gap_diagnostic_enabled,false,p.diagnostic_record_every, ...
            p.exact_sync_diagnostic_every,p.rolling_window,false,inf, ...
            1,community,rng_cfg,avg_cfg,1);
    end
end
