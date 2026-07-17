function summary = run_S3_checkpoint_resume_test(checkpoint_file)
% 用已有 P1 checkpoint 恢复到第 2 次真实外层更新，验证不重复 local-call count。

if nargin < 1 || isempty(checkpoint_file)
    files = dir(fullfile(pwd,'checkpoints_pilotP01_*','S3_checkpoint_*.mat'));
    if isempty(files)
        error('run_S3_checkpoint_resume_test:noCheckpoint', ...
            '请先运行 run_S3_resource_pilot(1)。');
    end
    [~,idx] = max([files.datenum]);
    checkpoint_file = fullfile(files(idx).folder,files(idx).name);
end
checkpoint = load(checkpoint_file,'outer_price_update_count','community_local_call_count', ...
    'checkpoint_next_k','run_tag');
P = load('params.mat');
if checkpoint.outer_price_update_count ~= 1 || checkpoint.checkpoint_next_k ~= 2
    error('run_S3_checkpoint_resume_test:invalidCheckpoint', ...
        '输入 checkpoint 不是 update=1 的可恢复状态。');
end

run_tag = ['resume_' datestr(now,'yyyymmdd_HHMMSS')];
setenv('HDEM_QRE_S3_MAX_OUTER_UPDATES','2');
setenv('HDEM_QRE_S3_CHECKPOINT_EVERY','2');
setenv('HDEM_QRE_S3_CHECKPOINT_DIR',fileparts(checkpoint_file));
setenv('HDEM_QRE_S3_RUN_TAG',run_tag);
setenv('HDEM_QRE_S3_RESUME','1');
setenv('HDEM_QRE_S3_CHECKPOINT',checkpoint_file);
cleanup_env = onCleanup(@clearEnvironment); %#ok<NASGU>
script_path = strrep(fullfile(pwd,'upperMarketS3.m'),'''','''''');
evalin('base',sprintf('run(''%s'')',script_path));
loaded = load('outer_data_3.mat','outer');
outer = loaded.outer;
resume_call_counts = outer.community_local_call_count;
call_count_delta = resume_call_counts - checkpoint.community_local_call_count;
new_local_calls = outer.num_local_call - sum(checkpoint.community_local_call_count);
summary = struct();
summary.checkpoint_file = checkpoint_file;
summary.checkpoint_outer_price_update_count = checkpoint.outer_price_update_count;
summary.resumed_outer_price_update_count = outer.outer_price_update_count;
summary.checkpoint_next_k = checkpoint.checkpoint_next_k;
summary.community_local_call_count_before = checkpoint.community_local_call_count;
summary.community_local_call_count_after = resume_call_counts;
summary.no_duplicate_call_count = all(call_count_delta >= 0 & call_count_delta <= 1) && ...
    sum(call_count_delta) == new_local_calls;
summary.call_count_delta = call_count_delta;
summary.price_update_progress_pass = outer.outer_price_update_count == 2;
summary.success_gate = summary.no_duplicate_call_count && summary.price_update_progress_pass;
summary.formal_parameters = struct('outer_avg_start_price_update',P.outer_avg_start_price_update, ...
    'max_outer_price_updates',P.max_outer_price_updates);
summary.effective_parameters = struct('outer_avg_start_price_update',P.outer_avg_start_price_update, ...
    'max_outer_price_updates',2);
summary.test_overrides = struct('outer_price_update_cap',2, ...
    'resume_from_checkpoint',true,'label','CHECKPOINT TEST ONLY; NOT FORMAL PROFILE');
save('S3_checkpoint_resume_test.mat','summary');
fprintf('S3 checkpoint/resume test: %s\n',ternary(summary.success_gate,'PASS','FAIL'));

    function clearEnvironment()
        setenv('HDEM_QRE_S3_MAX_OUTER_UPDATES','');
        setenv('HDEM_QRE_S3_CHECKPOINT_EVERY','');
        setenv('HDEM_QRE_S3_CHECKPOINT_DIR','');
        setenv('HDEM_QRE_S3_RUN_TAG','');
        setenv('HDEM_QRE_S3_RESUME','');
        setenv('HDEM_QRE_S3_CHECKPOINT','');
    end
    function out = ternary(condition,yes_value,no_value)
        if condition, out = yes_value; else, out = no_value; end
    end
end
