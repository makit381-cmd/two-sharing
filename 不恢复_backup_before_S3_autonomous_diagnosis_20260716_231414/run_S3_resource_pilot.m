function pilot_summary = run_S3_resource_pilot(pilot_outer_price_updates)
% 执行 S3 资源试跑；正式 params.mat 不被修改。
% 运行档位：P1/P5/P20。脚本 upperMarketS3.m 通过环境变量读取本次 cap。

if nargin < 1 || isempty(pilot_outer_price_updates)
    pilot_outer_price_updates = 1;
end
validateattributes(pilot_outer_price_updates,{'numeric'}, ...
    {'scalar','real','finite','positive','integer'});
P = load('params.mat');
run_tag = sprintf('pilotP%02d_%s',pilot_outer_price_updates,datestr(now,'yyyymmdd_HHMMSS'));
checkpoint_dir = fullfile(pwd,['checkpoints_' run_tag]);
if ~isfolder(checkpoint_dir), mkdir(checkpoint_dir); end
checkpoint_every = min(max(pilot_outer_price_updates,1),5);
formal_parameters = struct( ...
    'outer_avg_start_price_update',P.outer_avg_start_price_update, ...
    'max_outer_price_updates',P.max_outer_price_updates, ...
    'checkpoint_every_outer_updates',P.checkpoint_every_outer_updates);
effective_parameters = struct( ...
    'outer_avg_start_price_update',P.outer_avg_start_price_update, ...
    'max_outer_price_updates',pilot_outer_price_updates, ...
    'checkpoint_every_outer_updates',checkpoint_every);
test_overrides = struct( ...
    'outer_price_update_cap',pilot_outer_price_updates, ...
    'checkpoint_every_outer_updates',checkpoint_every, ...
    'label','RESOURCE PILOT ONLY; NOT FORMAL PROFILE');

setenv('HDEM_QRE_S3_MAX_OUTER_UPDATES',num2str(pilot_outer_price_updates));
setenv('HDEM_QRE_S3_CHECKPOINT_EVERY',num2str(checkpoint_every));
setenv('HDEM_QRE_S3_CHECKPOINT_DIR',checkpoint_dir);
setenv('HDEM_QRE_S3_RUN_TAG',run_tag);
setenv('HDEM_QRE_S3_RESUME','0');
setenv('HDEM_QRE_S3_CHECKPOINT','');
cleanup_env = onCleanup(@clearPilotEnvironment); %#ok<NASGU>

mem_before = captureMemory();
run_tic = tic;
script_path = strrep(fullfile(pwd,'upperMarketS3.m'),'''','''''');
evalin('base',sprintf('run(''%s'')',script_path));
total_runtime_s = toc(run_tic);
mem_after = captureMemory();

if ~isfile('outer_data_3.mat')
    error('run_S3_resource_pilot:missingResult','upperMarketS3 did not produce outer_data_3.mat.');
end
loaded = load('outer_data_3.mat','outer');
outer = loaded.outer;
pilot_outer_file = sprintf('S3_resource_pilot_P%02d_outer.mat',pilot_outer_price_updates);
save(pilot_outer_file,'outer','formal_parameters','effective_parameters','test_overrides');
checkpoint_files = dir(fullfile(checkpoint_dir,'S3_checkpoint_*.mat'));
checkpoint_bytes = sum([checkpoint_files.bytes]);
local_runtime = outer.local_call_runtime_history(:);
outer_runtime = outer.outer_update_runtime_history(:);
pilot_summary = struct();
pilot_summary.profile = sprintf('P%d',pilot_outer_price_updates);
pilot_summary.total_runtime_s = total_runtime_s;
pilot_summary.outer_price_updates = outer.outer_price_update_count;
pilot_summary.outer_wall_iter = outer.outer_wall_iter;
pilot_summary.num_local_call = outer.num_local_call;
pilot_summary.total_inner_price_updates = outer.total_inner_price_updates;
pilot_summary.total_inner_wall_iter = outer.total_inner_wall_iter;
pilot_summary.local_call_runtime_p50_s = percentileOrNaN(local_runtime,50);
pilot_summary.local_call_runtime_p90_s = percentileOrNaN(local_runtime,90);
pilot_summary.local_call_runtime_p95_s = percentileOrNaN(local_runtime,95);
pilot_summary.local_call_runtime_max_s = maxOrNaN(local_runtime);
pilot_summary.outer_update_runtime_history_s = outer_runtime;
pilot_summary.outer_update_runtime_p50_s = percentileOrNaN(outer_runtime,50);
pilot_summary.outer_update_runtime_max_s = maxOrNaN(outer_runtime);
pilot_summary.pb_cv_kW = outer.PB_history_kW(:);
pilot_summary.max_finite_line_cv_kW = outer.max_finite_line_cv_history_kW(:);
pilot_summary.outer_price_step = outer.outer_price_step(:);
pilot_summary.final_outer_avg_cnt = outer.outer_avg_cnt;
pilot_summary.checkpoint_dir = checkpoint_dir;
pilot_summary.checkpoint_count = numel(checkpoint_files);
pilot_summary.checkpoint_bytes = checkpoint_bytes;
pilot_summary.memory_before = mem_before;
pilot_summary.memory_after = mem_after;
pilot_summary.formal_parameters = formal_parameters;
pilot_summary.effective_parameters = effective_parameters;
pilot_summary.test_overrides = test_overrides;
pilot_summary.result_file = pilot_outer_file;
pilot_summary.success_gate = outer.outer_price_update_count == pilot_outer_price_updates && ...
    P.num_LESMs == 123 && outer.num_local_call > 0 && ...
    all(isfinite(outer.PB_history_kW));
summary_file = sprintf('S3_resource_pilot_P%02d_summary.mat',pilot_outer_price_updates);
save(summary_file,'pilot_summary');
fprintf('S3 resource pilot P%d: %s, %.2f s, %d local calls, %d checkpoints\n', ...
    pilot_outer_price_updates,ternary(pilot_summary.success_gate,'PASS','FAIL'), ...
    total_runtime_s,outer.num_local_call,numel(checkpoint_files));

    function clearPilotEnvironment()
        setenv('HDEM_QRE_S3_MAX_OUTER_UPDATES','');
        setenv('HDEM_QRE_S3_CHECKPOINT_EVERY','');
        setenv('HDEM_QRE_S3_CHECKPOINT_DIR','');
        setenv('HDEM_QRE_S3_RUN_TAG','');
        setenv('HDEM_QRE_S3_RESUME','');
        setenv('HDEM_QRE_S3_CHECKPOINT','');
    end
    function out = captureMemory()
        out = struct();
        try
            [user_view,system_view] = memory;
            out.MemUsedMATLAB = user_view.MemUsedMATLAB;
            out.MaxPossibleArrayBytes = user_view.MaxPossibleArrayBytes;
            out.VirtualAddressSpace = system_view.VirtualAddressSpace;
        catch
            out = struct('unavailable',true);
        end
    end
    function out = percentileOrNaN(v,p)
        if isempty(v), out = NaN; else, out = prctile(v,p); end
    end
    function out = maxOrNaN(v)
        if isempty(v), out = NaN; else, out = max(v); end
    end
    function out = ternary(condition,yes_value,no_value)
        if condition, out = yes_value; else, out = no_value; end
    end
end
