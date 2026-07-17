%% S0--S3：统一口径的外层 RE 与总 CV 对比图
% 运行前请依次运行 upperMarketS0、upperMarketS1、upperMarketS2、upperMarketS3。
% 本脚本只读取 outer_data_0.mat 至 outer_data_3.mat 并显示两张图，不保存任何图片。

clearvars;
close all;

data_dir = fileparts(mfilename('fullpath'));
scenario_names = {'S0','S1','S2','S3'};
scenario_desc = {'同步等待，瞬时输出', '异步保持，瞬时输出', ...
                 '同步等待，遍历平均', '异步保持，遍历平均'};
% 高对比、色盲友好：黑、朱红、蓝、蓝绿（Okabe--Ito 调色板）。
colors = [0.0000 0.0000 0.0000; ...  % S0: black
          0.8353 0.3686 0.0000; ...  % S1: vermilion
          0.0000 0.4471 0.6980; ...  % S2: blue
          0.0000 0.6196 0.4510];     % S3: bluish green
line_styles = {'-','--','-.',':'};
legend_labels = cell(1,4);
for s = 1:4
    legend_labels{s} = sprintf('%s — %s', scenario_names{s}, scenario_desc{s});
end

outer = cell(4,1);
missing_files = strings(0,1);
for s = 0:3
    file_name = fullfile(data_dir, sprintf('outer_data_%d.mat', s));
    if ~isfile(file_name)
        missing_files(end+1,1) = string(file_name); %#ok<SAGROW>
        continue;
    end
    S = load(file_name);
    assert(isfield(S,'outer'), '%s 中没有 outer 变量。', file_name);
    assert(isfield(S.outer,'RE_totcost') && isfield(S.outer,'CV_l2_history_kW') && ...
        isfield(S.outer,'system_balance_scale_kW'), ...
        '%s 不具备统一绘图字段；请用当前版本 upperMarketS%d 重新运行。', file_name, s);
    outer{s+1} = S.outer;
end
if ~isempty(missing_files)
    error('缺少结果文件。请先运行对应工况：\n%s', strjoin(missing_files,newline));
end

% 四个工况必须来自同一物理实例；此处检查总负荷基准是否一致。
load_scales = cellfun(@(x) x.system_balance_scale_kW, outer);
assert(all(isfinite(load_scales)) && all(load_scales > 0), '总负荷基准必须为正且有限。');
if max(abs(load_scales - load_scales(1))) > 1e-9 * load_scales(1)
    warning('四个 MAT 的总负荷基准不一致；CV p.u. 已分别按各自基准归一化，不宜作严格横向比较。');
end

%% 图 1：外层目标相对误差
figure('Name','S0--S3 outer RE','Color','w');
hold on;
for s = 1:4
    [k,re_percent] = getHistory(outer{s}, 'RE_totcost');
    plot(k, re_percent, 'LineWidth',1.5, 'LineStyle',line_styles{s}, ...
        'Color',colors(s,:));
end
grid on;
box on;
xlabel('Outer iteration k');
ylabel('RE (%)');
title('Outer relative objective error');
legend(legend_labels, 'Location','northeast');

%% 图 2：总约束违反的 p.u.
% 定义：||[PB; v^+; v^-]||_2 / sum_j D_j。
figure('Name','S0--S3 total CV','Color','w');
hold on;
for s = 1:4
    [k,cv_kW] = getHistory(outer{s}, 'CV_l2_history_kW');
    cv_pu = cv_kW ./ outer{s}.system_balance_scale_kW;
    plot(k, cv_pu, 'LineWidth',1.5, 'LineStyle',line_styles{s}, ...
        'Color',colors(s,:));
end
grid on;
box on;
xlabel('Outer iteration k');
ylabel('Total CV (p.u.)');
title('Total constraint violation normalized by total load');
legend(legend_labels, 'Location','northeast');

%% 命令行汇总（不写入文件）
fprintf('\n%-3s  %-9s  %-14s  %-16s\n', 'case','iterations','final RE (%)','final total CV (p.u.)');
for s = 1:4
    [k,re_percent] = getHistory(outer{s}, 'RE_totcost');
    [~,cv_kW] = getHistory(outer{s}, 'CV_l2_history_kW');
    fprintf('%-3s  %-9d  %-14.6f  %-16.8f\n', ...
        scenario_names{s}, numel(k), re_percent(end), cv_kW(end)/outer{s}.system_balance_scale_kW);
end

function [k,value] = getHistory(out, field_name)
value = out.(field_name)(:);
if isfield(out,'k') && numel(out.k) == numel(value)
    k = out.k(:);
else
    k = (1:numel(value)).';
end
assert(~isempty(value) && all(isfinite(value)), '字段 %s 为空或含 NaN/Inf。', field_name);
end
