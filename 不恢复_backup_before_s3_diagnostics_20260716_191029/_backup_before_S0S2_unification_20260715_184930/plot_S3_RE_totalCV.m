%% S3：外层 RE 与总约束违反（CV）绘图
% 直接在 F:\Comeon\不恢复 中运行本脚本。
% 只读取 MAT 文件和显示图窗；不会保存或覆盖任何图片/结果文件。

clearvars;
close all;

data_dir = fileparts(mfilename('fullpath'));
outer_file = fullfile(data_dir, 'outer_data_3.mat');
inner_file = fullfile(data_dir, 'inner_k_data_3.mat');

assert(isfile(outer_file), '找不到文件：%s', outer_file);
S_outer = load(outer_file);
assert(isfield(S_outer, 'outer'), 'outer_data_3.mat 中没有 outer 变量。');
outer = S_outer.outer;

required_fields = {'RE_totcost', 'CV_l2_history_kW', 'system_balance_scale_kW'};
for f = 1:numel(required_fields)
    assert(isfield(outer, required_fields{f}), ...
        'outer_data_3.mat 缺少字段：%s', required_fields{f});
end

re_outer_percent = outer.RE_totcost(:);       % 当前代码中已是百分数，例如 0.985 表示 0.985%%
cv_total_kW = outer.CV_l2_history_kW(:);      % ||[PB; v^+; v^-]||_2，单位 kW
load_scale_kW = outer.system_balance_scale_kW; % sum_j D_j；当前数据中 D_j 均为正

assert(isscalar(load_scale_kW) && isfinite(load_scale_kW) && load_scale_kW > 0, ...
    'system_balance_scale_kW 必须是正的有限标量。');
assert(numel(re_outer_percent) == numel(cv_total_kW), ...
    'RE_totcost 与 CV_l2_history_kW 的长度不一致。');

k_outer = (1:numel(re_outer_percent)).';
cv_total_pu = cv_total_kW / load_scale_kW;

figure('Name', 'S3 outer: RE and total CV', 'Color', 'w');
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot(k_outer, re_outer_percent, 'LineWidth', 1.25, 'Color', [0 0.4470 0.7410]);
grid on;
box on;
xlabel('Outer iteration k');
ylabel('RE (%)');
title('Outer relative objective error');

nexttile;
plot(k_outer, cv_total_pu, 'LineWidth', 1.25, 'Color', [0.8500 0.3250 0.0980]);
grid on;
box on;
xlabel('Outer iteration k');
ylabel('Total CV (p.u.)');
title('Total constraint violation normalized by total load');

fprintf('外层记录点数：%d\n', numel(k_outer));
fprintf('总负荷基准：%.6f kW\n', load_scale_kW);
fprintf('最终 RE：%.6f %%\n', re_outer_percent(end));
fprintf('最终总 CV：%.6f kW = %.8f p.u.\n', cv_total_kW(end), cv_total_pu(end));

% 若 parameters.m 开启了内层记录，inner_k_data_3.mat 会有数据；此时额外画内层诊断图。
if isfile(inner_file)
    S_inner = load(inner_file);
    if isfield(S_inner, 'inner_k') && isfield(S_inner.inner_k, 'h') && ~isempty(S_inner.inner_k.h)
        inner = S_inner.inner_k;
        assert(isfield(inner, 'RE_obj') && isfield(inner, 'val_violation'), ...
            'inner_k_data_3.mat 缺少 RE_obj 或 val_violation 字段。');

        h_inner = inner.h(:);
        re_inner_percent = inner.RE_obj(:);     % 当前代码中已是百分数
        cv_inner_kW = inner.val_violation(:);
        assert(numel(h_inner) == numel(re_inner_percent) && numel(h_inner) == numel(cv_inner_kW), ...
            'inner_k 的 h、RE_obj、val_violation 长度不一致。');

        figure('Name', 'S3 inner: RE and CV', 'Color', 'w');
        tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

        nexttile;
        plot(h_inner, re_inner_percent, 'LineWidth', 1.25, 'Color', [0.4660 0.6740 0.1880]);
        grid on;
        box on;
        xlabel('Inner iteration h');
        ylabel('RE (%)');
        title('Inner relative objective error');

        nexttile;
        plot(h_inner, cv_inner_kW, 'LineWidth', 1.25, 'Color', [0.4940 0.1840 0.5560]);
        grid on;
        box on;
        xlabel('Inner iteration h');
        ylabel('CV (kW)');
        title('Inner constraint violation');
    else
        fprintf('inner_k_data_3.mat 当前没有内层记录，未绘制内层图。\n');
    end
end
