function plot_four_scenarios_view_only()
% Display only: no figure, table, or directory is written to disk.
root_dir = fileparts(mfilename('fullpath'));
labels = categorical({'S0';'S1';'S2';'S3'},{'S0','S1','S2','S3'});
outer = cell(4,1);
final_k = zeros(4,1);
re = zeros(4,1);
cv = zeros(4,1);
success = false(4,1);

for s = 0:3
    file = fullfile(root_dir, sprintf('outer_data_%d.mat',s));
    assert(isfile(file), '缺少 %s，请先运行相应工况。', file);
    data = load(file,'outer');
    outer{s+1} = data.outer;
    final_k(s+1) = data.outer.iter(end);
    re(s+1) = data.outer.final_RE_totcost;
    cv(s+1) = data.outer.CV_MW;
    success(s+1) = logical(data.outer.state.success);
end

summary = table(labels,final_k,success,re,cv, ...
    'VariableNames',{'工况','最终迭代次数','外层成功','总成本相对误差_RE_百分比','最终约束违反_CV_MW'});
disp(summary);

figure('Name','四种工况对比（仅显示）','Color','w');
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile;
bar(labels,re); grid on;
ylabel('RE (%)'); title('最终总成本相对误差');
nexttile;
bar(labels,cv); grid on;
ylabel('CV (MW)'); title('最终约束违反');
nexttile; hold on;
for s = 1:4
    semilogy(outer{s}.iter,max(outer{s}.RE_totcost,eps),'LineWidth',1.1);
end
grid on; xlabel('外层迭代 k'); ylabel('RE (%)'); title('RE 轨迹');
legend({'S0','S1','S2','S3'},'Location','best');
nexttile; hold on;
for s = 1:4
    semilogy(outer{s}.iter,max(outer{s}.val_violation,eps),'LineWidth',1.1);
end
grid on; xlabel('外层迭代 k'); ylabel('CV (MW)'); title('CV 轨迹');
legend({'S0','S1','S2','S3'},'Location','best');
end
