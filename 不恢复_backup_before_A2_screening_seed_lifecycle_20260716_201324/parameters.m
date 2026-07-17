clc;clear;
seed_cost_load = 112;
seed_community_cost = 123;
seed_pg_qre = 12;
seed_final_state = 115;
rng_seeds = struct( ...
    'cost_load', seed_cost_load, ...
    'community_cost', seed_community_cost, ...
    'pg_qre', seed_pg_qre, ...
    'final_state', seed_final_state, ...
    'outer_schedule', 999, ...
    'inner_delay', 31001, ...
    'qre_noise', 31002, ...
    'qre_audit', 31003, ...
    'outer_delay', 31004);
%将11250个产消者分配给123个社区
u = [85, 85, 50, 85, 85, 50, 70, 50, 100, 85, ...
     85, 120, 50, 50, 50, 120, 120, 50, 70, 120, ...
     80, 80, 90, 90, 90, 50, 120, 110, 80, 78, ...
     80, 80, 80, 85, 70, 80, 80, 58, 80, 69, ...
     50, 85, 50, 80, 80, 80, 110, 80, 400, 68, ...
     90, 90, 90, 90, 90, 90, 90, 90, 90, 50, ...
     70, 90, 120, 90, 420, 100, 50, 80, 80, 80, ...
     80, 90, 80, 200, 80, 377, 80, 80, 90, 60, ...
     85, 90, 85, 80, 80, 80, 80, 120, 80, 80, ...
     70, 90, 90, 80, 80, 50, 50, 120, 120, 120, ...
     70, 70, 70, 70, 110, 80, 80, 100, 100, 80, ...
     80, 90, 80, 80, 50, 100, 100, 80, 120, 90, ...
     90, 90, 90];
%定义常量参数
pi_max = 0.2;
pi_min = 0.05;

num_LESMs = 123;
num_prosumers = 11250;
% pg_ranges = [35 50; 20 35; 15 25; 5 10; 0 5];
% pg_max = zeros(num_prosumers,1);
% pg_min = zeros(num_prosumers,1);

%%
c = zeros(num_prosumers,1);
b = zeros(num_prosumers,1);                                                                                                                  
D = zeros(num_prosumers,1);
a = zeros(num_LESMs,1);
%power flow
PTDF = PTDF_matrix();

F_l = ones(num_LESMs-1,1) .* inf;%网络约束congestion       ---先不设置网络约束  122*1                                                                                     
congested_lines_idx = [
    find_line_index(10,63);
    find_line_index(35,36);
    find_line_index(4,74);
    find_line_index(8,9);
    find_line_index(15,16);
    find_line_index(78, 79);
    find_line_index(95, 96);
    ];
% line_list = [61,77,13,55,101,43,24];
lines_capacity = [0.5,8,1.5,21,9,10,15] .* 1e3;
F_l(congested_lines_idx) = lines_capacity;
end_idx = cumsum(u);
start_idx = [1,end_idx(1:end-1)+1];


init_pes = zeros(num_prosumers,1);
%% 
% init_pes = initial_pes;
% init_pb = initial_pb;
% init_ps = initial_ps;
% init_pi = initial_pi;
% init_pi_0 = initial_pi_0;
% init_pi_PB = 0.1769;
% init_pi_l_pos = initial_pi_l_pos;
% init_pi_l_neg = initial_pi_l_neg;
%%
init_pesc = zeros(num_LESMs,1);

% 迭代与 QRE 参数：所有场景从此文件生成正式 params.mat。

qre_epsilon = 3e-2;    % 单个 EU 局部拉格朗日目标的绝对误差上界
qre_z_cap = 3;
qre_backoff_factor = 0.5;
qre_max_backoffs = 16; % 固定随机方向下最多减半次数；失败则回退精确解
run_profile = getenv('HDEM_QRE_RUN_PROFILE');
if isempty(run_profile)
    run_profile = 'formal_audit';
end
qre_noise_enabled = true;
qre_audit_enabled = true;
qre_audit_rate = 0.01;
qre_audit_trace_community_full = true;
qre_audit_seed = 20260716;
diagnostic_record_every = 10;
exact_sync_diagnostic_every = 100;
rolling_window = 200;
agg_epsilon_margin = 3e-2;
agg_epsilon_i = u(:) .* qre_epsilon + agg_epsilon_margin;
agg_epsilon = agg_epsilon_i; % 兼容旧入口；正式语义是 123 维社区预算
agg_cert_tol = 1e-10;
qre_certificate_enabled = strcmp(run_profile,'certified_validation');
agg_certificate_enabled = false;
agg_gap_diagnostic_enabled = true;
qre_error_rel = qre_epsilon; % 仅供尚未迁移的 S0--S2 遗留入口读取；正式 S3 不使用该名称
% 理论步长版本关闭额外惯性项。
outer_momentum = 0;
inner_momentum = 0;

stable_inner_window = 5;
stable_outer_window = 5;
max_inner_iter = 20000; % 动态后缀平均后的有限安全保护，不是理论收敛条件
min_inner_iter = 20;
inner_cv_ratio = 0.01;
inner_price_scale_kW = 5;
inner_cv_stop_enabled = true;
% 动态稳定触发的工程性内层后缀平均；inner_avg_start_iter 仅保留兼容字段。
inner_avg_start_iter = 2;
inner_avg_policy = 'dynamic_stable_start';
inner_avg_start_cv_factor = 5;
inner_avg_start_price_factor = 10;
inner_avg_start_stable_window = 20;
inner_avg_min_price_updates = 100;
inner_avg_min_samples = 200;
max_outer_iter = 3000; % 仅安全保护
min_outer_iter = 100;
outer_cv_ratio = 0.01;
outer_pb_cv_ratio = 0.01;
outer_line_cv_ratio = 0.01;
outer_cv_stop_enabled = true;
outer_cv_l2_stop_enabled = false;
outer_cv_l2_tol_kW = 200;
% 工程性外层后缀平均：前 39 个外层更新仅作 burn-in，k=40 起统计上传输出。
outer_avg_start_iter = 100;
% 本轮不记录内层轨迹，避免无关的 RE/CV 计算与存档。
record_inner_comm = NaN;
record_k = [];
progress_print_enabled = true;
progress_print_inner_every = 1000;
progress_print_outer_every = 1;
% 每次 S3 仅覆盖固定结果文件，不再以时间戳新增历史 MAT。
save_s3_results = true;
run_tag = ['s3_kw_rawagg_inner50_price1em6_outer_history_', char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'))];
result_unit = 'kW';
rng(seed_cost_load, 'twister');
for j = 1:num_prosumers
    % c_min = 0.5e-3;
    % c_max = 1e-3;
    c_min = 2e-3;
    c_max = 6e-3;
    % b_min = 0.01;
    % b_max = 0.05;
    b_min = 0.03;
    b_max = 0.08;
    D_min = 0;
    D_max = 40;

    c(j) = c_min + (c_max - c_min) .* rand(1);
    b(j) = b_min + (b_max - b_min) .* rand(1);
    D(j) = D_min + (D_max - D_min) .* rand(1);
end

% 记录当前物理实例的社区负荷尺度，单位均为 kW。
community_load_scale_kW = zeros(num_LESMs,1);
for i = 1:num_LESMs
    community_load_scale_kW(i) = sum(abs(D(start_idx(i):end_idx(i))));
end
if strcmp(run_profile,'certified_validation')
    inner_cv_policy = 'absolute_5kW';
    inner_cv_tol_kW = 5 * ones(num_LESMs,1);
else
    inner_cv_policy = 'relative_1pct_community_load';
    inner_cv_tol_kW = inner_cv_ratio .* community_load_scale_kW;
end
% 外层诊断尺度：总平衡相对系统绝对负荷；线路违反相对各自有限容量。
system_balance_scale_kW = sum(abs(D));
outer_pb_cv_tol_kW = outer_pb_cv_ratio .* system_balance_scale_kW;
outer_line_cv_tol_kW = outer_line_cv_ratio .* F_l;
cv_policy = 'relative_componentwise_pb_and_finite_lines';

a_extend = zeros(num_prosumers,1);
init_xj2 = zeros(num_LESMs,1);
rng(seed_community_cost, 'twister');
for i = 1:num_LESMs
    init_pesc(i) = sum(init_pes(start_idx(i):end_idx(i)));
    a_min = (2.5e-3 / u(i));
    a_max = (5e-3 / u(i));
    a(i) = a_min + (a_max - a_min) .* rand(1);
    a_extend(start_idx(i):end_idx(i)) = a(i);
    init_xj2(i) = sum(init_pes(start_idx(i):end_idx(i)) .^ 2);
end

init_pi_l_pos = zeros(num_LESMs-1,1);
init_pi_l_neg = zeros(num_LESMs-1,1);
init_pi_PB = 0.1;
init_pi_0 = init_pi_PB + PTDF * (init_pi_l_pos - init_pi_l_neg);
%% 
% 定义 5 种最大发电功率区间 (单位: kW)
pg_ranges = [
    35, 50;  % Range 1: 超高发电 (盈余区主力)
    20, 35;  % Range 2: 高发电
    15, 25;  % Range 3: 中等 (平衡区主力)
     5, 10;  % Range 4: 低发电
     0,  5   % Range 5: 极低/不发电 (缺额区主力)
];

% 初始化变量
pg_max = zeros(num_prosumers, 1);
pg_min = zeros(num_prosumers, 1); 

% --- 步骤 A: 硬编码社区类型 (基于 Figure 3) ---
% Type 1 = 盈余 (Energy-surplus, 绿色区域)
% Type 2 = 平衡 (Energy-balance, 橙色区域)
% Type 3 = 缺额 (Energy-deficit, 蓝色区域)

% 初始化所有社区为 Type 2 (平衡)，然后覆盖盈余和缺额
community_types_list = 2 * ones(num_LESMs, 1); 

% === 1. 定义盈余区节点 (绿色区域) ===
% 从图中提取的节点编号
surplus_nodes = [ ...
    36,35,47,48,55,34,56,57,58,12,59,60,61,62,13,30,31,32,33,14,15,27,...
    28,29,16,25,26,17,19,18,20,24,21,22,23,93,90,91,92,75,76,77,89,78,88,...
    79,86,87,85,80,81,82,83,84,63,64,65,66,67
];
community_types_list(surplus_nodes) = 1;

% === 2. 定义缺额区节点 (蓝色区域) ===
% 从图中提取的节点编号
deficit_nodes = [ ...
    1,2,3,4,5,119,120,121,122,118,110,111,112,113,114,123,115,116,117,...
    42,41,43,40,39,44,38,45,37,46,95,106,96,101,104,105,102,103,97,98,99,100
];
community_types_list(deficit_nodes) = 3;

% === 3. 定义平衡区节点 (橙色区域) ===

% --- 步骤 B: 为每个产消者分配 pg_max (规则同前) ---
rng(seed_pg_qre, 'twister');

big = 0.35;
small = 0.2;
min_x_scale = 5;

rho_qre = zeros(num_prosumers, 1);
sigma_qre = zeros(num_prosumers, 1);
beta_qre = zeros(num_prosumers, 1);

for j = 1:num_prosumers

    comm_id = find(end_idx >= j, 1, 'first');
    c_type = community_types_list(comm_id);

    val_rand = rand(1);
    selected_range_idx = 0;

    if c_type == 1
        if val_rand < 0.8
            selected_range_idx = 1;
        else
            selected_range_idx = 2;
        end

    elseif c_type == 2
        selected_range_idx = randi(5);

    elseif c_type == 3
        if val_rand < 0.8
            selected_range_idx = 5;
        else
            selected_range_idx = 4;
        end
    end

    p_limits = pg_ranges(selected_range_idx, :);
    pg_max(j) = p_limits(1) + (p_limits(2) - p_limits(1)) * rand(1);
    % pg_max(j) = inf;
    pg_min(j) = 0;

    % rho_qre(j) = small + (big - small) * rand(1);
    % x_scale = max(pg_max(j) - pg_min(j), min_x_scale);
    % sigma_qre(j) = rho_qre(j) * x_scale;
    % beta_qre(j) = 1 / (a_extend(j) * sigma_qre(j)^2);
    beta_qre(j) = 2e3 + 8e5 * rand(1);
end

% 原冷启动：零交换、局部功率平衡且外层网络初始可行。
% 不使用随机 RE=30%% 数值初值。
init_mode = 'cold_feasible_original';
init_target_RE = NaN;
init_reference_obj = NaN;
init_pg_cold = min(max(D + init_pes, pg_min), pg_max);
init_pg = init_pg_cold;
init_pes = zeros(num_prosumers,1);
init_pb = max(D + init_pes - init_pg, 0);
init_ps = max(init_pg - D - init_pes, 0);
for i = 1:num_LESMs
    init_pesc(i) = sum(init_pes(start_idx(i):end_idx(i)));
    init_xj2(i) = sum(init_pes(start_idx(i):end_idx(i)).^2);
end
init_pi = init_pi_0 - a .* init_pesc;

%% 设置内外层最大延迟
n_comm = end_idx - start_idx + 1;
outer_max_delay = 12;
k0_max = max(min(outer_max_delay,ceil(0.5* sqrt(n_comm))),6);
h0_max = 6;
%% 内外层最大利普西茨常数
LD_in_i = (n_comm(:) + 1) ./ a;
LD_in = max(LD_in_i);
active_lines = isfinite(F_l);
% 只有容量有限的线路才属于实际外层不等式约束矩阵。
P_active = PTDF(:, active_lines);

c_out = a .* (1 + 1 ./ u(:));
LD_out = sum((1 + 2 * sum(P_active.^2, 2)) ./ c_out);

%% 内外层迭代步长
% 外层按与当前单向 held-information 假设对应的统一标量步长更新。
% PB 平衡块和正/负线路约束块均使用同一 alpha，直接对应 LD_out 的保守界。
k0_out = max(k0_max);
alpha_out_limit = min(1 / ((k0_out + 0.5) * LD_out), 1 / (4 * LD_out));
outer_step_safety = 0.8;

% 内层的 123 个社区是独立对偶块；按各自 L_D,i 取单向 held-information 理论步长。
% h0=6 仅描述 EU 原始变量 held，不使用论文双向异步的 2*h0+1/2 系数。
inner_step_safety = 0.90;
alpha_in_limit_i = min(1 ./ ((h0_max + 0.5) .* LD_in_i), 1 ./ (4 .* LD_in_i));
alpha_in = inner_step_safety .* alpha_in_limit_i;
alpha_in_limit_user = alpha_in_limit_i; % 遗留 MAT 字段名：保存各社区理论上界。
init_rho = alpha_in ./ a;
alpha_PB = outer_step_safety * alpha_out_limit;
alpha_l = alpha_PB;
%% 内外层容错
% 内层停止同时检查价格稳定、正式 CV、社区预算和 EU 逐次认证。
errTol_LESMs = inner_price_scale_kW .* alpha_in;
% 外层价格稳定窗口与正式平均输出 CV 共同决定停止。
errTol_UESM = 1e-7;

rng(seed_final_state, 'twister');

fprintf('LD_out = %.16e\n',LD_out);
fprintf('alpha_out_limit = %.16e\n',alpha_out_limit);
fprintf('alpha_PB = alpha_l = %.16e\n',alpha_PB);
fprintf('alpha_out ratio = %.6f\n',alpha_PB / alpha_out_limit);
fprintf('alpha_in range = [%.16e, %.16e]\n',min(alpha_in),max(alpha_in));
fprintf('LD_in range = [%.16e, %.16e]\n',min(LD_in_i),max(LD_in_i));
fprintf('k0_max range = [%d, %d]\n',min(k0_max),max(k0_max));
fprintf('beta_qre range = [%.16e, %.16e]\n',min(beta_qre),max(beta_qre));

assert(all(alpha_in > 0));
assert(all(alpha_in < alpha_in_limit_i));
assert(max(abs(init_rho - alpha_in ./ a)) <= 1e-14);
assert(alpha_PB > 0 && alpha_PB < alpha_out_limit);
assert(abs(alpha_l - alpha_PB) <= eps(max(1,abs(alpha_PB))));
assert(inner_momentum == 0 && outer_momentum == 0);
assert(qre_audit_enabled && agg_gap_diagnostic_enabled);
assert(numel(agg_epsilon_i) == num_LESMs && all(agg_epsilon_i > 0));

initmu_PB = zeros(num_prosumers, 1);
save('params.mat','a','b','c','D','start_idx','end_idx', ...
    'init_pes', 'init_pesc', 'init_pi_l_pos', 'init_pi_l_neg', ...
    'init_pb','init_ps','init_pg','init_pi','init_pi_0','init_pi_PB', 'init_xj2', ...
    'init_mode','init_target_RE','init_reference_obj','init_pg_cold', ...
    'num_LESMs','num_prosumers','pi_max','pi_min','u','F_l', ...
    'alpha_l','alpha_PB','PTDF','init_rho','init_pi_PB', ...
    'init_pi_0','init_pi','errTol_LESMs','errTol_UESM', ...
    'pg_max','pg_min','a_extend','beta_qre', ...
     'run_profile','qre_noise_enabled','qre_audit_enabled','qre_audit_rate','qre_audit_trace_community_full','qre_audit_seed', ...
     'diagnostic_record_every','exact_sync_diagnostic_every','rolling_window', ...
     'qre_epsilon','agg_epsilon','agg_epsilon_i','agg_epsilon_margin','agg_cert_tol','qre_z_cap','qre_backoff_factor','qre_max_backoffs','qre_certificate_enabled','agg_certificate_enabled','agg_gap_diagnostic_enabled','qre_error_rel','outer_momentum','inner_momentum','max_inner_iter', ...
    'min_inner_iter','inner_cv_ratio','inner_cv_policy','inner_cv_stop_enabled','stable_inner_window','stable_outer_window','community_load_scale_kW','inner_cv_tol_kW','inner_avg_start_iter','inner_avg_policy','inner_avg_start_cv_factor','inner_avg_start_price_factor','inner_avg_start_stable_window','inner_avg_min_price_updates','inner_avg_min_samples', ...
     'max_outer_iter','min_outer_iter','outer_cv_ratio','outer_pb_cv_ratio','outer_line_cv_ratio','outer_cv_stop_enabled','outer_cv_l2_stop_enabled','outer_cv_l2_tol_kW','system_balance_scale_kW','outer_pb_cv_tol_kW','outer_line_cv_tol_kW','cv_policy', ...
    'outer_avg_start_iter','record_inner_comm','record_k','progress_print_enabled','progress_print_inner_every','progress_print_outer_every','save_s3_results','run_tag','result_unit', ...
    'rho_qre','sigma_qre','small','big','min_x_scale', ...
    'rng_seeds', 'h0_max','k0_max','outer_max_delay','k0_out', ...
    'LD_in_i','LD_in','LD_out','active_lines','c_out', ...
    'alpha_in','alpha_in_limit_i','alpha_in_limit_user','inner_step_safety','inner_price_scale_kW','alpha_out_limit','outer_step_safety');
% 'init_pb','init_ps','init_pi','init_pi_0','init_pi_PB','init_pi_l_pos','init_pi_l_neg', ...
