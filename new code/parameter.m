function parameter()
% Formal modified-b/c parameter generator using value_initial.mat as initial state.
root_dir = fileparts(mfilename('fullpath'));
addpath(root_dir);

seed_cost_load = 112;
seed_community_cost = 123;
seed_pg_qre = 12;
seed_final_state = 115;
rng_seeds = struct('cost_load',seed_cost_load,'community_cost',seed_community_cost, ...
    'pg_qre',seed_pg_qre,'final_state',seed_final_state);

u = [85,85,50,85,85,50,70,50,100,85,85,120,50,50,50,120,120,50,70,120, ...
    80,80,90,90,90,50,120,110,80,78,80,80,80,85,70,80,80,58,80,69, ...
    50,85,50,80,80,80,110,80,400,68,90,90,90,90,90,90,90,90,90,50, ...
    70,90,120,90,420,100,50,80,80,80,80,90,80,200,80,377,80,80,90,60, ...
    85,90,85,80,80,80,80,120,80,80,70,90,90,80,80,50,50,120,120,120, ...
    70,70,70,70,110,80,80,100,100,80,80,90,80,80,50,100,100,80,120,90, ...
    90,90,90];
num_LESMs = 123;
num_prosumers = 11250;
pi_max = 0.2;
pi_min = 0.05;
PTDF = PTDF_matrix();
alpha_PB = 5e-7;
alpha_l = 1e-6;
F_l = ones(num_LESMs - 1,1) .* 999e3;
congested_lines_idx = [find_line_index(10,63); find_line_index(35,36); ...
    find_line_index(4,74); find_line_index(8,9); find_line_index(15,16); ...
    find_line_index(78,79); find_line_index(95,96)];
F_l(congested_lines_idx) = [0.5,8,1.5,21,9,10,15]' .* 1e3;

end_idx = cumsum(u);
start_idx = [1,end_idx(1:end-1)+1];
c = zeros(num_prosumers,1);
b = zeros(num_prosumers,1);
D = zeros(num_prosumers,1);
a = zeros(num_LESMs,1);

% Current formal ranges selected for the modified-b/c model.
rng(seed_cost_load, 'twister');
for j = 1:num_prosumers
    c(j) = 2e-3 + (6e-3 - 2e-3) * rand();
    b(j) = 0.03 + (0.08 - 0.03) * rand();
    D(j) = 40 * rand();
end

a_extend = zeros(num_prosumers,1);
rng(seed_community_cost, 'twister');
for i = 1:num_LESMs
    a(i) = (2.5e-3 / u(i)) + (2.5e-3 / u(i)) * rand();
    a_extend(start_idx(i):end_idx(i)) = a(i);
end

pg_ranges = [35,50;20,35;15,25;5,10;0,5];
community_types = 2 * ones(num_LESMs,1);
surplus_nodes = [36,35,47,48,55,34,56,57,58,12,59,60,61,62,13,30,31,32,33,14,15,27, ...
    28,29,16,25,26,17,19,18,20,24,21,22,23,93,90,91,92,75,76,77,89,78,88, ...
    79,86,87,85,80,81,82,83,84,63,64,65,66,67];
deficit_nodes = [1,2,3,4,5,119,120,121,122,118,110,111,112,113,114,123,115,116,117, ...
    42,41,43,40,39,44,38,45,37,46,95,106,96,101,104,105,102,103,97,98,99,100];
community_types(surplus_nodes) = 1;
community_types(deficit_nodes) = 3;

pg_max = zeros(num_prosumers,1);
pg_min = zeros(num_prosumers,1);
beta_qre = zeros(num_prosumers,1);
rng(seed_pg_qre, 'twister');
for j = 1:num_prosumers
    i = find(end_idx >= j,1,'first');
    r = rand();
    if community_types(i) == 1
        range_id = 1 + (r >= 0.8);
    elseif community_types(i) == 3
        range_id = 5 - (r >= 0.8);
    else
        range_id = randi(5);
    end
    limits = pg_ranges(range_id,:);
    pg_max(j) = limits(1) + (limits(2) - limits(1)) * rand();
    pg_min(j) = 0;
    sigma_qre = (0.05 + 0.25 * rand()) * max(pg_max(j),5);
    beta_qre(j) = 1 / (a_extend(j) * sigma_qre^2);
end

value_file = fullfile(root_dir, 'value_initial.mat');
assert(isfile(value_file), '请先运行 initial_data.m 生成 value_initial.mat。');
value = load(value_file, 'initial_pg','initial_pes','initial_pb','initial_ps', ...
    'initial_pesc','initial_pi','initial_pi_0','initial_pi_PB', ...
    'initial_pi_l_pos','initial_pi_l_neg','preheat_metadata');

init_pg = value.initial_pg(:);
init_pes = value.initial_pes(:);
init_pb = value.initial_pb(:);
init_ps = value.initial_ps(:);
init_pesc = value.initial_pesc(:);
init_pi = value.initial_pi(:);
init_pi_0 = value.initial_pi_0(:);
init_pi_PB = value.initial_pi_PB;
init_pi_l_pos = value.initial_pi_l_pos(:);
init_pi_l_neg = value.initial_pi_l_neg(:);
init_pi_l = init_pi_l_pos;

assert(numel(init_pg) == num_prosumers && numel(init_pes) == num_prosumers && ...
    numel(init_pb) == num_prosumers && numel(init_ps) == num_prosumers, ...
    'value_initial 的产消者初值维度不正确。');
assert(numel(init_pesc) == num_LESMs && numel(init_pi) == num_LESMs && ...
    numel(init_pi_0) == num_LESMs, 'value_initial 的社区初值维度不正确。');
assert(norm(init_pi_0 - (init_pi_PB + PTDF * ...
    (init_pi_l_pos - init_pi_l_neg)), inf) <= 1e-10, ...
    'value_initial 的价格状态不一致。');
for i = 1:num_LESMs
    assert(abs(init_pesc(i) - sum(init_pes(start_idx(i):end_idx(i)))) <= 1e-10, ...
        'value_initial 的社区耦合状态不一致。');
end

init_rho = 0.4 * ones(num_LESMs,1);
errTol_LESMs = 5e-6 * ones(num_LESMs,1);
errTol_UESM = 1e-4;
rng(seed_final_state, 'twister');
formal_metadata = struct('initial_source','value_initial.mat', ...
    'b_range',[0.03,0.08],'c_range',[2e-3,6e-3], ...
    'outer_delay_factor',0.5,'outer_delay_cap',10, ...
    'outer_avg_start_iter',1,'min_outer_iter',100,'max_outer_iter',1500);
save(fullfile(root_dir, 'param.mat'), ...
    'a','b','c','D','start_idx','end_idx','init_pg','init_pes','init_pb','init_ps', ...
    'init_pesc','num_LESMs','num_prosumers','pi_max','pi_min','u','F_l', ...
    'alpha_l','alpha_PB','PTDF','init_rho','init_pi_PB','init_pi_l', ...
    'init_pi_l_pos','init_pi_l_neg','init_pi_0','init_pi','errTol_LESMs','errTol_UESM', ...
    'pg_max','pg_min','a_extend','beta_qre','rng_seeds','formal_metadata');
fprintf('param.mat generated from value_initial.mat with modified b/c parameters.\n');
end
