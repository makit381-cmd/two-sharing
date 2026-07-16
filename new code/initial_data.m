function initial_data()
% Cold-start original parameters + S3 k=200 preheat in one callable file.
root_dir = fileparts(mfilename('fullpath'));
addpath(root_dir);

seed_cost_load = 112;
seed_community_cost = 123;
seed_pg_qre = 12;
seed_final_state = 115;

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
F_l = ones(num_LESMs - 1, 1) .* 999e3;
congested_lines_idx = [find_line_index(10,63); find_line_index(35,36); ...
    find_line_index(4,74); find_line_index(8,9); find_line_index(15,16); ...
    find_line_index(78,79); find_line_index(95,96)];
F_l(congested_lines_idx) = [0.5,8,1.5,21,9,10,15]' .* 1e3;

end_idx = cumsum(u);
start_idx = [1, end_idx(1:end-1) + 1];
c = zeros(num_prosumers,1);
b = zeros(num_prosumers,1);
D = zeros(num_prosumers,1);
a = zeros(num_LESMs,1);

% Original cold-start cost ranges.
rng(seed_cost_load, 'twister');
for j = 1:num_prosumers
    c(j) = 0.5e-3 + (1e-3 - 0.5e-3) * rand();
    b(j) = 0.01 + (0.05 - 0.01) * rand();
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
    i = find(end_idx >= j, 1, 'first');
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
    sigma_qre = (0.05 + (0.30 - 0.05) * rand()) * max(pg_max(j), 5);
    beta_qre(j) = 1 / (a_extend(j) * sigma_qre^2);
end

init_rho = 0.4 * ones(num_LESMs,1);
init_pes = zeros(num_prosumers,1);
init_pesc = zeros(num_LESMs,1);
init_pi_PB = 0.1;
init_pi_l_pos = zeros(num_LESMs - 1,1);
init_pi_l_neg = zeros(num_LESMs - 1,1);
init_pi_0 = init_pi_PB + PTDF * (init_pi_l_pos - init_pi_l_neg);
init_pi = init_pi_0 - a .* init_pesc;
init_pg = min(max(D + init_pes, pg_min), pg_max);
init_pb = max(D + init_pes - init_pg, 0);
init_ps = max(init_pg - D - init_pes, 0);

rng(seed_final_state, 'twister');
[initial_pg, initial_pes, initial_pb, initial_ps, initial_pesc, initial_pi, ...
    initial_pi_0, initial_pi_PB, initial_pi_l_pos, initial_pi_l_neg] = ...
    preheat_s3_200(c,b,a,D,init_rho,pg_max,pg_min,beta_qre, ...
    init_pg,init_pes,init_pb,init_ps,init_pesc,init_pi,init_pi_0, ...
    init_pi_PB,init_pi_l_pos,init_pi_l_neg,start_idx,end_idx,PTDF, ...
    alpha_PB,alpha_l,F_l,pi_max,pi_min);

assert(norm(initial_pi_0 - (initial_pi_PB + PTDF * ...
    (initial_pi_l_pos - initial_pi_l_neg)), inf) <= 1e-10, ...
    'value_initial 的价格状态不一致。');
preheat_metadata = struct('outer_iter',200,'outer_delay_factor',0.5, ...
    'outer_delay_cap',10,'b_range',[0.01,0.05],'c_range',[0.5e-3,1e-3]);
save(fullfile(root_dir, 'value_initial.mat'), ...
    'initial_pg','initial_pes','initial_pb','initial_ps','initial_pesc', ...
    'initial_pi','initial_pi_0','initial_pi_PB','initial_pi_l_pos', ...
    'initial_pi_l_neg','preheat_metadata');
fprintf('value_initial.mat generated from original cold-start S3 at k = 200.\n');
end

function [pg,pes,pb,ps,pesc,pi,pi_0,pi_PB,pi_l_pos,pi_l_neg] = ...
    preheat_s3_200(c,b,a,D,init_rho,pg_max,pg_min,beta_qre, ...
    pg,pes,pb,ps,pesc,pi,pi_0,pi_PB,pi_l_pos,pi_l_neg, ...
    start_idx,end_idx,PTDF,alpha_PB,alpha_l,F_l,pi_max,pi_min)

num_LESMs = numel(a);
n_comm = end_idx - start_idx + 1;
k0 = min(10, ceil(0.5 * sqrt(n_comm)));
delay = zeros(num_LESMs,1);
count = zeros(num_LESMs,1);
rng(999, 'twister');
for i = 1:num_LESMs
    delay(i) = randi([0,k0(i)]);
end

for k = 1:200
    for i = 1:num_LESMs
        idx = start_idx(i):end_idx(i);
        if count(i) < delay(i)
            count(i) = count(i) + 1;
        else
            res = localMarketS3(c(idx),b(idx),a(i),D(idx),init_rho(i), ...
                pg_max(idx),pg_min(idx),pi(i),pes(idx),pg(idx),pb(idx),ps(idx), ...
                pesc(i),5e-6,pi_0(i),pi_max,pi_min,false,beta_qre(idx));
            pi(i) = res.pi;
            pes(idx) = res.pes;
            pg(idx) = res.pg;
            pb(idx) = res.pb;
            ps(idx) = res.ps;
            pesc(i) = res.pesc;
            count(i) = 0;
            delay(i) = randi([0,k0(i)]);
        end
    end

    pi_PB = pi_PB - alpha_PB * sum(pesc);
    pi_l_pos = min(pi_l_pos - alpha_l * (PTDF' * pesc - F_l), 0);
    pi_l_neg = min(pi_l_neg - alpha_l * (-PTDF' * pesc - F_l), 0);
    pi_0 = pi_PB + PTDF * (pi_l_pos - pi_l_neg);
end
end
