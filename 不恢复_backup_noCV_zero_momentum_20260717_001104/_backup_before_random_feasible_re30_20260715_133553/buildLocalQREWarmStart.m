function state = buildLocalQREWarmStart(a,c,b,D,pg_min,pg_max,start_idx,end_idx, ...
    pi_0,pi_max,pi_min,beta_qre,z_pg,z_pes,noise_scale,root_iterations)
% 固定 QRE 随机方向下的社区 KKT 一致初值。
% 对每个社区二分 lambda，使 raw Agg 响应等于带误差 EU 响应之和。
% 本函数只保证 EU 可行性和社区内部一致性；外层平衡/线路由调用方审计。

num_communities = numel(a);
num_prosumers = numel(c);
if nargin < 16 || isempty(root_iterations)
    root_iterations = 40;
end
if isscalar(pi_0)
    pi_0 = repmat(pi_0,num_communities,1);
else
    pi_0 = pi_0(:);
end

state.pg = zeros(num_prosumers,1);
state.pes = zeros(num_prosumers,1);
state.pb = zeros(num_prosumers,1);
state.ps = zeros(num_prosumers,1);
state.pesc = zeros(num_communities,1);
state.lambda = zeros(num_communities,1);
state.pi_0 = pi_0;
state.inner_residual = zeros(num_communities,1);
state.noise_scale = noise_scale;

for i = 1:num_communities
    idx = start_idx(i):end_idx(i);
    [lambda_i,pg_i,pes_i,pb_i,ps_i] = solveCommunityQRE( ...
        a(i),c(idx),b(idx),D(idx),pg_min(idx),pg_max(idx),pi_0(i), ...
        pi_max,pi_min,beta_qre(idx),z_pg(idx),z_pes(idx),noise_scale,root_iterations);
    state.pg(idx) = pg_i;
    state.pes(idx) = pes_i;
    state.pb(idx) = pb_i;
    state.ps(idx) = ps_i;
    state.lambda(i) = lambda_i;
    state.pesc(i) = (pi_0(i) - lambda_i) / a(i);
    state.inner_residual(i) = state.pesc(i) - sum(pes_i);
end
end

function [lambda,pg,pes,pb,ps] = solveCommunityQRE(a,c,b,D,pg_min,pg_max,pi_0, ...
    pi_max,pi_min,beta,z_pg,z_pes,noise_scale,root_iterations)
lambda_lo = -1;
lambda_hi = 1;
r_lo = qreResidual(lambda_lo);
r_hi = qreResidual(lambda_hi);
expansion = 0;
while ~(r_lo >= 0 && r_hi <= 0) && expansion < 80
    lambda_lo = 2 * lambda_lo - 1;
    lambda_hi = 2 * lambda_hi + 1;
    r_lo = qreResidual(lambda_lo);
    r_hi = qreResidual(lambda_hi);
    expansion = expansion + 1;
end
if ~(r_lo >= 0 && r_hi <= 0)
    error('buildLocalQREWarmStart:BracketFailed', ...
        'Failed to bracket a fixed-noise QRE community dual root.');
end
for q = 1:root_iterations
    lambda_mid = 0.5 * (lambda_lo + lambda_hi);
    if qreResidual(lambda_mid) >= 0
        lambda_lo = lambda_mid;
    else
        lambda_hi = lambda_mid;
    end
end
lambda = 0.5 * (lambda_lo + lambda_hi);
[pg,pes,pb,ps] = qreResponse(lambda,a,c,b,D,pg_min,pg_max, ...
    pi_max,pi_min,beta,z_pg,z_pes,noise_scale);

    function r = qreResidual(lambda_v)
        [~,pes_v] = qreResponse(lambda_v,a,c,b,D,pg_min,pg_max, ...
            pi_max,pi_min,beta,z_pg,z_pes,noise_scale);
        r = (pi_0 - lambda_v) / a - sum(pes_v);
    end
end

function [pg,pes,pb,ps] = qreResponse(lambda,a,c,b,D,pg_min,pg_max, ...
    pi_max,pi_min,beta,z_pg,z_pes,noise_scale)
z_pg = noise_scale .* z_pg;
z_pes = noise_scale .* z_pes;
z_middle = (sqrt(c) .* z_pg + sqrt(a) .* z_pes) ./ sqrt(c + a);
pg_buyer = min(max((pi_max - b) ./ c + sqrt(1 ./ (beta .* c)) .* z_pg,pg_min),pg_max);
pes_buyer = (lambda - pi_max) ./ a + sqrt(1 ./ (beta .* a)) .* z_pes;
F_buyer = D + pes_buyer - pg_buyer;
pg_seller = min(max((pi_min - b) ./ c + sqrt(1 ./ (beta .* c)) .* z_pg,pg_min),pg_max);
pes_seller = (lambda - pi_min) ./ a + sqrt(1 ./ (beta .* a)) .* z_pes;
F_seller = D + pes_seller - pg_seller;
buyer = F_buyer > 1e-10;
seller = ~buyer & F_seller < -1e-10;
middle = ~(buyer | seller);
pg = zeros(numel(c),1);
pes = zeros(numel(c),1);
pb = zeros(numel(c),1);
ps = zeros(numel(c),1);
pg(buyer) = pg_buyer(buyer);
pes(buyer) = pes_buyer(buyer);
pb(buyer) = F_buyer(buyer);
pg(seller) = pg_seller(seller);
pes(seller) = pes_seller(seller);
ps(seller) = -F_seller(seller);
pes_middle = (lambda - b(middle) - c(middle) .* D(middle)) ./ (c(middle) + a) + ...
    sqrt(1 ./ (beta(middle) .* (c(middle) + a))) .* z_middle(middle);
pes(middle) = min(max(pes_middle,pg_min(middle) - D(middle)),pg_max(middle) - D(middle));
pg(middle) = D(middle) + pes(middle);
end
