function sol = solveFastQREResponse( ...
    c_j,b_j,a_i,D_j,pg_min_j,pg_max_j,lambda,pi_max,pi_min,beta_j, ...
    qre_z_cap,qre_noise_enabled,qre_stream)
%SOLVEFASTQRERESPONSE Fast projected QRE response without hard certification.
% The random direction is sampled once; this function is used by the formal
% large-scale profile.  Exact objective auditing is performed separately.

if nargin < 13 || isempty(qre_stream)
    qre_stream = RandStream('mt19937ar','Seed',31002);
end

if qre_noise_enabled
    z_pg_raw = max(min(randn(qre_stream,1),qre_z_cap),-qre_z_cap);
    z_pes_raw = max(min(randn(qre_stream,1),qre_z_cap),-qre_z_cap);
else
    z_pg_raw = 0;
    z_pes_raw = 0;
end

gamma = 1;
pg_buyer = min(max((pi_max-b_j)/c_j + ...
    gamma*sqrt(1/(beta_j*c_j))*z_pg_raw,pg_min_j),pg_max_j);
pes_buyer = (lambda-pi_max)/a_i + ...
    gamma*sqrt(1/(beta_j*a_i))*z_pes_raw;
F_buyer = D_j + pes_buyer - pg_buyer;

pg_seller = min(max((pi_min-b_j)/c_j + ...
    gamma*sqrt(1/(beta_j*c_j))*z_pg_raw,pg_min_j),pg_max_j);
pes_seller = (lambda-pi_min)/a_i + ...
    gamma*sqrt(1/(beta_j*a_i))*z_pes_raw;
F_seller = D_j + pes_seller - pg_seller;

if F_buyer > 1e-10
    sol.pg = pg_buyer;
    sol.pes = pes_buyer;
    sol.pb = F_buyer;
    sol.ps = 0;
    sol.branch = 1;
elseif F_seller < -1e-10
    sol.pg = pg_seller;
    sol.pes = pes_seller;
    sol.pb = 0;
    sol.ps = -F_seller;
    sol.branch = 2;
else
    z_middle = (sqrt(c_j)*z_pg_raw + sqrt(a_i)*z_pes_raw) / ...
        sqrt(c_j+a_i);
    pes_middle = (lambda-b_j-c_j*D_j)/(c_j+a_i) + ...
        gamma*sqrt(1/(beta_j*(c_j+a_i)))*z_middle;
    sol.pes = min(max(pes_middle,pg_min_j-D_j),pg_max_j-D_j);
    sol.pg = D_j + sol.pes;
    sol.pb = 0;
    sol.ps = 0;
    sol.branch = 3;
end

sol.obj = euLocalLagrangian(sol.pg,sol.pes,sol.pb,sol.ps, ...
    c_j,b_j,a_i,lambda,pi_max,pi_min);
sol.exact_obj = NaN;
sol.epsilon_gap = NaN;
sol.gamma = gamma;
sol.backoff_count = 0;
sol.fallback_to_exact = false;
sol.z_pg_raw = z_pg_raw;
sol.z_pes_raw = z_pes_raw;
sol.balance_residual = D_j + sol.pes + sol.ps - sol.pg - sol.pb;
end
