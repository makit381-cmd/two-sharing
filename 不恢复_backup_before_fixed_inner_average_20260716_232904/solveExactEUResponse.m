function sol = solveExactEUResponse( ...
    c_j,b_j,a_i,D_j,pg_min_j,pg_max_j,lambda,pi_max,pi_min)
%SOLVEEXACTEURESPONSE Exact response for one EU at a fixed lambda.
% The three feasible piecewise-quadratic candidates are enumerated and
% compared using the complete local Lagrangian objective.

validateattributes(a_i, {'numeric'}, {'scalar','real','finite','positive'});
validateattributes(c_j, {'numeric'}, {'scalar','real','finite','positive'});
validateattributes(pg_min_j, {'numeric'}, {'scalar','real','finite'});
validateattributes(pg_max_j, {'numeric'}, {'scalar','real','finite'});
if pg_min_j > pg_max_j
    error('solveExactEUResponse: invalid generator bounds.');
end

feas_tol = 1e-8;
obj_tie_tol = 1e-10;

pg_cand = nan(3,1);
pes_cand = nan(3,1);
pb_cand = nan(3,1);
ps_cand = nan(3,1);
obj_cand = inf(3,1);
feasible = false(3,1);

% Buyer candidate: pb >= 0, ps = 0.
pg_cand(1) = min(max((pi_max - b_j) / c_j, pg_min_j), pg_max_j);
pes_cand(1) = (lambda - pi_max) / a_i;
pb_cand(1) = D_j + pes_cand(1) - pg_cand(1);
ps_cand(1) = 0;
if pb_cand(1) >= -feas_tol
    pb_cand(1) = max(pb_cand(1), 0);
    feasible(1) = true;
end

% Seller candidate: ps >= 0, pb = 0.
pg_cand(2) = min(max((pi_min - b_j) / c_j, pg_min_j), pg_max_j);
pes_cand(2) = (lambda - pi_min) / a_i;
pb_cand(2) = 0;
ps_cand(2) = pg_cand(2) - D_j - pes_cand(2);
if ps_cand(2) >= -feas_tol
    ps_cand(2) = max(ps_cand(2), 0);
    feasible(2) = true;
end

% Middle candidate: pb = ps = 0 and pes is projected to the physical set.
pes_cand(3) = min(max((lambda - b_j - c_j * D_j) / (c_j + a_i), ...
    pg_min_j - D_j), pg_max_j - D_j);
pg_cand(3) = D_j + pes_cand(3);
pb_cand(3) = 0;
ps_cand(3) = 0;
feasible(3) = true;

for branch = 1:3
    if feasible(branch)
        obj_cand(branch) = euLocalLagrangian( ...
            pg_cand(branch), pes_cand(branch), pb_cand(branch), ...
            ps_cand(branch), c_j, b_j, a_i, lambda, pi_max, pi_min);
    end
end

best_obj = min(obj_cand);
if ~isfinite(best_obj)
    error('solveExactEUResponse:noFeasibleCandidate', ...
        'No feasible piecewise candidate was found.');
end

% Fixed branch order resolves numerical ties deterministically.
best_branch = find(obj_cand <= best_obj + obj_tie_tol, 1, 'first');

sol = struct();
sol.pg = pg_cand(best_branch);
sol.pes = pes_cand(best_branch);
sol.pb = pb_cand(best_branch);
sol.ps = ps_cand(best_branch);
sol.obj = obj_cand(best_branch);
sol.branch = best_branch;
sol.balance_residual = D_j + sol.pes + sol.ps - sol.pg - sol.pb;
sol.feasibility_tolerance = feas_tol;
sol.objective_tie_tolerance = obj_tie_tol;
sol.candidate_objectives = obj_cand;
sol.candidate_feasible = feasible;
end
