function obj_star = exactLowerObj(c,b,a,D,pg_max,pg_min,pi_0,pi_max,pi_min)

N = length(b);

pg  = sdpvar(N,1);
pb  = sdpvar(N,1);
ps  = sdpvar(N,1);
pes = sdpvar(N,1);
y   = sdpvar(1,1);

Objective = sum(0.5 .* c(:) .* pg.^2 ...
    + b(:) .* pg ...
    + pi_max .* pb ...
    - pi_min .* ps ...
    + 0.5 .* a .* pes.^2) ...
    + 0.5 * a * y^2 ...
    - pi_0 * y;

Cons = [];

Cons = [Cons;
    pg_min(:) <= pg;
    pg <= pg_max(:);
    pb >= 0;
    ps >= 0;
    D(:) + pes + ps == pg + pb;
    y == sum(pes)];

ops = sdpsettings('solver','gurobi','verbose',0);
sol = optimize(Cons,Objective,ops);

if sol.problem == 0
    obj_star = value(Objective);
else
    error('exactLowerObj failed: %s', sol.info);
end

end