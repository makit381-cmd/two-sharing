function obj = exactUpperObj(a,c,b,pi_max,pi_min,D,F_l,pai_il, ...
    start_idx,end_idx,pg_max,pg_min,a_extend)

N = length(c);
num_LESMs = length(a);

pg   = sdpvar(N,1);
pb   = sdpvar(N,1);
ps   = sdpvar(N,1);
pes  = sdpvar(N,1);
pesc = sdpvar(num_LESMs,1);

Objective = sum(0.5 .* c .* pg.^2 ...
    + b .* pg ...
    + pi_max .* pb ...
    - pi_min .* ps ...
    + 0.5 .* a_extend .* pes.^2) ...    
    + sum(0.5 .* a .* pesc.^2);

Cons = [];

Cons = [Cons;
    pg_min <= pg;
    pg <= pg_max;
    pb >= 0;
    ps >= 0;
    D + pes + ps == pg + pb;
    sum(pesc) == 0;
    -F_l <= pai_il' * pesc;
    pai_il' * pesc <= F_l];

for i = 1:num_LESMs
    Cons = [Cons;
        pesc(i) == sum(pes(start_idx(i):end_idx(i)))];
end

ops = sdpsettings('solver','gurobi','verbose',0);
sol = optimize(Cons,Objective,ops);

if sol.problem == 0
    obj = value(Objective);
else
    error('exactUpperObj failed: %s', sol.info);
end

end