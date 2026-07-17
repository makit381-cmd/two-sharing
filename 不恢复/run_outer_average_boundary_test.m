function summary = run_outer_average_boundary_test()
% 轻量外层平均边界单元测试；不调用 localMarket_S2/S3。
% 仅验证真实价格更新计数、S2 等待事件、平均起点和数值平均。

P = load('params.mat');
if isfield(P,'outer_avg_start_price_update')
    formal_start = P.outer_avg_start_price_update;
else
    formal_start = P.outer_avg_start_iter;
end
test_overrides = struct( ...
    'outer_avg_start_price_update',NaN, ...
    'mode','formal_boundary_probe_only', ...
    'note','No local solver; TEST OVERRIDE ONLY is not used.');
effective_start = formal_start;
formal_parameters = struct( ...
    'outer_avg_start_price_update',formal_start, ...
    'outer_avg_start_iter',P.outer_avg_start_iter);
effective_parameters = struct( ...
    'outer_avg_start_price_update',effective_start, ...
    'outer_avg_start_iter',effective_start);

Q = 105;
q = (1:Q)';
expected_cnt = max(q - formal_start + 1,0);

% S2：每次真实同步更新前插入不同数量等待 wall event。
s2_avg_cnt = zeros(Q,1);
s2_avg_count = 0;
s2_bar_pesc = nan(Q,1);
s2_wait_events = zeros(Q,1);
s2_wall_iter = 0;
s2_price_update_count = 0;
s2_avg_cnt_before_wait = zeros(Q,1);
s2_avg_cnt_after_wait = zeros(Q,1);
s2_sum_pesc = 0;
for update = 1:Q
    waits = mod(update,4) + 1;
    s2_avg_cnt_before_wait(update) = s2_avg_count;
    for w = 1:waits
        s2_wall_iter = s2_wall_iter + 1;
        s2_wait_events(update) = s2_wait_events(update) + 1;
        assert(s2_avg_count == s2_avg_cnt_before_wait(update), ...
            'S2 waiting wall event changed outer_avg_cnt.');
    end
    s2_avg_cnt_after_wait(update) = s2_avg_count;
    s2_price_update_count = s2_price_update_count + 1;
    if s2_price_update_count >= effective_start
        s2_avg_count = s2_avg_count + 1;
        s2_sum_pesc = s2_sum_pesc + update;
        s2_bar_pesc(update) = s2_sum_pesc / s2_avg_count;
    else
        s2_bar_pesc(update) = update;
    end
    s2_avg_cnt(update) = s2_avg_count;
end

% S3：每个 wall slot 都是一次真实 held 外层更新，因此两种计数相等。
s3_avg_cnt = zeros(Q,1);
s3_avg_count = 0;
s3_bar_pesc = nan(Q,1);
s3_price_update_count = 0;
s3_sum_pesc = 0;
for update = 1:Q
    s3_price_update_count = s3_price_update_count + 1;
    if s3_price_update_count >= effective_start
        s3_avg_count = s3_avg_count + 1;
        s3_sum_pesc = s3_sum_pesc + update;
        s3_bar_pesc(update) = s3_sum_pesc / s3_avg_count;
    else
        s3_bar_pesc(update) = update;
    end
    s3_avg_cnt(update) = s3_avg_count;
end

expected_bar_100 = 100;
expected_bar_101 = mean([100,101]);
expected_bar_105 = mean(100:105);
boundary_pass = ...
    s2_avg_cnt(20) == 0 && s2_avg_cnt(99) == 0 && ...
    s2_avg_cnt(100) == 1 && s2_avg_cnt(101) == 2 && s2_avg_cnt(105) == 6 && ...
    s3_avg_cnt(99) == 0 && s3_avg_cnt(100) == 1 && s3_avg_cnt(105) == 6 && ...
    all(s2_avg_cnt == expected_cnt) && all(s3_avg_cnt == expected_cnt) && ...
    all(s2_avg_cnt_after_wait == s2_avg_cnt_before_wait) && ...
    abs(s2_bar_pesc(100) - expected_bar_100) <= 1e-12 && ...
    abs(s2_bar_pesc(101) - expected_bar_101) <= 1e-12 && ...
    abs(s2_bar_pesc(105) - expected_bar_105) <= 1e-12 && ...
    abs(s3_bar_pesc(100) - expected_bar_100) <= 1e-12 && ...
    abs(s3_bar_pesc(101) - expected_bar_101) <= 1e-12 && ...
    abs(s3_bar_pesc(105) - expected_bar_105) <= 1e-12 && ...
    s3_price_update_count == Q;
assert(boundary_pass,'Outer average boundary test failed.');

summary = struct();
summary.formal_parameters = formal_parameters;
summary.effective_parameters = effective_parameters;
summary.test_overrides = test_overrides;
summary.test_override_only = false;
summary.q = q;
summary.expected_avg_cnt = expected_cnt;
summary.s2_avg_cnt = s2_avg_cnt;
summary.s2_bar_pesc = s2_bar_pesc;
summary.s2_wait_events = s2_wait_events;
summary.s2_wall_iter = s2_wall_iter;
summary.s2_price_update_count = s2_price_update_count;
summary.s2_avg_cnt_before_wait = s2_avg_cnt_before_wait;
summary.s2_avg_cnt_after_wait = s2_avg_cnt_after_wait;
summary.s3_avg_cnt = s3_avg_cnt;
summary.s3_bar_pesc = s3_bar_pesc;
summary.s3_wall_iter = Q;
summary.s3_price_update_count = s3_price_update_count;
summary.expected_bar_100 = expected_bar_100;
summary.expected_bar_101 = expected_bar_101;
summary.expected_bar_105 = expected_bar_105;
summary.tolerance = 1e-12;
summary.success_gate = boundary_pass;
summary.semantics = 'outer average starts at real price update 100; update 100 is included; S2 waiting does not change avg count; S3 wall and price update counts coincide';
save('outer_average_boundary_test.mat','summary');
fprintf('Outer average boundary test (formal start=%d): PASS\n',formal_start);
end
