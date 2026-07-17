function summary = run_outer_dual_direction_test()
%RUN_OUTER_DUAL_DIRECTION_TEST  Test the signs of the outer dual updates.
% This is a synthetic test and does not call any community solver.

alpha = 0.1;
ptdf = [1; -1];
F = 1;
tol = 1e-12;

% The controlled response is increasing in pi_0, as in pesc=(pi_0-lambda)/a.
chi = 1;

% PB: positive residual must lower pi_PB and lower the next response.
pb_pos = [1; 1];
pb_neg = -pb_pos;
[pb_pos_before,pb_pos_after,pb_pos_pi] = pbStep(pb_pos,0,alpha,chi);
[pb_neg_before,pb_neg_after,pb_neg_pi] = pbStep(pb_neg,0,alpha,chi);

% Positive line violation: y_pos is stored as a non-positive multiplier.
line_pos = [2; -2];
flow_pos_before = ptdf' * line_pos;
g_pos = flow_pos_before - F;
y_pos = min(0 - alpha*g_pos,0);
pi0_pos = ptdf*y_pos;
line_pos_after = line_pos + chi*pi0_pos;
flow_pos_after = ptdf' * line_pos_after;

% Negative line violation: y_neg is also stored as a non-positive multiplier.
line_neg = -line_pos;
flow_neg_before = ptdf' * line_neg;
g_neg = -flow_neg_before - F;
y_neg = min(0 - alpha*g_neg,0);
pi0_neg = ptdf*(-y_neg);
line_neg_after = line_neg + chi*pi0_neg;
flow_neg_after = ptdf' * line_neg_after;

summary = struct();
summary.pb_positive = struct('residual_before',pb_pos_before, ...
    'residual_after',pb_pos_after,'pi_after',pb_pos_pi);
summary.pb_negative = struct('residual_before',pb_neg_before, ...
    'residual_after',pb_neg_after,'pi_after',pb_neg_pi);
summary.line_positive = struct('flow_before',flow_pos_before, ...
    'flow_after',flow_pos_after,'violation_before',max(flow_pos_before-F,0), ...
    'violation_after',max(flow_pos_after-F,0),'dual_after',y_pos, ...
    'pi0_component',pi0_pos);
summary.line_negative = struct('flow_before',flow_neg_before, ...
    'flow_after',flow_neg_after,'violation_before',max(-flow_neg_before-F,0), ...
    'violation_after',max(-flow_neg_after-F,0),'dual_after',y_neg, ...
    'pi0_component',pi0_neg);
summary.pb_positive_pass = pb_pos_after < pb_pos_before - tol;
summary.pb_negative_pass = abs(pb_neg_after) < abs(pb_neg_before) - tol;
summary.line_positive_pass = summary.line_positive.violation_after < ...
    summary.line_positive.violation_before - tol;
summary.line_negative_pass = summary.line_negative.violation_after < ...
    summary.line_negative.violation_before - tol;
summary.success_gate = summary.pb_positive_pass && summary.pb_negative_pass && ...
    summary.line_positive_pass && summary.line_negative_pass;
save('outer_dual_direction_test.mat','summary');
fprintf('Outer dual direction test: %s\n',ternary(summary.success_gate,'PASS','FAIL'));

    function [before,after,pi_next] = pbStep(base,pi_current,alpha_local,chi_local)
        before = abs(sum(base + chi_local*pi_current));
        pi_next = pi_current - alpha_local*sum(base + chi_local*pi_current);
        after = abs(sum(base + chi_local*pi_next));
    end
    function out = ternary(condition,yes_value,no_value)
        if condition, out = yes_value; else, out = no_value; end
    end
end
