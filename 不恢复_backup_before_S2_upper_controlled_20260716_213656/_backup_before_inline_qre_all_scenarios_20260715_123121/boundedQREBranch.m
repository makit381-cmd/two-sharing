function [pg, pes, pb, ps, cert] = boundedQREBranch(c, b, a, D, pg_min, pg_max, lambda, pi_max, pi_min, beta_qre, ~, z_cap, ~)
% 未认证 QRE 快速响应：保留一次截断正态扰动和完整物理可行构造，
% 但不再为每个 EU 解计算 epsilon 目标差、回退或重抽样。
% 仅限速度/行为试验；不能作为 HDEM A2/A5 误差证书。

z_pg = max(min(randn(1), z_cap), -z_cap);
z_pes = max(min(randn(1), z_cap), -z_cap);
[pg, pes, pb, ps, boundary_hit] = response(z_pg, z_pes);
if nargout >= 5
    cert = struct('certificate_passed',false, 'perturbed_accepted',true, ...
        'used_fallback',false, 'backoff_steps',0, 'boundary_hit',boundary_hit, ...
        'objective_gap',NaN, 'epsilon',NaN, 'z_pg',z_pg, 'z_pes',z_pes, ...
        'mode','unchecked_qre');
end

    function [pg_v, pes_v, pb_v, ps_v, boundary_hit] = response(zpg, zpes)
        tol_F = 1e-8;
        z_middle = (sqrt(c) * zpg + sqrt(a) * zpes) / sqrt(c + a);

        pg_buyer_unclipped = (pi_max - b) / c + sqrt(1 / (beta_qre * c)) * zpg;
        pg_buyer = min(max(pg_buyer_unclipped, pg_min), pg_max);
        pes_buyer = (lambda - pi_max) / a + sqrt(1 / (beta_qre * a)) * zpes;
        F_buyer = D + pes_buyer - pg_buyer;

        pg_seller_unclipped = (pi_min - b) / c + sqrt(1 / (beta_qre * c)) * zpg;
        pg_seller = min(max(pg_seller_unclipped, pg_min), pg_max);
        pes_seller = (lambda - pi_min) / a + sqrt(1 / (beta_qre * a)) * zpes;
        F_seller = D + pes_seller - pg_seller;

        if F_buyer > tol_F
            pg_v = pg_buyer;
            pes_v = pes_buyer;
            pb_v = F_buyer;
            ps_v = 0;
        elseif F_seller < -tol_F
            pg_v = pg_seller;
            pes_v = pes_seller;
            pb_v = 0;
            ps_v = -F_seller;
        else
            pes_unclipped = (lambda - b - c * D) / (c + a) + ...
                sqrt(1 / (beta_qre * (c + a))) * z_middle;
            pes_v = min(max(pes_unclipped, pg_min - D), pg_max - D);
            pg_v = D + pes_v;
            pb_v = 0;
            ps_v = 0;
        end

        boundary_tol = 1e-10 * max(1, max(abs([pg_min, pg_max])));
        boundary_hit = abs(pg_v - pg_min) <= boundary_tol || ...
            abs(pg_v - pg_max) <= boundary_tol;
    end

end
