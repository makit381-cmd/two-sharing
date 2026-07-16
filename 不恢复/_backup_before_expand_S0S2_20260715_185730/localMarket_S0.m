function res = localMarket_S0(varargin)
% 兼容入口：S0 内层同步等待、无遍历平均。
res = localMarketScenario(varargin{:}, false, false);
end
