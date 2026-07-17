function res = localMarket_S1(varargin)
% 兼容入口：S1 内层异步保持、无遍历平均。
res = localMarketScenario(varargin{:}, true, false);
end
