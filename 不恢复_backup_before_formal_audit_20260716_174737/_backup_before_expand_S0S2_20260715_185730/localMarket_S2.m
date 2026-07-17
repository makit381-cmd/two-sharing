function res = localMarket_S2(varargin)
% 兼容入口：S2 内层同步等待、遍历平均。
res = localMarketScenario(varargin{:}, false, true);
end
