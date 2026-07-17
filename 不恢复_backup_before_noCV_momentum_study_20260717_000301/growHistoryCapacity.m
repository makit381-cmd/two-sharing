function nextCapacity = growHistoryCapacity(currentCapacity, neededIndex, maxCapacity)
%GROWHISTORYCAPACITY 按需几何扩容一维历史；支持 Inf 作为无显式上限。
% 不参与任何市场更新、随机数或停止判断，仅决定已发生历史的存储容量。

validateattributes(currentCapacity, {'numeric'}, {'scalar','real','finite','integer','positive'});
validateattributes(neededIndex, {'numeric'}, {'scalar','real','finite','integer','positive'});
validateattributes(maxCapacity, {'numeric'}, {'scalar','real','positive'});

if neededIndex <= currentCapacity
    nextCapacity = currentCapacity;
    return;
end
if isfinite(maxCapacity) && neededIndex > maxCapacity
    error('Requested history index exceeds its configured maximum iteration capacity.');
end

nextCapacity = max(neededIndex, 2 * currentCapacity);
if isfinite(maxCapacity)
    nextCapacity = min(maxCapacity, nextCapacity);
end
end
