function nextCapacity = growHistoryCapacity(currentCapacity, neededIndex, maxCapacity)
%GROWHISTORYCAPACITY Return a bounded, geometrically grown history capacity.
%   Histories start small and grow only when an actually reached iteration
%   needs more storage.  maxCapacity remains an iteration safety limit,
%   not an eager memory allocation request.

validateattributes(currentCapacity, {'numeric'}, {'scalar','integer','positive'});
validateattributes(neededIndex, {'numeric'}, {'scalar','integer','positive'});
validateattributes(maxCapacity, {'numeric'}, {'scalar','integer','positive'});

if neededIndex <= currentCapacity
    nextCapacity = currentCapacity;
    return;
end
if neededIndex > maxCapacity
    error('Requested history index exceeds its configured maximum iteration capacity.');
end

nextCapacity = min(maxCapacity, max(neededIndex, 2 * currentCapacity));
end
