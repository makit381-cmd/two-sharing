function seed = deriveDeterministicSeed(base_seed,community_id,local_call_count,stream_tag)
%DERIVEDETERMINISTICSEED Derive a reproducible keyed MATLAB RandStream seed.

validateattributes(base_seed,{'numeric'},{'scalar','real','finite'});
validateattributes(community_id,{'numeric'},{'scalar','real','finite','nonnegative','integer'});
validateattributes(local_call_count,{'numeric'},{'scalar','real','finite','positive','integer'});
validateattributes(stream_tag,{'numeric'},{'scalar','real','finite','positive','integer'});

modulus = 2^31 - 1;
seed = mod(double(base_seed) ...
    + 104729 * double(community_id) ...
    + 13007 * double(local_call_count) ...
    + 1009 * double(stream_tag),modulus);
if seed <= 0
    seed = seed + 1;
end
seed = floor(seed);
end
