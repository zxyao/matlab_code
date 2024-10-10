function result = longestIncreasingSubsequence(nums)
    n = length(nums);
    dp = ones(1, n);
    prev = zeros(1, n);
    for i = 1:n
        for j = 1:i
            if nums(j) < nums(i) && dp(j) + 1 > dp(i)
                dp(i) = dp(j) + 1;
                prev(i) = j;
            end
        end
    end
    [~, endIdx] = max(dp);
    sequence = zeros(1, max(dp));
    idx = endIdx;
    for i = max(dp):-1:1
        sequence(i) = nums(idx);
        idx = prev(idx);
    end
    result = sequence;
end