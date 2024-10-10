function result = longestIncreasingSubsequence(nums)
    n = length(nums);
    dp = ones(1, n);
    for i = 1:n
        for j = 1:i
            if nums(j) < nums(i)
                dp(i) = max(dp(i), dp(j) + 1);
            end
        end
    end
    maxLength = max(dp);
    result = maxLength;
end

