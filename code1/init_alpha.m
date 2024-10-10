function [alpha] = init_alpha(M,N)
    alpha = zeros(M,N); % initalization matrix A
    for i=1:M
        for j=1:N
            alpha(i,j)=0;
        end
    end
    for i = 1:M
        for j = 1:N
            % generate random numbers
            r = rand();
            % if the current row-column sum is greater than 1, set the
            % element to 0
            if sum(alpha(i, 1:j-1)) + r + sum(alpha(1:i-1, j)) >= 1
                alpha(i, j) = 0.01;
            % otherwise set the current element to a random number
            else
                alpha(i, j) = r;
            end
        end
    end
end

