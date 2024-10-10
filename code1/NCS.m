%no-relaxed channel sharing
function [alpha]=NCS(M,N)
%     alpha = randi([0, 1], 10, 5);
% 
%     for i = 1:M
%         while sum(alpha(i, :)) > 1
%             idx = find(alpha(i, :) == 1);
%             if ~isempty(idx)
%                 replace_idx = randi(length(idx));
%                 alpha(i, idx(replace_idx)) = 0;
%             else
%                 one_idx = find(alpha(i, :) == 0);
%                 replace_idx = randi(length(one_idx));
%                 alpha(i, one_idx(replace_idx)) = 1;
%             end
%         end
%     end
% 
%     for j = 1:N
%         while sum(alpha(:, j)) > 1
%             idx = find(alpha(:, j) == 1);
%             if ~isempty(idx)
%                 replace_idx = randi(length(idx));
%                 alpha(idx(replace_idx), j) = 0;
%             else
%                 one_idx = find(alpha(:, j) == 0);
%                 replace_idx = randi(length(one_idx));
%                 alpha(one_idx(replace_idx), j) = 1;
%             end
%         end
%     end
    alpha=zeros(M,N);
    i=randi([1,M]);
    j=randi([1,N]);
    alpha(i,j)=0.1;
end