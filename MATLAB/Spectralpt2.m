function [C] = Spectralpt2(U, k, Type)
U=U(:,1:k);
U=abs(U);
% in case of the Jordan-Weiss algorithm, we need to normalize
% the eigenvectors row-wise
if Type == 3
    U = bsxfun(@rdivide, U, sqrt(sum(U.^2, 2)));
end

% now use the k-means algorithm to cluster U row-wise
% C will be a n-by-1 matrix containing the cluster number for
% each data point
% C = kmeans(U, k, 'start', 'cluster', ...
%                  'EmptyAction', 'singleton');
C=kmeans(U,k);
% now convert C to a n-by-k matrix containing the k indicator
% vectors as columns
%C = sparse(1:size(D, 1), C, 1);

end