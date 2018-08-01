% Consensus clustering
%
% given a set of clusterings (as a matrix with each object as a row, and
% each clustering as a column), find a consensus clustering. The cluster
% assignments in clusterres_ext are used as feature vectors and
% hierarchically clustered into a new set of clusters. 
%
% [input] clusterres_ext: n X m matrix with n data points, and m clustering
% assignments
% [input] bestcluster: 1xm list of number of clusters used for each column
% (clustering) in clusterres_ext
% [output] consclust: consensus clustering assignments

function [consclust, kcons] = consensus_clustering(clusterres_ext, bestcluster)
kcons = mode(bestcluster);
disp(kcons);
cons_dists = pdist(clusterres_ext,'hamming');
cons_dists = squareform(cons_dists);
linkages=linkage(cons_dists,'ward');
consclust=cluster(linkages,'maxclust',kcons);
end