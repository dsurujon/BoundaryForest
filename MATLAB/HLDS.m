function HLI=HLDS(dists,clusters)   
%%%Hubert-Levin index for clustering 
% HLI=(S-Smin)/(Smax-Smin) where S=sum of distances within clusters (l occurences)
% Smin= sum of the l minimum distances in the complete distance matrix (lower triangular)
% Smax= sum of the l maximum distances in the complete distance matrix (lower triangular)
% distances need to be lower triangular, 
% sorted_dists need to be the sorted form of distances
dists=tril(dists);
sorted_dists=sort(squareform(dists));
uclusters=unique(clusters);
uclusters(isnan(uclusters))=[];
i=length(uclusters);
S=0;
l=0;
% for each cluster get the number of within cluster comparisons, 
% add to l, and also add the distance sum to S
for j=1:i
	cluster=uclusters(j);
	indx=find(clusters==cluster);
	n=length(indx);
	l=l+n*(n-1)/2;
	within_clust=dists(indx,indx);
	S=S+sum(sum(within_clust));
end

numdist=length(sorted_dists);
Smin=sum(sorted_dists(1:l));
Smax=sum(sorted_dists(numdist-l:numdist));

HLI=(S-Smin)/(Smax-Smin);

end