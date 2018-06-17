function sse=SSEDS(distances,clusters)   
%%%SSE: sum of squared error for clustering quality
% SSE=sum(over k) (sum(over i) (x_i-c_k)^2)
% with (c_k) is the medioid of cluster k

uclusters=unique(clusters);
uclusters(isnan(uclusters))=[];
i=length(uclusters);
medioids=zeros(1,i);
sse=0;
% find the indices of medioids of each cluster, and the diameter 
% defined as the avg distance to the medioid
for j=1:i
	cluster=uclusters(j);
	indx=find(clusters==cluster);
	within_clusters=distances(indx,indx);
	% sum of distances to the medioid, and the index of the medioid in the indx list
	[cluster_sum,medioid_ix]=min(sum(within_clusters));
	d=within_clusters(medioid_ix,:);
	sse=sse+sum(d.^2);
end


end