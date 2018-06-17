function DI=dunnsDS(distances,clusters)   
%%%Dunn's index for clustering compactness and separation measurement
% min(intercluster distance)/max(cluster diameter)
uclusters=unique(clusters);
uclusters(isnan(uclusters))=[];
i=length(uclusters);
icl_dist=5000000000;
clus_diam=0;
for jj=1:i
    cluster=uclusters(jj);
    indx=(clusters==cluster);
    indy=(clusters~=cluster);
	
    bet_clusters=distances(indx,indy);
	icl_dist=min(icl_dist,min(min(bet_clusters)));
	
	within_clusters=distances(indx,indx);
	clus_diam=max(clus_diam,max(max(within_clusters)));
	
end

DI=icl_dist/clus_diam;
end