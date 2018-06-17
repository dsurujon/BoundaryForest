function within_max=withinMax(distances,clusters)
    within_max=0;
    
	uclusters=unique(clusters);
    uclusters(isnan(uclusters))=[];
    i=length(uclusters);
    for j=1:i
        cluster=uclusters(j);
        indx=find(clusters==cluster);
        
        within_clusters=distances(indx,indx);
        within_max=max(within_max,max(max(within_clusters)));

    end

end