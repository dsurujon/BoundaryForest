function between_min=betweenMin(distances,clusters)
    between_min=5000000000;
    
	uclusters=unique(clusters);
    uclusters(isnan(uclusters))=[];
    i=length(uclusters);
    for j=1:i
        cluster=uclusters(j);
        indx=find(clusters==cluster);
        indy=find(clusters~=cluster);
        
        bet_clusters=distances(indx,indy);
        between_min=min(between_min,min(min(bet_clusters)));

    end

end