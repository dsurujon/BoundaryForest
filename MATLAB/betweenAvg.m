function between_avg=betweenAvg(distances,clusters)
    between_sum=0;
    between_count=0;

	uclusters=unique(clusters);
    uclusters(isnan(uclusters))=[];
    i=length(uclusters);
    for j=1:i
        cluster=uclusters(j);
        indx=find(clusters==cluster);
        indy=find(clusters~=cluster);
        
        bet_clusters=distances(indx,indy);
        between_count=between_count+length(indy)^2;
        between_sum=between_sum+sum(sum(bet_clusters));

    end
    between_avg=between_sum/between_count;
end