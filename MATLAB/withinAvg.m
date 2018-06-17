function within_avg=withinAvg(distances,clusters)
    within_sum=0;
    within_count=0;
    
	uclusters=unique(clusters);
    uclusters(isnan(uclusters))=[];
    i=length(uclusters);
    for j=1:i
        cluster=uclusters(j);
        indx=find(clusters==cluster);

        within_clusters=distances(indx,indx);
        within_count=within_count+length(indx)^2;
        within_sum=within_sum+sum(sum(within_clusters));
    end
    within_avg=within_sum/within_count;
end