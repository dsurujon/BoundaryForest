function DBI=DBDS(distances,clusters)   
%%%Davies-Bouldin index for clustering compactness and scatter
% DBI=(1/k)Sum_{k}(max_{k~=k'} (S(c_k)+S(c_k'))/S(c_k,c_k'))
% with S(c_k) is the average distance to medioid c_k for 
% cluster k and S(c_k,c_k') is the distance between medioids 
% c_k and c_k' for clusters k and k' respectively

uclusters=unique(clusters);
uclusters(isnan(uclusters))=[];
i=length(uclusters);
medioids=zeros(1,i);
cluster_diams=zeros(1,i);
% find the indices of medioids of each cluster, and the diameter 
% defined as the avg distance to the medioid
for j=1:i
	cluster=uclusters(j);
    indx=find(clusters==cluster);
	within_clusters=distances(indx,indx);
	% sum of distances to the medioid, and the index of the medioid in the indx list
	[cluster_sum,medioid_ix]=min(sum(within_clusters));
	%diameter of the cluster = S(c_k)
	cluster_diams(j)=cluster_sum/length(indx);
	%disp([cluster,length(indx),medioid_ix]);
	medioids(j)=indx(medioid_ix);
end

numerator=(cluster_diams+cluster_diams');
denominator=distances(medioids,medioids);
numerator=numerator-diag(diag(numerator));

S=numerator./denominator;
%get rid of the diagonal
S=S-diag(diag(S));
%disp(sum(sum(isinf(S))));
%get rid of inf's
S(isinf(S))=0;
DBI=sum(max(S))/i;

end