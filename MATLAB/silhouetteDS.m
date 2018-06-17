%silhouette value from distance matrix and cluster assignments
%si=Si = (bi-ai)/ max(ai,bi)
%where ai is the average distance from the ith point to the 
%other points in the same cluster as i, 
%and bi is the minimum average distance from the ith point 
%to points in a different cluster, minimized over clusters.

function sil = silhouetteDS(distances,clusters)
    n=length(clusters);
    s=zeros(length(n));
    uniqueclusts=unique(clusters);
    uniqueclusts(isnan(uniqueclusts))=[];
    bi=50000000;
    for i=1:n
       ai=avg_dist_pt_to_clust(i,clusters(i),distances,clusters);
       otherclusts=setdiff(uniqueclusts,clusters(i));
       for kk=1:length(otherclusts)
           j=otherclusts(kk);
           bs=avg_dist_pt_to_clust(i,j,distances,clusters);
           if(bs<bi)
               bi=bs;
           end
       end
       bi=min(bs);
       s(i)=(bi-ai)/ max(ai,bi);
    end
    sil=sum(s)/length(s);
end

%calculate average distance from point i to points in 
%cluster j
function d = avg_dist_pt_to_clust(i,j,distances,clusters)
    thisclusterix=(clusters==j);
    d=sum(distances(i,thisclusterix))/sum(thisclusterix);
end