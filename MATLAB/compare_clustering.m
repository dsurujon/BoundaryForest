%DS May 2017
%compare two clustering assignments
%[input] dists: distance matrix for nodes in the boundary tree (size n')
%[input] realclusts_tree: 1Xn' cluster assignments 
%[input] newclusts_tree: 1Xn' second set of cluster assignments
%[input] realclusts_all: 1Xn cluster assignments 
%[input] newclusts_all: 1Xn second set of cluster assignments
%[output] H: nXn comparison matrix.
%[output] er: error rate

function [H,er]=compare_clustering(dists,realclusts_tree,newclusts_tree,realclusts_all,newclusts_all)
    n=length(realclusts_all);
    comp=zeros(n,n,2);
    for i=1:n
        for j=1:i-1
            comp(i,j,1)=(realclusts_all(i)==realclusts_all(j));
            comp(i,j,2)=(newclusts_all(i)==newclusts_all(j));
        end
    end
    mycomp=comp(:,:,1)+comp(:,:,2)';
    %lower half of comp has the realclusts assignments, and the upper half
    %has newclusts assignments. 
    H=mycomp;
    er=2*sum(sum(comp(:,:,1)~=comp(:,:,2)))/(n*n);
    
end
