% given a list of errors corresponding to clustering outputs at different 
% numbers of clusters, and a threshold for error, find the correct number  
% clusters defined as the first time error is below the threshold. 
% [input] errorlist: list of error values
% [input] n: number of clusters corresponding to the error values
% [input] e: error threshold
% [input] s: whether to smooth the errorlist 
% [output] k: correct number of clusters
% [output] i: index of correct number of clusters

function [k,i] = find_best_cluster(errorlist,n,e,s)
    if (s)
        errorlist=smooth(errorlist,3);
    end
    
    %normalize to a range of [0,1]
    errorlist=errorlist/max(errorlist);

    i=1;
    while(errorlist(i)>e && i<length(errorlist))
       i=i+1;
    end
    k=n(i);
end