% given a list of errors corresponding to clustering outputs at different 
% numbers of clusters, find the best number of clusters. This will be the
% point with the maximal curvature (looks like an elbow if you plot the
% errors)
% [input] errorlist: list of error values
% [input] n: number of clusters corresponding to the error values
% [output] k: correct number of clusters
% [output] i: index of correct number of clusters

function [k,i] = find_elbow(errorlist,n)
    p = polyfit(n,errorlist,4);
    pd = polyder(p);
    pdd = polyder(pd);
    [k,i] = max(polyval(pdd, n));
    k=n(i);
end