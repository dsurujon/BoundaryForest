% given a list of errors corresponding to clustering outputs at different 
% numbers of clusters, find the best number of clusters. This will be the
% point with the maximal curvature (looks like an elbow if you plot the
% errors)
% [input] errorlist: list of error values
% [input] n: number of clusters corresponding to the error values
% [output] k: correct number of clusters
% [output] i: index of correct number of clusters

function [k,i] = find_elbow_spline(errorlist,n)

    errorsmooth=smooth(errorlist,3);
    
    [pp, ~] = csaps( n(1,:),errorsmooth ,0.9 );
    p_der2 = fnder(pp,2);
    [k,i] = max(p_der2.coefs(:,2));
    
    k=n(i);
end