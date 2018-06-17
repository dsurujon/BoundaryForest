% Compare two clustering assignments, possibly of different numbers 
% of clusters.
% [input] A: First clustering assignment
% [input] B: Secon clustering assignment
% [output] cost: cost value

function cost=computeCostNew(A,B)
    %make the cost matrix
    [cmx,ma,mb]=makeCostMatrix(A,B);
    %assignment of clusters
    [rowsol,cost,v,u,costMat]=lapjv(cmx);

end