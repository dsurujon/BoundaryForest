%DS May 2017
%take a distance matrix, recover a vectorized form of the underlying data
%using classical euclidean embedding
%[input] distances: nXn distance matrix
%[output] X: vectorized form of data, nXm matrix
function X = vectorize_dm_fn(distances)
    %make sure it's the full distance matrix
    if istril(distances)
        distances=distances+distances';
    end

    n=length(distances);
    %retain eigenvalues that will capture 99% of variabilty
    evalthreshold=0.30;
    M=zeros(n);
    for i=1:n
        for j=1:n
            M(i,j)=0.5*((distances(1,j)^2)+(distances(i,1)^2)-(distances(i,j)^2));
        end
    end
    [V,D]=eig(M);
    evals=diag(D);
    cumevals=cumsum(evals)/sum(evals);
    
    X=real(V*(D.^0.5));
    %only get the "nonzero" columns
    X=X(:,cumevals>evalthreshold);
    %X=X(:,1:5);
end

