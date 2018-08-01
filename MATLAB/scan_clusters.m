%DS May 2017
%Perform cluster analysis
%[input] tree: output of boundary tree
%[input] tref: output of boundary tree
%[input] dm: distance matrix
%[input] method: clustering method ('km': k-means, 'hi': hierarchical, 'sp':
%spectral)
%[input] krange: number of clusters (as a list)
%[input - optional] vectorizeDM: whether or not to use the vectorized form
%of the input data. 
%   1: vectorize the data from the distance matrix (default)
%   0: use the distance matrix as input for k-means
%[input - optional] linkage_method: method to be used for linkage
%calculation in hierarchical clustering
%[input - optional] norm_type: Spectral clustering normalization type. 
%	1: Unnormalized
%	2: Normalized according to Shi and Malik (2000)
%	3: Normalized according to Jordan and Weiss (2002)
%[input - optional] sigma: sigma to be used for Gaussian similarity matrix
%computation in spectral clustering
%[output] clusters_tree: cluster asignments on the tree nodes
%[output] krange: adjusted krange (if necessary)

function [clusters_tree, krange] = scan_clusters(tree,krange,dm,method,varargin)
    % avoid k's that are larger than the size of the tree
    krange = krange(krange<size(tree,2));

    lenk = length(krange);
    clusters_tree = cell(1,lenk);
    
    p=inputParser;
    addParameter(p,'m',30);
    addParameter(p,'linkage_method','ward');
    addParameter(p,'vectorizeDM',1);
    addParameter(p,'norm_type',1);
    addParameter(p,'sigma',0.1);
    p.parse(varargin{:});
    
    if method=='km'
        if p.Results.vectorizeDM==1
            dmv=vectorize_dm_fn(dm);
        else
            dmv=dm;
        end
        for kix=1:lenk
            k=krange(kix);
            [clusters,centroids]=kmeans(dmv,k);
            clusters_tree{kix} = clusters;
        end
    end

    if method=='hi'
        linkages=linkage(dm,p.Results.linkage_method);
        for kix=1:lenk
            k=krange(kix);
            clusters=cluster(linkages,'maxclust',k);
            clusters_tree{kix} = clusters;
        end
    end
    if method=='sp'
        %W=simGaussian(dm,p.Results.sigma)>0.9;
        W = exp(-dm.^2)>0.9;
        [L,U,V] = Spectralpt1(W, p.Results.norm_type);
        for kix=1:lenk
            k=krange(kix);
            clusters=Spectralpt2(U,V,k, p.Results.norm_type);
            clusters_tree{kix} = clusters;
        end
    end
    
end