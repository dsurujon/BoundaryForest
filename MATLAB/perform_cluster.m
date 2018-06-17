%DS May 2017
%Perform cluster analysis
%[input] tree: output of boundary tree
%[input] tref: output of boundary tree
%[input] dm: distance matrix
%[input] method: clustering method ('km': k-means, 'hi': hierarchical, 'sp':
%spectral)
%[input - optional] k: number of clusters for k-means or spectral
%[input - optional] vectorizeDM: whether or not to use the vectorized form
%of the input data. 
%   1: vectorize the data from the distance matrix (default)
%   0: use the distance matrix as input for k-means
%[input - optional] m: max number of clusters for hierarchical clustering
%[input - optional] linkage_method: method to be used for linkage
%calculation in hierarchical clustering
%[input - optional] norm_type: Spectral clustering normalization type. 
%	1: Unnormalized
%	2: Normalized according to Shi and Malik (2000)
%	3: Normalized according to Jordan and Weiss (2002)
%[input - optional] sigma: sigma to be used for Gaussian similarity matrix
%computation in spectral clustering
%[output] centroids: cluster centroids
%[output] clusters: cluster asignments on the tree nodes
%[output] all_clusters: cluster assignments on all input data

function [centroids, clusters, all_clusters] = perform_cluster(tree,tref,dm,method,varargin)
    centroids=[];
    clusters=[];
    all_clusters=zeros(length(tref),1);
    
    p=inputParser;
    addParameter(p,'k',30);
    addParameter(p,'m',30);
    addParameter(p,'linkage_method','ward');
    addParameter(p,'vectorizeDM',1);
    addParameter(p,'norm_type',1);
    addParameter(p,'sigma',0.1);
    p.parse(varargin{:});
    
    if method=='km'
        %%%%%%%find best k
        k=p.Results.k;
        if p.Results.vectorizeDM==1
            dmv=vectorize_dm_fn(dm);
        else
            dmv=dm;
        end
        [clusters,centroids]=kmeans(dmv,k);
    end

    if method=='hi'
        m=p.Results.m;
        linkages=linkage(dm,p.Results.linkage_method);
        clusters=cluster(linkages,'maxclust',m);
    end
    if method=='sp'
        W=simGaussian(dm,p.Results.sigma);
        [clusters, L,U] = SpectralClustering(W,p.Results.k,p.Results.norm_type);
    end
    
    for i=1:length(clusters)
        tree_ix=tree{i}{1};
        all_clusters(tref==tree_ix)=clusters(i);
    end
end