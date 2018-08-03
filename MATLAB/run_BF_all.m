function run_BF_all(fastafile, krange, replicates, outdir)

%clear('all');

addpath('BF_clustering');

%fastafile = '.\fasta\T4first10.fasta';
eps = 0.1;
max_deg = 10;

%replicates = 10;
%krange = 3:15;

[filepath,name,ext] = fileparts(fastafile);
outfilename = strcat(name,'.mat');
outfilename = fullfile(outdir, outfilename);

parobj = parpool(replicates);

%make trees (n replicates)
seqs = fastaread(fastafile);
trees = cell(1,replicates);
tree_node_refs = cell(1,replicates);
data_order_ixs = cell(1,replicates);

disp('Making Boundary Forest');
parfor i=1:replicates
    disp(i);
    [tree,tree_node_ref,data_order_ix] = boundary_tree(seqs, eps, max_deg);
    trees{i} = tree;
    tree_node_refs{i} = tree_node_ref;
    data_order_ixs{i} = data_order_ix;
end

save(outfilename, 'trees','tree_node_refs','data_order_ixs','-v7.3');

%make dms (n replicates)
disp('Calculating distance matrices');
dms = cell(1,replicates);
parfor i=1:replicates
    disp(i);
    dms{i} = pairwise_distances(trees{i}, seqs);
end
save(outfilename, 'dms','-append');


% cluster
disp('Starting Clustering, scanning the range:');
disp(krange);

ALL_methods = {'Ward','Kmeans','Kmeans vectorized','Spectral NN',...
    'Spectral SM','Spectral JW'};
ALL_clusters = cell(6,length(krange), replicates);

for i=1:replicates
    [clusters_tree_HI,krangenew] = scan_clusters(trees{i},...
        krange,dms{i},'hi');
    [clusters_tree_KM,krangenew] = scan_clusters(trees{i},...
        krange,dms{i},'km','vectorizeDM',0);
    [clusters_tree_KMV,krangenew] = scan_clusters(trees{i},...
        krange,dms{i},'km','vectorizeDM',1);
    [clusters_tree_SP1,krangenew] = scan_clusters(trees{i},...
        krange,dms{i},'sp','norm_type',1,'sigma',0.1);
    [clusters_tree_SP2,krangenew] = scan_clusters(trees{i},...
        krange,dms{i},'sp','norm_type',2,'sigma',0.1);
    [clusters_tree_SP3,krangenew] = scan_clusters(trees{i},...
        krange,dms{i},'sp','norm_type',3,'sigma',0.1);
    for j=1:length(krangenew)
       ALL_clusters{1,j,i} = clusters_tree_HI{j};
       ALL_clusters{2,j,i} = clusters_tree_KM{j};
       ALL_clusters{3,j,i} = clusters_tree_KMV{j};
       ALL_clusters{4,j,i} = clusters_tree_SP1{j}; 
       ALL_clusters{5,j,i} = clusters_tree_SP2{j}; 
       ALL_clusters{6,j,i} = clusters_tree_SP3{j}; 
    end
end
save(outfilename, 'ALL_clusters','-append');


% pick best cluster - finding elbow
disp('Picking best clusters');
SSE = zeros(6, length(krange),replicates);
bestcluster_ix = zeros(6,replicates);
bestcluster = zeros(6,replicates);
bestSSE = zeros(6,replicates);

for method = 1:6
    for i=1:replicates
       for j=1:length(krange) 
           disp([i,j]);
           SSE(method,j,i) = SSEDS(dms{i},ALL_clusters{method,j,i});
       end
       [bestk, bestk_ix]=find_elbow(SSE(method,:,i),krange);
       bestcluster_ix(method,i) = bestk_ix;
       bestcluster(method,i) = bestk;
       bestSSE(method,i) = SSE(method,bestk_ix,i);
    end
end
save(outfilename, 'SSE','bestcluster_ix','bestcluster','bestSSE','-append');

%fprintf('Tree\tK\tSSE\n')
%disp([[1:replicates]' bestcluster' bestSSE']);

% plot SSE and best clusters determined
figfilename = strcat(name,'.png');
figfilename = fullfile(outdir,figfilename);
figure('visible','off');
for method = 1:6
    subplot(2,3,method);
    plot(krange, squeeze(SSE(method,:,:)));hold on;
    xlabel('number of clusters');ylabel('SSE');title(ALL_methods(method));
    scatter(krange(bestcluster_ix(method,:)),bestSSE(method,:),'filled','k');
end
saveas(gcf,figfilename);

disp('Extending best clusters');
% cluster results extended (n X replicates)
clusterres_ext = cell(1,6);
for method = 1:6
    clusres_thismethod = zeros(length(seqs),replicates);
    for i=1:replicates
       clusres_thismethod(:,i) = extendBF_DS(ALL_clusters{method,bestcluster_ix(i),i},...
            trees{i},tree_node_refs{i});
    end
    clusterres_ext{method} = clusres_thismethod;
end

save(outfilename, 'clusterres_ext','-append');

%consensus clustering
disp('Consensus clustering using hierarchical-ward');
consclust = cell(1,6);
kcons = zeros(1,6);

for method = 1:6
    [this_consclust,this_kcons] = consensus_clustering(clusterres_ext{method},bestcluster(method,:));
    consclust{method} = this_consclust;
    kcons(method) = this_kcons;
end

disp('Saving all data');
% save all data into a file

save(outfilename,'fastafile',...
    'krange',...
    'consclust',...
    'kcons','-append')

end
