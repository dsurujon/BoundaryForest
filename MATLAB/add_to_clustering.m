function [new_consclust, new_clusterres_ext] = add_to_clustering(newseqsfile, treeseqsfile, clusterdatafile)

newseqs = fastaread(newseqsfile);
treeseqs = fastaread(treeseqsfile);
load(clusterdatafile, 'trees', 'clusterres_ext', 'consclust');

new_consclust = cell(1,6);
% assign to a representative for each method, each BT
new_clusterres_ext = cell(1,6); 
for method_ix = 1:6
    ext_thismethod = zeros(length(newseqs), ntree);
   for tree_ix = 1:ntree 
        [rep_seq, rep_seq_tree_ix] = find_closest_on_BF(newseqs,trees{tree_ix},treeseqs);
        clusters_thistree = clusterres_ext{method_ix}(rep_seq, tree_ix);
        ext_thismethod(:, tree_ix)=clusters_thistree;
        
   end
   new_clusterres_ext{method_ix} = ext_thismethod;
   % find closest neighbor on consensus. 
    Idx = knnsearch(clusterres_ext{method_ix}, new_clusterres_ext{method_ix});
    new_consclust{method_ix} = consclust{method_ix}(Idx);
end


end