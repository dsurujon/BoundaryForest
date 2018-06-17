% 052918 DS
% Run each point through a BF to determine its closest representative
% Traverse the entire tree, always going to the closest child node
% Keep track of nodes traversed
% Assign node on the path that is closest to current item

function X = extendBF_v3(items,tree)
    n = length(items);
    X = zeros(n,1);
    % for each new item
    for i=1:n
        curr_tree_node_ix = 1;
        curr_node_all_children_node_ix = tree{curr_tree_node_ix}{2};
        
        curr_node_dist = seqpdist([items(tree{curr_tree_node_ix}{1}),...
            items(i)],'Method','alignment-score',...
            'ScoringMatrix','BLOSUM62','GapOpen',10,'ExtendGap',0.5,...
            'PairwiseAlignment',true);
        
        path_ixs = [curr_tree_node_ix];
        path_dists = [curr_node_dist];
        while ~isempty(curr_node_all_children_node_ix)
            smallest_dist = 1000000;
            for child_node_ix = curr_node_all_children_node_ix  
                dis = seqpdist([items(tree{child_node_ix}{1}),...
                    items(i)],'Method','alignment-score',...
                    'ScoringMatrix','BLOSUM62','GapOpen',10,'ExtendGap',0.5,...
                    'PairwiseAlignment',true);
                
                if (dis < smallest_dist)
                    smallest_dist = dis;
                    best_child_node_ix = child_node_ix;
                end    
            end
            path_ixs = [path_ixs best_child_node_ix];
            path_dists = [path_dists smallest_dist];
            curr_tree_node_ix = best_child_node_ix;
            curr_node_all_children_node_ix = tree{curr_tree_node_ix}{2};
        end
        
        [minpath, argminpath] = min(path_dists);
        X(i) = tree{path_ixs(argminpath)}{1};
    end
end