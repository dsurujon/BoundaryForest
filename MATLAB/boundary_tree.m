%DS May 2017
%Boundary Forest to reduce input data size - Jose Bento
%[input] seqs: sequences (output from fastaread)
%[input] dist: distance method
%   sw: Smith-Waterman
%   en: relative entropy
%[input] q: word size for entropy calculation
%[input] eps: similarity threshold
%[input] max_deg: maximum number of children a node is allowed to have
%[output] tree: final tree with few redundancies
%[output] tree_node_ref: record of sequences and the tree nodes they were 
%assigned to, if the sequence was not included in the tree.
%[output] data_order_ix: the order data was read in

function [tree,tree_node_ref,data_order_ix] = boundary_tree(seqs, eps, max_deg,varargin)
    n=length(seqs);
    data_order_ix=randperm(n);
    tree={};
    p=inputParser;
    addParameter(p,'dist','sw');
    addParameter(p,'q',10);
    p.parse(varargin{:});
    %keep a record of which sequences are within the similarity threshold
    %of which tree node
    tree_node_ref=zeros(1,n);
    num_tree_nodes=0;
    % the first two nodes always go to the tree
    tree{1} = {data_order_ix(1),[2]};
    tree{2} = {data_order_ix(2),[]};
    num_tree_nodes = 2;
    tree_node_ref(data_order_ix(1))=data_order_ix(1);
    tree_node_ref(data_order_ix(2))=data_order_ix(2);
    for i = 3:n
        ele_being_proc_data_ix = data_order_ix(i);
        curr_tree_node_ix = 1;
        while(1)
            curr_node_data_ix = tree{curr_tree_node_ix}{1};
            curr_node_all_children_node_ix = tree{curr_tree_node_ix}{2};
            % this finds the child closest to the current point being processed
            smallest_dist = inf;
            best_child_node_ix = -1;
            for child_node_ix = curr_node_all_children_node_ix  
                
                if p.Results.dist=='sw'
                    dis = seqpdist([seqs(tree{child_node_ix}{1}),seqs(ele_being_proc_data_ix)],'Method','alignment-score','ScoringMatrix','BLOSUM62','GapOpen',10,'ExtendGap',0.5,'PairwiseAlignment',true);
                elseif p.Results.dist=='en'
                    dis=rel_entropy(seqs(tree{child_node_ix}{1}),seqs(ele_being_proc_data_ix),p.Results.q);
                end
                if (dis < smallest_dist)
                    smallest_dist = dis;
                    best_child_node_ix = child_node_ix;
                end        
            end
            % distance between the current data point and the current father
            if p.Results.dist=='sw'
                dis2 = seqpdist([seqs(curr_node_data_ix),seqs(ele_being_proc_data_ix)],'Method','alignment-score','ScoringMatrix','BLOSUM62','GapOpen',10,'ExtendGap',0.5,'PairwiseAlignment',true);
            elseif p.Results.dist=='en'
                dis2=rel_entropy(seqs(curr_node_data_ix),seqs(ele_being_proc_data_ix),p.Results.q);
            end

            % we add the new_point as a new child of the current father only if the father is
            % not saturated with children AND the father-new_point distance is
            % otherwise we move on to the next closest children
            if (dis2 < smallest_dist && length(tree{curr_tree_node_ix}{2}) < max_deg)
                % we only add the point if the point is not very close to the current father
                if (dis2 > eps)
                    num_tree_nodes = num_tree_nodes + 1;
                    tree{curr_tree_node_ix}{2} = [tree{curr_tree_node_ix}{2} , num_tree_nodes ];
                    tree{num_tree_nodes} = {ele_being_proc_data_ix, [] };
                    tree_node_ref(ele_being_proc_data_ix)=ele_being_proc_data_ix;
                % otherwise we put the current father as its reference
                else
                     tree_node_ref(ele_being_proc_data_ix)=tree{curr_tree_node_ix}{1};
                end
                break;
            else
                curr_tree_node_ix = best_child_node_ix;
            end
        end    
    end

end