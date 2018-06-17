%DS May 2017
%From a set of sequences, get pairwise Smith-Waterman distances
%[input] tree: output of bounady_tree
%[input] seqs: sequences, output of fastaread
%[input] metric: distance metric
%   'sw': smith-waterman sequence distance
%   'en': relative entropy
%[input] q: word size for entropy
%[output] dm: pairwise distance matrix

function dm = pairwise_distances(tree,seqs,varargin)
    n=length(tree);
    
    tree_idx=zeros(1,n);
    for i=1:n
       tree_idx(i)=tree{i}{1}; 
    end
    select_seqs=seqs(tree_idx);
    p=inputParser;
    addParameter(p,'metric','sw');
    addParameter(p,'q',10);
    p.parse(varargin{:});
    met=p.Results.metric;
    if met=='sw'
        dm=seqpdist(select_seqs,...
                'Method','alignment-score','ScoringMatrix','BLOSUM62',...
                'GapOpen',10,'ExtendGap',0.5,'PairwiseAlignment',true);
        dm=squareform(dm);
    elseif met=='en'
        dm=rel_entropy_alt(select_seqs,p.Results.q);

    dm=dm+dm';
    
end