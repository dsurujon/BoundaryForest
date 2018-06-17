% [input] A: cluster assignments A
% [input] B: cluster assignments B
% [output] costs: cost matrix of all possible matchings between A and B
% cost(i,j) = 1-(|Ai %% Bj|/|Ai U Bj|)

function costs = makeCostMatrixNorm2(A,B)
	if isrow(A)==0
		A=A';
	end
	if isrow(B)==0
		B=B';
	end
	n = length(A);
	ma = length(unique(A));
	mb = length(unique(B));
	costs = zeros(ma,mb);
    jj=1;
	for i=unique(A)
        %intersection
        in=sum((repmat(A==i,mb,1))+(B==unique(B)')==2,2);
        %union
        un=sum((repmat(A==i,mb,1))+(B==unique(B)')>0,2);
        q=1-in./un;
        costs(jj,:)=q';
        jj=jj+1;
        %disp(i/m);
    end
end