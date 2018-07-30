% [input] A: cluster assignments A with ma unique assignments
% [input] B: cluster assignments B with mb unique assignments
% [output] costs: cost matrix (size (ma+mb)X(ma+mb)) of all 
% possible matchings between A and B including virtual matchings. 
% cost(i,j) = |Ai || Bj| -|Ai %% Bj| if 1<=i<=ma and 1<=j<=mb (real-to-real match)
% cost = |Real Cluster| for real-to-virtual matches
% cost = 0 for virtual-to-virtual matches

function [costs,ma,mb] = makeCostMatrix(A,B)
	if isrow(A)==0
		A=A';
	end
	if isrow(B)==0
		B=B';
	end
	n = length(A);
	ma = length(unique(A)) - sum(isnan(A));
	mb = length(unique(B)) - sum(isnan(B));
	costs = zeros(ma+mb,ma+mb);
    jj=1;
	for i=unique(A)
        %intersection
        in=sum((repmat(A==i,mb,1))+(B==unique(B)')==2,2);
        %union
        un=sum((repmat(A==i,mb,1))+(B==unique(B)')>0,2);
        q=un-in;
        costs(jj,1:mb)=q';
        jj=jj+1;
        %disp(i/m);
    end
    disp([ma,mb,size(costs(1:ma,mb+1:ma+mb)),size(repmat(histcounts(A)',1,ma))]);
    costs(ma+1:ma+mb,1:mb)=repmat(histcounts(B,mb),mb,1);
    costs(1:ma,mb+1:ma+mb)=repmat(histcounts(A,ma)',1,ma);
end