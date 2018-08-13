% Writes the clustering (consensus) output for Spectral SM
% into a csv file with three columns: 
% 1. Cluster ID
% 2. Gene (Locus Tag)
% 3. Strain
%
% [input] fastafile: sequences used for clustering
% [input] clusterfile: cluster data file
% [input] outfile: name of the output csv file

function make_cluster_csv(fastafile, clusterfile, outfile)


seqs = fastaread(fastafile);
load(clusterfile,'consclust');

mymat = strings(length(seqs),3);
for i = 1:length(seqs)
    mymat(i,3) = extractBetween(seqs(i).Header,'|','.');
    mymat(i,1) = string(consclust{5}(i));
    try
        mymat(i,2) = extractBetween(seqs(i).Header,'[gene=',']');
    catch
        mymat(i,2) = 'X';
    end
end

mymatheaders = string({'Cluster','Gene','Strain'});
mymat  = [mymatheaders;mymat];

cell2csv(outfile,mymat);



end