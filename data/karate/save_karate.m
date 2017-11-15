
fid = fopen('karate.paj.txt');
A = textscan(fid, '%f%f');
fclose(fid);
ind1=A{1};
ind2=A{2};
num_events=length(A{1});
num_nodes=numel(unique([A{2};unique(A{1})]));  
G = sparse([ind1;ind2], [ind2;ind1], ones(2*num_events,1));
G=G|G';
meta.ind1=ind1;
meta.ind2=ind2;
save karate.mat G meta
