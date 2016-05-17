filename = 'bter_input.mat'
load(filename)
addpath('~/Downloads/feastpack_v1.1')
[ccd,gcc] = ccperdeg(G);
nd = accumarray(nonzeros(sum(G,2)),1);
[E1,E2] = bter(nd,ccd);
G_bter = bter_edges2graph(E1,E2);
save('bter_output', 'G_bter')
exit
