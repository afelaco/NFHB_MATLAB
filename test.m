clc
clear all
close all

%%
L = 11;

%% Global.
l = 3;
m = 2;

%%
[Grid, new_order] = GridLebedev(2*L);

c = CG(l, m, 1, 0, l, m);
Y = squeeze(TSH(l, m, l, Grid.Theta, Grid.Phi));