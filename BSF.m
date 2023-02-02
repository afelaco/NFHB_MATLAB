function J = BSF(l, m, r, t, p)

%% Global.
global k

%% HSF.
J = 1i.^(-l).*SSH(l, m, t, p).*sphbsl(l, 1, k.*r);

end