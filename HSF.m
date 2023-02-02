function H = HSF(l, m, r, t, p)

%% Global.
global k

%% HSF.
H = 1i.^(-l).*SSH(l, m, t, p).*sphhnk(l, 2, k.*r);

end