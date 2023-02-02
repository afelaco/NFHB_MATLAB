function R = RSF(l, m, r, t, p)

%% Global.
global k

%% RSF.

if l == 0
    
    R = 1i.^(-l).*(sphbsl(l, 1, k.*r) - 1i.*sphbsl(l+1, 1, k.*r)).*SSH(l, m, t, p);
    
else
    
    R = 1i.^(-l).*(sphbsl(l, 1, k.*r) + 1i.*sphbsl(l-1, 1, k.*r)).*SSH(l, m, t, p);
    
end