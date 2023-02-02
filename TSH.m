function Y = TSH(j, m, l, t, p)

%% Initialization
Y = zeros([size(t), 3]);

%% Closed form.
Y(:,:,1) = CG(l, m-1, 1, 1, j, m).*SSH(l, m-1, t, p);
Y(:,:,2) = CG(l, m+1, 1, -1, j, m).*SSH(l, m+1, t, p);
Y(:,:,3) = CG(l, m, 1, 0, j, m).*SSH(l, m, t, p);

end