function [e, m] = VSHA(F, L, quadrature)

%% Quadrature.
if quadrature == "gauss"
    
    [t, p, w] = gauss(L);

elseif quadrature == "lebedev"
    
    [t, p, w] = lebedev(L);
    
else
    
    return

end

%% VSH.
tag = strcat('lib\', quadrature, '\', sprintf('VSH_%d.mat', L));

if isfile(tag) == 0
    
    for l = 0 : L
        for m = -l : l
            
            [E{lm2i(l, m)}, M{lm2i(l, m)}] = VSH(l, m, t, p);
            
        end
    end
    
    save(tag, 'E', 'M')
    
else
    
    load(tag)
    
end

E = cat(4, E{:});
M = cat(4, M{:});

%% Interpolate data.
F_int(:,:,1) = interp2(F.p, F.t, F.F(:,:,1), p, t, 'spline');
F_int(:,:,2) = interp2(F.p, F.t, F.F(:,:,2), p, t, 'spline');

%% Polar-to-spherical.
F_sph(:,:,1) = exp(-1i.*p)./sqrt(2).*(1i.*F_int(:,:,2) - cos(t).*F_int(:,:,1));
F_sph(:,:,2) = exp(1i.*p)./sqrt(2).*(1i.*F_int(:,:,2) + cos(t).*F_int(:,:,1));
F_sph(:,:,3) = -sin(t).*F_int(:,:,1);

%% Projection.
e = squeeze(sum(w.*F_sph.*conj(E), [1 2 3]));
m = squeeze(sum(w.*F_sph.*conj(M), [1 2 3]));

end