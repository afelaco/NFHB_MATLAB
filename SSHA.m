function f = SSHA(F, L, quadrature, basis)

%% Quadrature.
if quadrature == "gauss"
    
    [t, p, w] = gauss(L);
    
elseif quadrature == "lebedev"
    
    [t, p, w] = lebedev(2*L);
    
else
    
    return
    
end

%% SSH.
tag = strcat('lib\', quadrature, '\', sprintf('SSH_%d.mat', L));

if isfile(tag) == 0
    
    for l = 0 : L
        for m = -l : l
            
            Y{lm2i(l, m)} = SSH(l, m, t, p);
            
        end
    end
    
    save(tag, 'Y')
    
else
    
    load(tag)
    
end

Y = cat(4, Y{:});

%% Interpolate data.
F_int(:,:,1) = interp2(F.p, F.t, F.F(:,:,1), p, t, 'spline');
F_int(:,:,2) = interp2(F.p, F.t, F.F(:,:,2), p, t, 'spline');

%% Basis transformation.
if basis == "cartesian"
    
    % Polar-to-cartesian.
    F_cart(:,:,1) = cos(t).*cos(p).*F_int(:,:,1) - sin(p).*F_int(:,:,2);
    F_cart(:,:,2) = cos(t).*sin(p).*F_int(:,:,1) + cos(p).*F_int(:,:,2);
    F_cart(:,:,3) = -sin(t).*F_int(:,:,1);
    
    % Projection.
    f = permute(squeeze(sum(w.*F_cart.*conj(Y), [1 2])), [2 1]);
    
elseif basis == "spherical"
    
    % Polar-to-spherical.
    F_sph(:,:,1) = exp(-1i.*p)./sqrt(2).*(1i.*F_int(:,:,2) - cos(t).*F_int(:,:,1));
    F_sph(:,:,2) = exp(1i.*p)./sqrt(2).*(1i.*F_int(:,:,2) + cos(t).*F_int(:,:,1));
    F_sph(:,:,3) = -sin(t).*F_int(:,:,1);
    
    % Projection.
    f = permute(squeeze(sum(w.*F_sph.*conj(Y), [1 2])), [2 1]);
    
else
    
    return
    
end

end