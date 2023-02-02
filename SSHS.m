function F = SSHS(f, t, p, basis)

%% Order
L = sqrt(length(f))-1;

f = reshape(f.', 1, 1, 3, []);

%% Synthesis.
for l = 0 : L
    for m = -l : l
        
        Y = SSH(l, m, t, p);
        
        F_cart{lm2i(l, m)} = f(1, 1, :, lm2i(l, m)).*Y;
        
    end
end

F_cart = sum(cat(4, F_cart{:}), 4);

if basis == "cartesian"
    
    F(:,:,1) = F_cart(:,:,1);
    F(:,:,2) = F_cart(:,:,2);
    F(:,:,3) = F_cart(:,:,3);
    
end

%% Cartesian-to-polar.
if basis == "polar"
    
    F(:,:,1) = cos(t).*cos(p).*F_cart(:,:,1) + cos(t).*sin(p).*F_cart(:,:,2) - sin(t).*F_cart(:,:,3);
    F(:,:,2) = cos(p).*F_cart(:,:,2) - sin(p).*F_cart(:,:,1);
    
end

%% Cartesian-to-spherical.
if basis == "spherical"
    
    F(:,:,1) = 1./sqrt(2).*(1i.*F_cart(:,:,2) - F_cart(:,:,1));
    F(:,:,2) = 1./sqrt(2).*(F_cart(:,:,1) + 1i.*F_cart(:,:,2));
    F(:,:,3) = F_cart(:,:,3);
    
end

end