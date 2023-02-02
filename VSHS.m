function F = VSHS(e, m, t, p, basis)

%% Order
L = sqrt(length(e))-1;

%% Synthesis.
for l = 0 : L
    for m_l = -l : l
        
        [E, M] = VSH(l, m_l, t, p);

        F_sph{lm2i(l, m_l)} = e(lm2i(l, m_l)).*E + m(lm2i(l, m_l)).*M;
        
    end
end

F_sph = sum(cat(4, F_sph{:}), 4);

if basis == "spherical"
    
    F(:,:,1) = F_sph(:,:,1);
    F(:,:,2) = F_sph(:,:,2);
    
end

%% Spherical-to-polar.
if basis == "polar"
    
    F(:,:,1) = 1./sqrt(2).*(cos(t).*(exp(-1i.*p).*F_sph(:,:,2) - exp(1i.*p).*F_sph(:,:,1)) - sin(t).*F_sph(:,:,3));
    F(:,:,2) = -1i.*(exp(1i.*p).*F_sph(:,:,1) + exp(-1i.*p).*F_sph(:,:,2));

end

%% Spherical-to-cartesian.
if basis == "cartesian"
    
    F(:,:,1) = 1./sqrt(2).*(F_sph(:,:,2) - F_sph(:,:,1));
    F(:,:,2) = -1i./sqrt(2).*(F_sph(:,:,2) + F_sph(:,:,1));
    F(:,:,3) = F_sph(:,:,3);

end

end