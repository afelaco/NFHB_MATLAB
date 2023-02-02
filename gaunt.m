function out = gaunt(j1, m1, j2, m2, j3)

% Gaunt coefficient (https://dlmf.nist.gov/34.3.E22).
out = (-1)^(m1+m2)*sqrt((2*j1+1)*(2*j2+1)*(2*j3+1)/(4*pi)).*...
    wigner3j(j1, 0, j2, 0, j3, 0)*wigner3j(j1, m1, j2, m2, j3, -(m1+m2));

end