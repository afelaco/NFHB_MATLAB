function P = ALF(l, m, x)

%% Initialization.
P = zeros(size(x));

%% Closed form (https://en.wikipedia.org/wiki/Associated_Legendre_polynomials).
for k = abs(m) : l
    
    P = P + (-1).^abs(m).*2.^l.*(1-x.^2).^(abs(m)/2) ...
        ./factorial(l-k)./factorial(k-abs(m)).*x.^(k-abs(m)) ...
        .*gamma((l+k-1)/2+1)./gamma((l+k-1)/2-l+1) ...
        .*((1+sign(m)) + (1-sign(m)).*(-1).^abs(m).*factorial(l-abs(m))./factorial(l+abs(m)))./2;
    
end

end