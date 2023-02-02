function output = sphbsl(l, k, x)

if k == 1
    
    output = sqrt(pi./(2.*x)).*besselj(l+0.5, x);
    
elseif k == 2
    
    output = sqrt(pi./(2.*x)).*bessely(l+0.5, x);
    
end