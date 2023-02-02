function output = sphhnk(l, k, x)

output = sqrt(pi./(2.*x)).*besselh(l+0.5, k, x);

end