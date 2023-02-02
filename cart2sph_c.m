function [r, t, p] = cart2sph_c(x, y, z)

r = sqrt(x.^2 + y.^2 + z.^2);
t = acos(z./r);

p = zeros(size(y));
p(y >= 0) = acos(x(y >= 0)./sqrt(x(y >= 0).^2+y(y >= 0).^2));
p(y < 0) = 2*pi - acos(x(y < 0)./sqrt(x(y < 0).^2+y(y < 0).^2));

t(isnan(t)) = 0;
p(isnan(p)) = 0;

end