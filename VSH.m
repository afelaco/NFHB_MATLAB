function [E, M] = VSH(l, m, t, p)

E = sqrt((l+1)/(2*l+1)).*TSH(l, m, l-1, t, p) + ...
        sqrt(l/(2*l+1)).*TSH(l, m, l+1, t, p);

M = TSH(l, m, l, t, p);

end