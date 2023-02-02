function E = SSWS(c, r, t, p)

%% Global.
global k

%% Order.
L = sqrt(size(c, 1))-1;

%% SSWS.
c = reshape(c.', 1, 1, 1, size(c, 2), []);

for l = 0 : L
    for m = -l : l
       
        E{lm2i(l, m)} = -1i.*k.*c(1, 1, 1, :, lm2i(l, m)).*HSF(l, m, r, t, p);
        
    end 
end

E = sum(cat(5, E{:}), 5);

end