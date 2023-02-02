function Y = SSH(l, m, t, p)

%% Order check.
if l < 0 || abs(m) > l
   
    Y = zeros(size(t));
    
    return
    
end

%% Normalization.
N = sqrt((2*l+1)/(4*pi)*factorial(l-m)/factorial(l+m));

%% Associated Legendre functions.
P = ALF(l, m, cos(t));

%% Phase.
E = exp(1i.*m.*p);

%% Closed form.
Y = N.*P.*E;

end