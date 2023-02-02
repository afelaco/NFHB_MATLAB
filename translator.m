function T = translator(r, t, p, L_t, L)

%% Global variables.
global k

%% Condition.
if r == 0
    
    T = eye((L_t+1)^2, (L+1)^2);
    
    return
    
end

%% Order and degree vectors.
i = (0 : (L_t+1)^2-1)';
la = floor(sqrt(i));
mu = i-la.*(la+1);

i = 0 : (L+1)^2-1;
l = floor(sqrt(i));
m = i-l.*(l+1);

q = reshape(0 : (L_t+L), 1, 1, []);

n = reshape(0 : (L_t+L), 1, 1, 1, []);

%% Library.
if isfile('lib\translator.mat') == 1
    
    load('lib\translator')
    
    for i = 1 : length(T_norm)
        
        T(:,:,i) = full(T_norm{i});
        
    end
    
    if size(T, 1) >= (L_t+1)^2 && size(T, 2) >= (L+1)^2
        
        T = T(1:(L_t+1)^2, 1:(L+1)^2, 1:(L_t+L+1));
        
    end
    
end
    
if isfile('lib\translator.mat') == 0 || size(T, 1) < (L_t+1)^2 || size(T, 2) < (L+1)^2
    
    clear T_norm

    % Normalization.
    B = abs(l-la) <= q & q <= l+la & abs(m-mu) <= q;                                                     % SSH.
    
    N = 4.*pi.*1i.^(-q).*(-1).^(la-l+mu).* ...
        (-1).^(m-mu).*sqrt((2.*l+1).*(2.*la+1).*(2.*q+1)./(4.*pi)).* ...
        sqrt((2.*q+1)/(4.*pi).*fct(q-m+mu)./fct(q+m-mu));
    
    N(~B) = 0;
    
    % Delta function (https://dlmf.nist.gov/34.2.E5).
    D = sqrt(fct(l+la-q).*fct(l-la+q).*fct(-l+la+q)./fct(l+la+q+1));
    
    % Wigner 3j symbols (https://dlmf.nist.gov/34.2.E4).
    B = n >= 0 & n >= la-q & n >= l-q & n <= l+la-q & n <= l & n <= la;
    
    W1 = (-1).^(l-la+n).*D.*fct(l).*fct(la).*fct(q)./ ...
        fct(n)./fct(q-la+n)./fct(q-l+n)./fct(l+la-q-n)./fct(l-n)./fct(la-n);
    
    W1(~B) = 0;
    
    W1 = sum(W1, 4);
    
    B = abs(m) <= l & abs(mu) <= la & abs(mu-m) <= q & ...
        n >= 0 & n >= la-q-m & n >= l-q-mu & n <= l+la-q & n <= l-m & n <= la-mu;
    
    W2 = (-1).^(l-la-mu+m+n).*D.* ...
        sqrt(fct(l+m).*fct(l-m).*fct(la-mu).*fct(la+mu).*fct(q+mu-m).*fct(q-mu+m))./ ...
        (fct(n).*fct(q-la+m+n).*fct(q-l+mu+n).*fct(l+la-q-n).*fct(l-m-n).*fct(la-mu-n));
    
    W2(~B) = 0;
    
    W2 = sum(W2, 4);
    
    % Precomputed translator.
    T = N.*W1.*W2;

    for i = 1 : size(T, 3)
        
        T_norm{i} = sparse(T(:,:,i));

    end
    
    save('lib\translator', 'T_norm')
    
end

%% Associated Legendre functions (https://en.wikipedia.org/wiki/Associated_Legendre_polynomials).
P = (-1).^abs(m-mu).*2.^(q)./fct(q-n)./fct(n-abs(m-mu)).*gamma((q+n-1)./2+1)./gamma((q+n-1)./2-q+1).* ...
    ((1+sign(m-mu)) + (1-sign(m-mu)).*(-1).^abs(m-mu).*fct(q-abs(m-mu))./fct(q+abs(m-mu)))./2.* ...
    cos(t).^(n-abs(m-mu)).*(1-cos(t).^2).^(abs(m-mu)/2);

P(n < abs(m-mu) | n > q) = 0;

P = sum(P, 4);

%% Phi phase term.
E = exp(1i.*(m-mu).*p);

%% Bessel function of first kind (https://dlmf.nist.gov/10.49.E7).
J = (-1i).^(n-q-1).*fct(q+n)./2.^(n)./fct(n)./fct(q-n).*exp(-1i.*k.*r)./(k.*r).^(n+1);

J(n > q) = 0;

J = real(sum(J, 4));

%% Translator.
T = sum(T.*P.*E.*J, 3);

end