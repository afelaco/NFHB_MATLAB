function [t, p, w] = gauss(L)

%% Quadrature degree.
N = floor((L+1)/2);

N = N-1;
N1 = N+1;
N2 = N+2;
xu = linspace(-1,1,N1)';

%% Initial guess
y = cos((2*(0:N)'+1)*pi/(2*N+2))+(0.27/N1)*sin(pi*xu*N/N2);

%% Legendre-Gauss Vandermonde Matrix
L = zeros(N1,N2);

%% Derivative of LGVM
Lp = zeros(N1,N2);

%% Compute the zeros of the N+1 Legendre Polynomial using the recursion relation and the Newton-Raphson method
y0 = 2;

%% Iterate until new points are uniformly within epsilon of old points
while max(abs(y-y0)) > eps
    
    
    L(:,1) = 1;
    Lp(:,1) = 0;
    
    L(:,2) = y;
    Lp(:,2) = 1;
    
    for k=2:N1
        
        L(:,k+1) = ((2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1))/k;
        
    end
    
    Lp = (N2)*( L(:,N1)-y.*L(:,N2) )./(1-y.^2);
    
    y0 = y;
    y = y0-L(:,N2)./Lp;
    
end

%% Elevation nodes and weigths.
% x = (a*(1-y)+b*(1+y))/2;
t = pi*(1+y)/2;

% w = (b-a)./((1-y.^2).*Lp.^2)*(N2/N1)^2;
w_t = pi./((1-y.^2).*Lp.^2)*(N2/N1)^2;

%% Azimuth nodes.
p = linspace(0, 2*pi, 2*N+1);
p = p(1:end-1);

w_p = pi/N;

%% Grid.
[t, p] = ndgrid(t, p);
w = w_t.*w_p.*sin(t);

end