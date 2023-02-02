function out = wigner3j(j1, m1, j2, m2, j3, m3)

if all(abs([m1, m2, m3]) <= [j1, j2, j3]) && m1+m2+m3 == 0 && abs(j1-j2) <= j3 && j3 <= j1+j2
    
    % Delta function (https://dlmf.nist.gov/34.2.E5).
    D = sqrt(factorial(j1+j2-j3)*factorial(j1-j2+j3)*factorial(-j1+j2+j3)/factorial(j1+j2+j3+1));
    
    % Wigner 3j symbol (https://dlmf.nist.gov/34.2.E4).
    idx_s = 1;
    
    for s = max([0, j2-j3-m1, j1-j3+m2]) : min([j1+j2-j3, j1-m1, j2+m2])
        
        out(idx_s) = (-1)^(j1-j2-m3+s)*D*...
            sqrt(factorial(j1+m1)*factorial(j1-m1)*factorial(j2+m2)*factorial(j2-m2)*factorial(j3+m3)*factorial(j3-m3))/...
            (factorial(s)*factorial(j1+j2-j3-s)*factorial(j1-m1-s)*factorial(j2+m2-s)*factorial(j3-j2+m1+s)*factorial(j3-j1-m2+s));
        
        idx_s = idx_s + 1;
        
    end
    
    out = sum(out);
    
else
    
    out = 0;
    
end

end