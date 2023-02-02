function output = GBC(a, k)

output = gamma(a+1)./gamma(a-k+1)./factorial(k);

end