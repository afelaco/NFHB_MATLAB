function viewCoefficients(c)

if iscell(c) == 1
   
    c = cell2mat(c);
    
end

%% Degree and order vectors.
L = sqrt(length(c))-1;
i = 0:(L+1)^2-1;
l = floor(sqrt(i));

dim = size(c, 2);

for i = 0 : L
    for j = 1 : dim
        
        C(i+1, L-i+1:L+i+1, j) = c(l == i, j);
        
    end
end

l = 0:L;
m = -L:L;

figure
for i = 1 : dim
    
    subplot(dim, 1, i)
    imagesc(m,l,mag2db(abs(C(:,:,i))),'AlphaData',l' >= abs(m))
    axis equal
    xlim([-L-0.5 L+0.5])
    ylim([-0.5 L+0.5])
    caxis([min(mag2db(abs(c)),[],'all') max(mag2db(abs(c)),[],'all')])
    colorbar
    set(gca, 'visible', 'off')
    
end

end