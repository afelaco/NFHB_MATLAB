function viewFarfield(F, t, p)

figure

for i = 1 : size(F, 3)
   
    maxF = max(abs(F), [], 'all');
    
end

for i = 1 : size(F, 3)
    
    subplot(2, size(F, 3), i)
    [x, y, z] = sph2cart_c(abs(F(:,:,i)), t, p);
    c = abs(F(:,:,i));
    surf(x, y, z, c, 'edgecolor', 'none')
    colorbar
    xlim([-maxF maxF])
    ylim([-maxF maxF])
    zlim([-maxF maxF])
    caxis([0 maxF])
    view([1 1 1])
    axis equal
    
    subplot(2, size(F, 3), size(F, 3)+i)
    [x, y, z] = sph2cart_c(1, t, p);
    c = rad2deg(angle(F(:,:,i)));
    surf(x, y, z, c, 'edgecolor', 'none')
    colorbar
    view([1 1 1])
    axis equal
    
end

end