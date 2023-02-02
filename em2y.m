function y = em2y(e, m)
%% Order.
L = sqrt(length(e))-1;

tag = 'lib\em2y';

if isfile(tag) == 0
    
    M = em2yMatrix(L);
    save(tag, 'M')
    
else
    
    load(tag)
    
    if M.L < L
        
        M = em2yMatrix(L);
        save(tag, 'M')
        
    elseif M.L > L
        
        M.e = M.e(1:(L+1)^2, 1:(L+1)^2, :);
        M.m = M.m(1:(L+1)^2, 1:(L+1)^2, :);
        
    end
    
end

y_sph = squeeze(sum(M.e.*e.' + M.m.*m.', 2));

y(:,1) = 1/sqrt(2).*(y_sph(:,2) - y_sph(:,1));
y(:,2) = -1i/sqrt(2).*(y_sph(:,1) + y_sph(:,2));
y(:,3) = y_sph(:,3);

end