clc
clear all
close all

%% Save option.
saveOption = questdlg('Save data?');

switch saveOption
    case 'Yes'
        saveOption = 1;
    case 'No'
        saveOption = 0;
    case 'Cancel'
        return
end

%% Global.
global k

acc = -40;

%% Import antenna.
Selpath = uigetdir('data');
addpath(Selpath)

tag = strsplit(Selpath, '\');
tag = strsplit(tag{end}, '_');

%% Import frequency.
f = str2double(tag{1});
w = physconst('Lightspeed')/f;
k = 2*pi/w;

%% Import array layout.
layout = stlread('layout.stl');
R = max(vecnorm(layout.Points, 2, 2));
L = ceil(k*R);

figure
hold on
trisurf(layout)
grid on
view([1 1 1])

%% Import ports.
import = xlsread('ports.xlsx');

N = length(import);

x_n = import(:, 1);
y_n = import(:, 2);
z_n = import(:, 3);

[r_n, t_n, p_n] = cart2sph_c(x_n, y_n, z_n);

for n = 1 : N
   
    R_n(n) = max(vecnorm(layout.Points - [x_n(n), y_n(n), z_n(n)], 2, 2));
    
end

L_n = ceil(k.*R_n);

%% Focus point.
t_f = 0;
p_f = 0;
r_f = 3*w;

if r_f < R
    
    error('Focusing too close!');
    
end

%% Polarization.
u = [0 1 0];

%% Focal shift.
a = pi.^2.*R.^4 - 45.*w.^2.*r_f.^2;
b = -4.*pi.^2.*R.^4.*r_f;
c = 3.*pi.^2.*R.^4.*r_f.^2;

r_f_eff = -(sqrt(b.^2 - 4.*a.*c) + b)./(2.*a);

[x_f_eff, y_f_eff, z_f_eff] = sph2cart_c(r_f_eff, t_f, p_f);

[x_f, y_f, z_f] = sph2cart_c(r_f, t_f, p_f);

title(strcat('\theta_f=', sprintf('%.0f°', rad2deg(t_f)), {' '}, ...
    '\phi_f=', sprintf('%.0f°', rad2deg(p_f)), {' '}, ...
    '{\it r}_{f,eff}=', sprintf('%.1f', r_f_eff/w), '\lambda', {' '}, ...
    '{\it r}_f=', sprintf('%.1f', r_f/w), '\lambda'))

scatter3(x_f_eff, y_f_eff, z_f_eff, 'filled', ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r')
axis equal

%%
% r_f = linspace(0, 10*w, 1000);
% 
% a = pi.^2.*R.^4 - 45.*w.^2.*r_f.^2;
% b = -4.*pi.^2.*R.^4.*r_f;
% c = 3.*pi.^2.*R.^4.*r_f.^2;
% 
% r_f_eff = -(sqrt(b.^2 - 4.*a.*c) + b)./(2.*a);
% 
% figure
% hold on
% plot(r_f./w, r_f_eff./w, 'LineWidth', 3)
% xline(2.5, 'k:', 'linewidth', 2)
% yline(.3053, 'k:', 'linewidth', 2)
% yline(pi*R^2/sqrt(15)/w^2, 'k:', 'linewidth', 2)
% hold off
% xlabel('$r_\text{f} [\lambda]$', 'Interpreter', 'none')
% ylabel('$r_\text{f,eff} [\lambda]$', 'Interpreter', 'none')
% grid on
% set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.5)
% 
% gca_hold = gca;
% gca_hold.XRuler.TickLabelGapOffset = 5;
% gca_hold.YRuler.TickLabelGapOffset = 10;

%% Closest element.
[~, idx] = min(sqrt((x_f-x_n).^2 + (y_f-y_n).^2 + (z_f-z_n).^2));

%% Circumscribing sphere.
[X, Y, Z] = sphere;
surf(X.*R, Y.*R, Z.*R, 'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.1)

clear X Y Z

tag = sprintf('x=%.3f_y=%.3f_z=%.3f', x_f_eff, y_f_eff, z_f_eff);

%% Import far-fields.
progressBar = waitbar(0, 'Importing far-fields... 0%', 'Position', [441  150  270   56.25]);

F = cell(N, 1);

for n = 1 : N
    
    % Import far-field.
    import = readmatrix(strcat(Selpath, strcat('\port_', sprintf('%.f', n), '.txt')));
    
    [F{n}.t, F{n}.p] = ndgrid(deg2rad(unique(import(:,1))), deg2rad(unique(import(:,2))));
    
    F{n}.F(:,:,1) = reshape(import(:,4).*exp(1i.*deg2rad(import(:,5))), size(F{n}.t));
    F{n}.F(:,:,2) = reshape(import(:,6).*exp(1i.*deg2rad(import(:,7))), size(F{n}.t));
    
    quiver3(x_n(n), y_n(n), z_n(n), x_f_eff-x_n(n), y_f_eff-y_n(n), z_f_eff-z_n(n), ...
        'r', 'LineStyle', ':', 'ShowArrowHead', 'off', 'Marker', 'o')
    text(x_n(n), y_n(n), z_n(n), sprintf('%.0f', n))
    
    waitbar(n/N, progressBar, strcat('Importing far-fields...', {' '}, sprintf('%0.f%%', n/N*100)));
    
end

hold off
close(progressBar)

clear import

%% Folders.
if saveOption == 1
    
    imagesFolder = strcat(Selpath, '\images\', tag);
    if ~exist(imagesFolder, 'dir')
        mkdir(imagesFolder)
    end
    
    resultsFolder = strcat(Selpath, '\results\', tag);
    if ~exist(resultsFolder, 'dir')
        mkdir(resultsFolder)
    end
    
end

%% Translation.
progressBar = waitbar(0, 'Translating fields... 0%', 'Position', [441 150 270 56.25]);

err = zeros(N, 1);
pow = zeros(N, 1);
pow_t = zeros(N, 1);
L_t = zeros(size(L_n));

figure
for n = 1 : N

    while err(n) > acc

        % SSHA.
        f = SSHA(F{n}, L_n(n), "lebedev", "cartesian");
        SFT{n} = f;
        
        % Error check.
        err(n) = max(mag2db(abs(F{n}.F - SSHS(f, F{n}.t, F{n}.p, "polar"))), [], 'all');
        
        waitbar(n/N, progressBar, strcat('Field', {' '}, sprintf('%d%', n), ...
            {','}, {' '},'L =', {' '}, sprintf('%0.f%', L_n(n)), {','}, {' '}, ...
            'err =  ', {' '}, sprintf('%.1f%', err(n)), {' '}, 'dB,', {' '}, ...
            'ratio =  ', {' '}, sprintf('%.1f%', pow_t(n)/pow(n)*100), {' '}, '%'));
        
        % Iteration.
        if err(n) > acc
            
            L_n(n) = L_n(n)+1;
        
        end
        
    end

    % Projection.
    f_co = sum(f.*conj(u), 2);
    
    % Power.
    pow(n) = sum(abs(f_co).^2);
    
    L_t(n) = L;
    
    % Translation.
    while pow_t(n) < pow(n)
        
        T = translator(r_n(n), t_n(n), p_n(n), L_t(n), L_n(n));
        
        f_co_t = T*f_co;
        
        pow_t(n) = sum(abs(f_co_t).^2);
        
        waitbar(n/N, progressBar, strcat('Field', {' '}, sprintf('%d%', n), ...
            {','}, {' '},'L =', {' '}, sprintf('%0.f%', L_n(n)), {','}, {' '}, ...
            'err =  ', {' '}, sprintf('%.1f%', err(n)), {' '}, 'dB,', {' '}, ...
            'ratio =  ', {' '}, sprintf('%.1f%', pow_t(n)/pow(n)*100), {' '}, '%'));
        
        if pow_t(n) < pow(n)
           
            L_t(n) = L_t(n)+1;
 
        end
        
    end
    
    if pow_t(n) > pow(n)
        
        L_t(n) = L_t(n)-1;
        
        T = translator(r_n(n), t_n(n), p_n(n), L_t(n), L_n(n));
        
        f_co_t = T*f_co;
        
        pow_t(n) = sum(abs(f_co_t).^2);
        
    end

    M(1:(length(f_co)), n) = f_co;
    M_t(1:(length(f_co_t)), n) = f_co_t;

    I(n,1) = exp(-1i.*angle(SSWS(f_co_t, r_f, t_f, p_f)));
    
    % Plot.
    subplot(2,1,1)
    imagesc(mag2db(abs(M)))
    
    subplot(2,1,2)
    imagesc(mag2db(abs(M_t)))
    
end

close(progressBar)

I = I./I(idx);

if saveOption == 1
    
    writematrix(round(rad2deg(angle(I)), 3), ...
        strcat(resultsFolder, '\current_', tag, '.xlsx'))
    
end

%% Current PCM.
I_PCM = exp(1i.*k.*(sqrt((x_f - x_n(:)).^2 + (y_f - y_n(:)).^2 + (z_f - z_n(:)).^2)));
I_PCM = I_PCM./I_PCM(idx);

if saveOption == 1
    
    writematrix(round(rad2deg(angle(I_PCM)), 3), ...
        strcat(resultsFolder, '\current_PCM_', tag, '.xlsx'))
    
end

%% Current QPA.
% I_QPA = exp(1i.*k./r_f.*(r_n.^2./2 - (x_n.*x_f + y_n.*y_f + z_n.*z_f)));
% I_QPA = I_QPA./I_QPA(idx);
% 
% if saveOption == 1
%     
%     writematrix(round(rad2deg(angle(I_QPA)), 3), ...
%         strcat(resultsFolder, '\current_QPA_', tag, '.xlsx'))
%     
% end

%% Current phase profile.
figure
hold on
plot(angle(I), ':o', 'LineWidth', 2)
plot(angle(I_PCM), ':o', 'LineWidth', 2)
hold off
xlim([1 N])
ylim([-pi pi])
xticks(1:N)
grid on
xlabel('Antennas')
ylabel('Phase')
yticks([-pi, -pi/2, 0, pi/2, pi])
yticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'})
legend('HPC', 'PC')

set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.5)

if saveOption == 1
    
    saveas(gcf, strcat(imagesFolder, '\current_phase.png'))
    
end

%% Solution.
viewCoefficients(M_t*I)