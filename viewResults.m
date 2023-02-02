clc
clear all
close all

%% Path.
selpath = uigetdir('data');
addpath(selpath)

%% Images folder.
imagesFolder = strsplit(selpath, '\');
imagesFolder{end-1} = 'images';
imagesFolder = strjoin(imagesFolder, '\');

if ~exist(imagesFolder, 'dir')
    mkdir(imagesFolder)
end

%% Import frequency.
f = strsplit(selpath, '\');
f = strsplit(f{end-2}, '_');
f = str2double(f{1});

w = physconst('Lightspeed')/f;
k = 2*pi/w;

%% Focus point.
tag = strsplit(selpath, '\');
tag = strsplit(tag{end}, '_');

for i = 1 : 3
   
    eval(tag{i});
    
end

%% Files list.
list = dir(fullfile(selpath, '*.m'));

%% Import fields.
for n = 1 : length(list)
    
    efield = importdata(list(n).name);

    E_co{n} = mag2db(abs(efield.field(:,:,:,2)));
    E_cross{n} = mag2db(sqrt(abs(efield.field(:,:,:,1)).^2 + abs(efield.field(:,:,:,3)).^2));
    
end

X = efield.X;
Y = efield.Y;
Z = efield.Z;

[~, idx_x] = min(abs(unique(X) - x));
[~, idx_y] = min(abs(unique(Y) - y));
[~, idx_z] = min(abs(unique(Z) - z));

idx_z = idx_z + 2;

clear efield

%% Trim.
[~, idx_trim] = min(abs(unique(Z) - 2.5*w));

X = X(:, :, 1:idx_trim+1);
Y = Y(:, :, 1:idx_trim+1);
Z = Z(:, :, 1:idx_trim+1);

for n = 1 : length(list)
    
    E_co{n} = E_co{n}(:, :, 1:idx_trim+1);
    E_cross{n} = E_cross{n}(:, :, 1:idx_trim+1);
    
end

%% Levels.
for n = 1 : length(list)
   
    E_min(n) = min(E_cross{n}(:,:,idx_z), [], 'all');
    E_max(n) = E_co{n}(idx_y, idx_x, idx_z);
    lev{n} = (floor((E_min(n)-E_max(n))/3)*3 : 3 : 0) + E_max(n);
  
end

%% Copolar slice plot.
for n = 1 : length(list)

    [~, name, ~] = fileparts(list(n).name);
    
    extension = '.svg';
    
    E_co{n}(E_co{n} < E_min(n)) = E_min(n);
    
    % XZ Cut.
    figure
    hold on
    [c, h] = contourf(squeeze(X(idx_y,:,:))./w, squeeze(Z(idx_y,:,:))./w, squeeze(E_co{n}(idx_y,:,:)), ...
        lev{n});
    h.LevelList = round(h.LevelList, 1);
%     clabel(c, h, 'manual')
    xline(X(1,idx_x,1)/w, ':', 'LineWidth', 3)
    yline(Z(1,1,idx_z)/w, ':', 'LineWidth', 3)
    ylim([min(Z(:))/w 2.5])
    hold off
    axis equal
    xlh = xlabel('$\small x [\lambda]$', 'interpreter', 'none');
    ylh = ylabel('$\small z [\lambda]$', 'interpreter', 'none');
    
    c = customcolormap([0 .25 .5 .75 1], [1 0 0; 1 1 0; 0 1 0; 0 1 1; 0 0 1]);
    colormap(c)
    caxis([min([lev{1}(:); lev{2}(:)]) max([lev{1}(:); lev{2}(:)])])
    
    gca_hold = gca;
    gca_hold.XRuler.TickLabelGapOffset = 5;
    gca_hold.YRuler.TickLabelGapOffset = 10;
    
    xlh.Position(1) = xlh.Position(1) + .4;
    ylh.Position(2) = ylh.Position(2) + .5;
    
    saveas(gcf, fullfile(imagesFolder, [name, '_co_xz', extension]))
    
    % YZ Cut.
    figure
    hold on
    [c, h] = contourf(squeeze(Y(:,idx_x,:))./w, squeeze(Z(:,idx_x,:))./w, squeeze(E_co{n}(:,idx_x,:)), ...
        lev{n});
    h.LevelList = round(h.LevelList, 1);
%     clabel(c, h, 'manual')
    xline(Y(idx_y,1,1)/w, ':', 'LineWidth', 3)
    yline(Z(1,1,idx_z)/w, ':', 'LineWidth', 3)
    ylim([min(Z(:))/w 2.5])
    hold off
    axis equal
    xlh = xlabel('$\small y [\lambda]$', 'interpreter', 'none');
    ylh = ylabel('$\small z [\lambda]$', 'interpreter', 'none');
    
    c = customcolormap([0 .25 .5 .75 1], [1 0 0; 1 1 0; 0 1 0; 0 1 1; 0 0 1]);
    colormap(c)
    caxis([min([lev{1}(:); lev{2}(:)]) max([lev{1}(:); lev{2}(:)])])
    
    gca_hold = gca;
    gca_hold.XRuler.TickLabelGapOffset = 5;
    gca_hold.YRuler.TickLabelGapOffset = 10;
    
    xlh.Position(1) = xlh.Position(1) + .4;
    ylh.Position(2) = ylh.Position(2) + .5;
    
    saveas(gcf, fullfile(imagesFolder, [name, '_co_yz', extension]))
    
    % XY Cut.
    figure
    hold on
    [c, h] = contourf(X(:,:,idx_z)./w, Y(:,:,idx_z)./w, E_co{n}(:,:,idx_z), ...
        lev{n}, 'ShowText', 'on');
    h.LevelList = round(h.LevelList, 1);
    clabel(c, h, lev{n}(end-3:end))
    [~, h] = contour(X(:,:,idx_z)./w, Y(:,:,idx_z)./w, E_co{n}(:,:,idx_z), ...
        [E_max(n)-3 E_max(n)-3], 'LineColor', 'k', 'LineStyle', ':', 'LineWidth', 5, 'ShowText', 'on');
    h.LevelList = round(h.LevelList, 1);
    xline(X(1,idx_x,1)/w, ':', 'LineWidth', 3)
    yline(Y(idx_y,1,1)/w, ':', 'LineWidth', 3)
    hold off
    axis equal
    xlh = xlabel('$\small x [\lambda]$', 'interpreter', 'none');
    ylh = ylabel('$\small y [\lambda]$', 'interpreter', 'none');
    
    c = customcolormap([0 .25 .5 .75 1], [1 0 0; 1 1 0; 0 1 0; 0 1 1; 0 0 1]);
    colormap(c)
    caxis([min([lev{1}(:); lev{2}(:)]) max([lev{1}(:); lev{2}(:)])])
    
    gca_hold = gca;
    gca_hold.XRuler.TickLabelGapOffset = 5;
    gca_hold.YRuler.TickLabelGapOffset = 10;
    
    xlh.Position(1) = xlh.Position(1) + .2;
    ylh.Position(2) = ylh.Position(2) + .2;
    
    saveas(gcf, fullfile(imagesFolder, [name, '_co_xy', extension]))
    
end

%% Crosspolar slice plot.
for n = 1 : length(list)

    [~, name, ~] = fileparts(list(n).name);
    
    extension = '.svg';

    % XY Cut.
    figure
    hold on
    [c, h] = contourf(X(:,:,idx_z)./w, Y(:,:,idx_z)./w, E_cross{n}(:,:,idx_z), ...
        lev{n}, 'ShowText', 'on');
    h.LevelList = round(h.LevelList, 1);
    clabel(c, h, lev{n}(1:2:end))
    xline(X(1,idx_x,1)/w, ':', 'LineWidth', 3)
    yline(Y(idx_y,1,1)/w, ':', 'LineWidth', 3)
    hold off
    axis equal
    xlh = xlabel('$\small x [\lambda]$', 'interpreter', 'none');
    ylh = ylabel('$\small y [\lambda]$', 'interpreter', 'none');
    
    c = customcolormap([0 .25 .5 .75 1], [1 0 0; 1 1 0; 0 1 0; 0 1 1; 0 0 1]);
    colormap(c)
    caxis([min([lev{1}(:); lev{2}(:)]) max([lev{1}(:); lev{2}(:)])])
    
    gca_hold = gca;
    gca_hold.XRuler.TickLabelGapOffset = 5;
    gca_hold.YRuler.TickLabelGapOffset = 10;
    
    xlh.Position(1) = xlh.Position(1) + .2;
    ylh.Position(2) = ylh.Position(2) + .2;
    
    saveas(gcf, fullfile(imagesFolder, [name, '_cross_xy', extension]))
    
end

%% Focal depth.
figure
hold on
plot(squeeze(Z(idx_y, idx_x, :))/w, squeeze(E_co{1}(idx_y, idx_x, :)), 'LineWidth', 3)
plot(squeeze(Z(idx_y, idx_x, :))/w, squeeze(E_co{2}(idx_y, idx_x, :)), 'LineWidth', 3)
xlabel('$z [\lambda]$', 'interpreter', 'none')
xlim([0.1 2.5])
ylabel('Co-polar field magnitude $[dB V\m]$', 'Interpreter', 'none')
legend('This work', 'Phase-conjugate')

grid on
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.5)

gca_hold = gca;
gca_hold.XRuler.TickLabelGapOffset = 5;
gca_hold.YRuler.TickLabelGapOffset = 10;