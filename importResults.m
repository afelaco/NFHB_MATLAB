clc
clear all
close all

%% Path.
path = uigetdir('data');
addpath(path)

%% Files list.
list = dir(fullfile(path, '*.txt'));

%% Import text files.
for n = 1 : length(list)
    
    import = readmatrix(list(n).name);
    [~, name, ~] = fileparts(list(n).name);
    
    % Grid.
    x = unique(import(:, 1))./1000;
    y = unique(import(:, 2))./1000;
    z = unique(import(:, 3))./1000;
    
    [X, Y, Z] = meshgrid(x, y, z);
    
    Ex = reshape(import(:, 4) + 1i.*import(:, 5), size(X));
    Ex = permute(Ex, [2 1 3]);
    
    Ey = reshape(import(:, 6) + 1i.*import(:, 7), size(X));
    Ey = permute(Ey, [2 1 3]);
    
    Ez = reshape(import(:, 8) + 1i.*import(:, 9), size(X));
    Ez = permute(Ez, [2 1 3]);
    
    field = cat(4, Ex, Ey, Ez);
    
    % Save.
    efield = struct('field', field, 'X', X, 'Y', Y, 'Z', Z);
    
    save(fullfile(path, [name, '.m']), 'efield');
    
end