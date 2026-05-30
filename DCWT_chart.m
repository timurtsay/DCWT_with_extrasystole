%% ============================================================
%% 3D ВИЗУАЛИЗАЦИЯ СЕГМЕНТА СКЕЙЛОГРАММЫ 
%% ============================================================
clear; clc; close all;

%%  1. ПАРАМЕТРЫ 
t_start_custom = 284;           
t_end_custom   = 384;           
freq_max_display = 10;          

%% 2. ЗАГРУЗКА ДАННЫХ 
fprintf('Загрузка данных...\n');
load('dcwt_power_202.mat');

t_full = t_dcwt;           
nu = nu_dcwt;                    
Power_full = Power_DW;   

%% 3. ОПРЕДЕЛЕНИЕ РАЗМЕРНОСТИ 
[n_rows, n_cols] = size(Power_full);
len_t = length(t_full);
len_nu = length(nu);

if len_t == n_rows && len_nu == n_cols
    Power_for_surf = Power_full;
elseif len_t == n_cols && len_nu == n_rows
    Power_for_surf = Power_full';
else
    error('Не удалось определить размерность данных!');
end

%% 4. ВЫДЕЛЕНИЕ СЕГМЕНТА 
idx_time = (t_full >= t_start_custom) & (t_full <= t_end_custom);
idx_freq = (nu <= freq_max_display);

t_segment = t_full(idx_time);
nu_viz = nu(idx_freq);
Power_viz = Power_for_surf(idx_time, :);
Power_viz_display = Power_viz(:, idx_freq);

%% 5. 3D ВИЗУАЛИЗАЦИЯ 
fig = figure('Units', 'centimeters', ...
             'Position', [1 1 14 12], ...
             'Color', 'white');

ax = axes('Parent', fig);
hold(ax, 'on');

[T_mesh, NU_mesh] = meshgrid(t_segment, nu_viz);

surf(ax, NU_mesh, T_mesh, Power_viz_display', ...
    'EdgeColor', 'none', ...
    'FaceLighting', 'gouraud');

% Подписи осей 
xlabel(ax, '\nu, Гц', 'FontSize', 11);
ylabel(ax, 't, с', 'FontSize', 11);
zlabel(ax, '|V(\nu,t)|^2', 'FontSize', 11);

% Настройка осей
ax.FontSize = 9;
ax.Box = 'on';
grid on;
ax.ZTickLabel = {};  % Убираем метки на оси Z

% Угол обзора - частота слева, время справа
view([60, 75])

xRange = range(nu_viz);
yRange = range(t_segment);
zRange = range(Power_viz_display(:));

xlim(ax, [min(nu_viz), max(nu_viz)]);
ylim(ax, [min(t_segment), max(t_segment)]);
zlim(ax, [min(Power_viz_display(:)), max(Power_viz_display(:))]);

% Освещение
camlight(ax, 'headlight');
lighting(ax, 'gouraud');
material(ax, 'dull');

% Автоматическая подстройка осей под данные
axis tight

% Сохранение
exportgraphics(fig, 'dcwt_202_284_384.png', 'Resolution', 600, 'BackgroundColor', 'white');

fprintf('Визуализация сохранена\n');