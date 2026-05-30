%% ============================================================
%% 3D ВИЗУАЛИЗАЦИЯ СЕГМЕНТА СКЕЙЛОГРАММЫ 
%% ============================================================
clear; clc; close all;

%% 1. ПАРАМЕТРЫ
t_start_custom = 527;           % Начало сегмента (с)
t_end_custom   = 627;           % Конец сегмента (с)
freq_max_display = 10;          % Макс. частота для отображения (Гц)
           
%% 2. ЗАГРУЗКА ДАННЫХ 
fprintf('Загрузка данных...\n');
load('dcwt_power_102.mat');

t_full = t_dcwt;           
nu = nu_dcwt;                    
Power_full = Power_DW;   

%% 3. ОПРЕДЕЛЕНИЕ РАЗМЕРНОСТИ
[n_rows, n_cols] = size(Power_full);
len_t = length(t_full);
len_nu = length(nu);

if len_t == n_rows && len_nu == n_cols
    Power_for_surf = Power_full;  % [time × freq]
    dim_time = 1;
elseif len_t == n_cols && len_nu == n_rows
    Power_for_surf = Power_full';  % [time × freq]
    dim_time = 2;
else
    error('Не удалось определить размерность данных!');
end

%% 4. ВЫДЕЛЕНИЕ СЕГМЕНТА 
idx_time = (t_full >= t_start_custom) & (t_full <= t_end_custom);
idx_freq = (nu <= freq_max_display);

if sum(idx_time) == 0
    error('Пустой временной сегмент. Доступно: [%.2f, %.2f] с', min(t_full), max(t_full));
end

t_segment = t_full(idx_time);
nu_viz = nu(idx_freq);
Power_viz = Power_for_surf(idx_time, :);
Power_viz_display = Power_viz(:, idx_freq);

fprintf('Сегмент: %.3f — %.3f с (%d отсчётов)\n', t_segment(1), t_segment(end), length(t_segment));
fprintf('Частоты: %.2f — %.2f Гц (%d точек)\n', nu_viz(1), nu_viz(end), length(nu_viz));

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
zlabel(ax, '|V_{DCWT}(\nu,t)|^2', 'FontSize', 11);


% Настройка осей
ax.FontSize = 9;
ax.Box = 'on';
ax.ZTickLabel = {};  % Убираем метки на оси Z

% Угол обзора
view([60, 85])
axis tight
% Отступы (как на предыдущих графиках)
ax.InnerPosition = [0.15 0.20 0.82 0.67];

% Освещение
camlight(ax, 'headlight');
lighting(ax, 'gouraud');
material(ax, 'dull');

% Colormap
colormap(ax, 'parula')


% Сетка
grid(ax, 'on');

% Сохранение 
exportgraphics(fig, 'dcwt_527_627_102.png', 'Resolution', 600, 'BackgroundColor', 'white');

fprintf('Визуализация сохранена: dcwt_527_627_102.png.png\n');