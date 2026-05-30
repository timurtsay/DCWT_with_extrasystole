clear; clc; close all;

%% 1. ЗАГРУЗКА ДАННЫХ
filename = 'peaks_123.txt';
data = readtable(filename);
t_raw = data{:, 1};
nu_HRmax_raw = data{:, 2};

%% 2. ВЫБОР ИНТЕРВАЛА
t_zoom_start = 413;    
t_zoom_end   = 433;  
idx_zoom = (t_raw >= t_zoom_start) & (t_raw <= t_zoom_end);
t_zoom = t_raw(idx_zoom);
nu_max_zoom = nu_HRmax_raw(idx_zoom);

%% 3. СОЗДАНИЕ ФИГУРЫ
fig = figure('Units', 'centimeters', ...
             'Position', [1 1 15 12], ...
             'Color', 'white');

ax = axes('Parent', fig);
hold(ax, 'on');

% Основная линия
plot(ax, t_zoom, nu_max_zoom, 'k-', 'LineWidth', 1.2, 'HandleVisibility', 'off');

% Вертикальная линия PVC
PVC = xline(ax, 334.635, '--k', 'LineWidth', 1.5, ...
            'HandleVisibility', 'on', 'DisplayName', 'PVC – 334.5с');

% Вместо ylabel используем text
text(ax, 0.03, 1, '\nu_{max}(t)', ...  % -0.05 выносит за пределы оси
     'Units', 'normalized', ...
     'HorizontalAlignment', 'right', ...
     'VerticalAlignment', 'bottom', ...
     'FontSize', 11, ...
     'Interpreter', 'tex');

% Подпись оси X - справа внизу
text(ax, 1, -0.05, 't, с', ...  % 1 - правый край, -0.05 - ниже оси
     'Units', 'normalized', ...
     'HorizontalAlignment', 'right', ...
     'VerticalAlignment', 'top', ...
     'FontSize', 11);

% Легенда
legend(ax, 'show', ...
       'Location', 'northeast', ...
       'FontSize', 11, ...              % Согласовано с основным текстом
       'Box', 'off', ...
       'Interpreter', 'tex');

% Оси
xlim(ax, [t_zoom(1), t_zoom(end)]);
ylim(ax, [0, max(nu_max_zoom) * 1.15]);

ax.Box = 'off';

% Настройка шрифтов
ax.FontSize = 9;                      % Подписи делений
ax.TickLength = [0.02 0.01];

% Отступы
ax.InnerPosition = [0.15 0.20 0.82 0.67];

% Сохранение
exportgraphics(fig, 'loqf_123_413_433.png', 'Resolution', 600, 'BackgroundColor', 'white');

%close(fig);