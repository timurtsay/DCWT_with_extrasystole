%% ===== ГРАФИК СИГНАЛА Z(t) =====
clear; clc; close all;

%% ===== ПАРАМЕТРЫ =====
t_L = [419.189, 420.464, 421.806, 423.008, 423.531, 424.269, 425.383, 426.542, 427.681, 428.878];

L_values = 1:length(t_L);
tau_L = 0.02 * ones(size(t_L));  % Все tau_L равны 0.02

%% ===== АНАЛИЗ РАЗНОСТЕЙ =====
differences = diff(t_L);

[max_diff, max_index] = max(differences);
[min_diff, min_index] = min(differences);

fprintf('Все разницы: %s\n', mat2str(differences, 3));
fprintf('Максимальная разница: %.2f\n', max_diff);
fprintf('Между элементами %.3f и %.3f\n', t_L(max_index), t_L(max_index + 1));
fprintf('Минимальная разница: %.2f\n', min_diff);
fprintf('Между элементами %.3f и %.3f\n', t_L(min_index), t_L(min_index + 1));

min_f_max = 1/max_diff;
max_f_max = 1/min_diff;

fprintf('Минимальная локальная частота: %.2f Гц\n', min_f_max);
fprintf('Максимальная локальная частота: %.2f Гц\n', max_f_max);

%% ===== ФУНКЦИЯ z_L(t) =====
z_L = @(t, t_L_j, tau_L_j) exp(-(t - t_L_j).^2 ./ (4 * tau_L_j^2));

%% ===== СОЗДАНИЕ МАССИВА ВРЕМЕНИ =====
t = linspace(419, 429, 3000);  % Диапазон времени

%% ===== ВЫЧИСЛЕНИЕ Z(t) =====
Z = zeros(size(t));
for j = 1:length(L_values)
    Z = Z + z_L(t, t_L(j), tau_L(j));
end

%% ===== ПОСТРОЕНИЕ ГРАФИКА =====
fig = figure('Units', 'centimeters', ...
             'Position', [1 1 9 7], ...
             'Color', 'white');

ax = axes('Parent', fig);
hold(ax, 'on');

plot(ax, t, Z, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Z(t)');

grid(ax, 'off');
xlabel(ax, 'Время t, с', 'FontSize', 11);
ylabel(ax, 'Z(t)', 'FontSize', 11);




ax.FontSize = 9;
ax.Box = 'off';
ax.InnerPosition = [0.15 0.20 0.82 0.67];

% Сохранение
%exportgraphics(fig, 'Z_signal.png', 'Resolution', 600, 'BackgroundColor', 'white');

fprintf('\u2713 Визуализация сохранена: Z_signal.png\n');