%% ============================================================
%% ГРАФИК ПАРАМЕТРОВ d НА ЗАДАННОМ ИНТЕРВАЛЕ
%% ============================================================
clear; clc; close all;

%% 1. НАСТРОЙКИ 
results_file = 'dcwt_HF_LF_VLF_d_results_105.mat';  % Файл с данными
t_zoom_start = 530;    % Начало интервала (с)
t_zoom_end   = 570;    % Конец интервала (с)
save_fig     = true;   % Сохранять ли графики

%эктрасистола
vline_x           = 550.96;       % позиция линии по X
vline_style       = '--k';         % стиль линии
vline_label       = 'PVC – 550.96c';% подпись в легенде

%% 2. ЗАГРУЗКА И ВЫБОРКА 
load(results_file);  % Загружает: d_*, t_dcwt, nu_dcwt, ...

idx_zoom = (t_dcwt >= t_zoom_start) & (t_dcwt <= t_zoom_end);
t_zoom = t_dcwt(idx_zoom);

% Параметры d для отдельных диапазонов
d_VLF_zoom = d_VLF(idx_zoom);
d_LF_zoom  = d_LF(idx_zoom);
d_HF_zoom  = d_HF(idx_zoom);

% Параметры отношений (если есть в файле)
if exist('d_VLF_HF', 'var')
    d_VLF_HF_zoom = d_VLF_HF(idx_zoom);
    d_VLF_LF_zoom = d_VLF_LF(idx_zoom);
    d_LF_HF_zoom  = d_LF_HF(idx_zoom);
    has_ratios = true;
else
    warning('Параметры отношений не найдены в файле %s', results_file);
    has_ratios = false;
end

%% ГРАФИК 1: d_{VLF}(t)
fig1 = figure('Name', 'd_{VLF}', 'Color', 'white', 'Position', [100 100 600 400]);
plot(t_zoom, d_VLF_zoom, 'g-', 'LineWidth', 1.5, 'DisplayName', 'd_{VLF}(t)');
grid on; xlabel('Время t, с'); ylabel('d_{VLF}(t)');
xline(vline_x, vline_style, 'LineWidth', 1.5, 'HandleVisibility', 'on', 'DisplayName', vline_label);
xlim([t_zoom(1), t_zoom(end)]);
legend('show');
%if save_fig, saveas(gcf, 'd_VLF_102.png'); end

%% ГРАФИК 2: d_{LF}(t)
fig2 = figure('Name', 'd_{LF}', 'Color', 'white', 'Position', [100 100 600 400]);
plot(t_zoom, d_LF_zoom, 'b-', 'LineWidth', 1.5, 'DisplayName', 'd_{LF}(t)');
grid on; xlabel('Время t, с'); ylabel('d_{LF}(t)');
xline(vline_x, vline_style, 'LineWidth', 1.5, 'HandleVisibility', 'on', 'DisplayName', vline_label);
xlim([t_zoom(1), t_zoom(end)]);
legend('show');
%if save_fig, saveas(gcf, 'd_LF_102.png'); end

%% ГРАФИК 3: d_{HF}(t)
fig3 = figure('Name', 'd_{HF}', 'Color', 'white', 'Position', [100 100 600 400]);
plot(t_zoom, d_HF_zoom, 'r-', 'LineWidth', 1.5, 'DisplayName', 'd_{HF}(t)');
grid on; xlabel('Время t, с'); ylabel('d_{HF}(t)');
xline(vline_x, vline_style, 'LineWidth', 1.5, 'HandleVisibility', 'on', 'DisplayName', vline_label);
xlim([t_zoom(1), t_zoom(end)]);
legend('show');
%if save_fig, saveas(gcf, 'd_HF_102.png'); end

%% НОВЫЙ ГРАФИК 4: все три параметра вместе (d_VLF, d_LF, d_HF)
fig4 = figure('Name', 'd_all', 'Color', 'white', 'Position', [100 100 600 400]);
hold on;
plot(t_zoom, d_VLF_zoom, 'g-', 'LineWidth', 1.5, 'DisplayName', 'd_{VLF}(t)');
plot(t_zoom, d_LF_zoom,  'b-', 'LineWidth', 1.5, 'DisplayName', 'd_{LF}(t)');
plot(t_zoom, d_HF_zoom,  'r-', 'LineWidth', 1.5, 'DisplayName', 'd_{HF}(t)');
grid on; xlabel('Время t, с'); ylabel('d(t)');
xline(vline_x, vline_style, 'LineWidth', 1.5, 'HandleVisibility', 'on', 'DisplayName', vline_label);
xlim([t_zoom(1), t_zoom(end)]);
legend('show', 'Location', 'best');
box on;
hold off;
if save_fig, saveas(gcf, 'd_all_105.png'); end

%% ГРАФИК 5: d_{LF/HF}(t) — если есть
if has_ratios
    fig5 = figure('Name', 'd_{LF/HF}', 'Color', 'white', 'Position', [100 100 600 400]);
    plot(t_zoom, d_LF_HF_zoom, '-', 'Color', '#A2142F', 'LineWidth', 1.5, 'DisplayName', 'd_{LF/HF}(t)');
    grid on; xlabel('Время t, с'); ylabel('d_{LF/HF}(t)');
    xline(vline_x, vline_style, 'LineWidth', 1.5, 'HandleVisibility', 'on', 'DisplayName', vline_label);
    yline(1, '--k', 'LineWidth', 1.5, 'HandleVisibility', 'on', 'DisplayName', '1');
    xlim([t_zoom(1), t_zoom(end)]);
    legend('show');
    %if save_fig, saveas(gcf, 'd_LF_HF_102.png'); end
end

%% ГРАФИК 6: d_{VLF/HF}(t) — если есть
if has_ratios
    fig6 = figure('Name', 'd_{VLF/HF}', 'Color', 'white', 'Position', [100 100 600 400]);
    plot(t_zoom, d_VLF_HF_zoom, '-', 'Color', '#A2142F', 'LineWidth', 1.5, 'DisplayName', 'd_{VLF/HF}(t)');
    grid on; xlabel('Время t, с'); ylabel('d_{VLF/HF}(t)');
    xline(vline_x, vline_style, 'LineWidth', 1.5, 'HandleVisibility', 'on', 'DisplayName', vline_label);
    yline(1, '--k', 'LineWidth', 1.5, 'HandleVisibility', 'on', 'DisplayName', '1');
    xlim([t_zoom(1), t_zoom(end)]);
    legend('show');
    %if save_fig, saveas(gcf, 'd_VLF_HF_102.png'); end
end

%% ГРАФИК 7: d_{VLF/LF}(t) — если есть
if has_ratios
    fig7 = figure('Name', 'd_{VLF/LF}', 'Color', 'white', 'Position', [100 100 600 400]);
    plot(t_zoom, d_VLF_LF_zoom, '-', 'Color', '#A2142F', 'LineWidth', 1.5, 'DisplayName', 'd_{VLF/LF}(t)');
    grid on; xlabel('Время t, с'); ylabel('d_{VLF/LF}(t)');
    xline(vline_x, vline_style, 'LineWidth', 1.5, 'HandleVisibility', 'on', 'DisplayName', vline_label);
    yline(1, '--k', 'LineWidth', 1.5, 'HandleVisibility', 'on', 'DisplayName', '1');
    xlim([t_zoom(1), t_zoom(end)]);
    legend('show');
    %if save_fig, saveas(gcf, 'd_VLF_LF_102.png'); end
end