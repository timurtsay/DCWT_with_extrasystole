%% =========================================================================
%% модуль 5 - DCWT 
%% =========================================================================

clear; clc; close all;

%% 0. ЗАГРУЗКА ТРАЕКТОРИИ ЧАСТОТЫ (ν_HRmax(t))
filename = 'peaks_105.txt';
data = readtable(filename);
t_raw = data{:, 1};
nu_HRmax_raw = data{:, 2};
% Сортировка и удаление дубликатов
if ~issorted(t_raw)
    [t_raw, idx] = sort(t_raw);
    nu_HRmax_raw = nu_HRmax_raw(idx);
end
[t_raw, unique_idx] = unique(t_raw);
nu_HRmax_raw = nu_HRmax_raw(unique_idx);

% Проверка равномерности сетки
dt_actual = median(diff(t_raw));
fprintf('Загружено %d точек, dt = %.3f с, T = %.1f с\n', ...
    length(t_raw), dt_actual, t_raw(end) - t_raw(1));

%% 1. ПАРАМЕТРЫ

fprintf('\n=== ПАРАМЕТРЫ DCWT ===\n');

% Вейвлет Морле
Omega0 = 2 * pi;
Dm = 1 / sqrt(sqrt(pi) * (1 - 2*exp(-3*Omega0^2/4) + exp(-Omega0^2)));
correction = exp(-0.5 * Omega0^2);

% Границы (формулы из статьи, стр. 36)
m = 7;
Delta_B = 1.5;
Delta_x = 0.707;

T_total = t_raw(end) - t_raw(1);
nu_min = (m + 4 * Delta_B) / T_total;
t_L_border = 2 * Delta_x / nu_min;
t_R_border = T_total - t_L_border;

fprintf('ν_min = %.5f Гц\n', nu_min);
fprintf('Границы: t_L = t_R = %.2f с\n', t_L_border);
fprintf('Рабочий интервал: [%.2f, %.2f] с (%.1f%%)\n\n', ...
    t_raw(1) + t_L_border, t_raw(1) + t_R_border, ...
    100*(t_R_border - t_L_border)/T_total);

% Логарифмическая частотная сетка

nu_max_dcwt = 0.4;
dlog_nu = 0.01;  % шаг по log10(ν), как в CWT
N_nu = round((log10(nu_max_dcwt) - log10(nu_min)) / dlog_nu) + 1;
nu_dcwt = logspace(log10(nu_min), log10(nu_max_dcwt), N_nu);  

% Временная сетка (используем шаг 0.6 с)
dt_dcwt = dt_actual;  % 0.6 с — не меняем
idx_work = (t_raw >= t_raw(1) + t_L_border) & ...
           (t_raw <= t_raw(1) + t_R_border);
t_dcwt = t_raw(idx_work);
nu_signal = nu_HRmax_raw(idx_work);

fprintf('Сетка DCWT: %d частот × %d временных точек (логарифмическая)\n', ...
    N_nu, length(t_dcwt));
fprintf('Шаг по времени: %.3f с\n\n', dt_dcwt);


% Инициализация параллельного пула (4 ядра) 

fprintf('=== ИНИЦИАЛИЗАЦИЯ ПАРАЛЛЕЛЬНЫХ ВЫЧИСЛЕНИЙ ===\n');
parpool('local', 4);  % Запуск пула на 4 ядра
fprintf('Параллельный пул активен: 4 рабочих\n\n');


%% 2. DCWT (БЕЗ КРАЕВЫХ ЭФФЕКТОВ)

fprintf('=== ВЫЧИСЛЕНИЕ DCWT (ПАРАЛЛЕЛЬНОЕ, 4 ЯДРА) ===\n');
tic;

% Используем полный временной массив для вычислений
N_t_full = length(t_raw);
V_DW_full = zeros(N_nu, N_t_full);

% parfor для параллелизации на 4 ядра 
parfor i = 1:N_nu
    nu_val = nu_dcwt(i);
    
    % Ширина вейвлета 5σ
    sigma_t = Delta_x / nu_val;
    n_sigma = 5;  % 5σ как в CWT
    if nu_val < 0.05
        n_sigma = 6;  % для низких частот ещё шире
    end
    window_width = 2 * n_sigma * sigma_t;
    
    % Пропускаем частоты, где вейвлет шире всего сигнала
    if window_width > (t_raw(end) - t_raw(1))
        continue;
    end
    
    V_col = zeros(N_t_full, 1);
    
    % Вычисляем для всез точек времени (включая краевые зоны)
    for j = 1:N_t_full
        t_center = t_raw(j);
        
        % Окно интегрирования ±5σ 
        t_win_min = t_center - n_sigma*sigma_t;
        t_win_max = t_center + n_sigma*sigma_t;
        
        % Находим точки данных в окне
        idx_win = (t_raw >= t_win_min) & (t_raw <= t_win_max);
        
        % Если точек мало - пропускаем (будет ноль)
        if sum(idx_win) < 3
            continue;
        end
        
        t_win = t_raw(idx_win);
        s_win = nu_HRmax_raw(idx_win);
        
        
        % Адаптивная интерполяция для интегрирования 
        
        % Оцениваем требуемую плотность сетки: >= 10 точек на периоде осцилляции
        T_osc = 2*pi / (Omega0 * nu_val);  % период осцилляции вейвлета
        N_min = ceil((t_win(end) - t_win(1)) / (T_osc / 10));
        N_interp = max(N_min, 8 * length(t_win));  % минимум 8× исходных точек
        
        t_fine = linspace(t_win(1), t_win(end), N_interp);
        s_fine = interp1(t_win, s_win, t_fine, 'pchip');  % PCHIP для гладкости
        
        % Вейвлет-функция на мелкой сетке
        x = nu_val * (t_fine - t_center);
        gauss_term = exp(-0.5 * x.^2);
        osc_term = exp(-1i * Omega0 * x) - correction;
        psi = Dm * sqrt(nu_val) * gauss_term .* osc_term;
        
        % Интегрирование: метод Симпсона для осциллирующих функций
        integrand = s_fine .* conj(psi);
        if mod(N_interp, 2) == 1 && N_interp >= 5
            % Симпсон: точность O(h^4)
            h_step = (t_fine(end) - t_fine(1)) / (N_interp - 1);
            S = integrand(1) + integrand(end) + ...
                4*sum(integrand(2:2:end-1)) + ...
                2*sum(integrand(3:2:end-2));
            V_col(j) = nu_val * (h_step/3) * S;
        else
            % Fallback на trapz, если чётное число точек
            V_col(j) = nu_val * trapz(t_fine, integrand);
        end
    end
    
    V_DW_full(i, :) = V_col.';
end


% 2.1. Устранение граничных эффектов (Обнуление краев)

% Обнуляем края после вычисления (как в статье, стр. 36)
V_DW = V_DW_full;  % Копируем рабочий вариант

% Находим индексы границ на полной сетке
idx_left  = (t_raw < t_raw(1) + t_L_border);
idx_right = (t_raw > t_raw(1) + t_R_border);

% Обнуляем значения в краевых зонах для ВСЕХ частот
V_DW(:, idx_left)  = 0;
V_DW(:, idx_right) = 0;

% Вырезаем рабочий интервал для визуализации (только для графиков)
idx_work = (t_raw >= t_raw(1) + t_L_border) & ...
           (t_raw <= t_raw(1) + t_R_border);
t_dcwt = t_raw(idx_work);
nu_signal = nu_HRmax_raw(idx_work);
V_DW_work = V_DW(:, idx_work);  % Используем обнуленную матрицу

elapsed_time = toc;
fprintf('Вычисление завершено за %.2f секунд\n', elapsed_time);


%% 3. МОЩНОСТЬ DCWT (без спектральных интегралов)

% Локальная плотность энергии спектра (Ур. 5, стр. 36)
% ε(ν,t) = (2/C_ψ) × |V(ν,t)|² / ν
Power_DW = abs(V_DW_work).^2;

% 3.1. СОХРАНЕНИЕ МАТРИЦЫ МОЩНОСТИ ДЛЯ СПЕКТРАЛЬНЫХ ИНТЕГРАЛОВ

fprintf('\n=== СОХРАНЕНИЕ ДАННЫХ ===\n');

% Сохраняем только необходимые данные для расчёта спектральных интегралов
save('dcwt_105.mat', ...
    'Power_DW', ...    % |V_DW|^2 — матрица мощности [N_nu × N_t]
    'nu_dcwt', ...     % частотная ось [1 × N_nu], Гц
    't_dcwt');         % временная ось [1 × N_t], с

fprintf('Сохранено: dcwt_105.mat (Power_DW, nu_dcwt, t_dcwt)\n');


%% 4. ВИЗУАЛИЗАЦИЯ (3D СКАЙЛОГРАММА DCWT)

fprintf('\n=== ПОСТРОЕНИЕ 3D СКАЙЛОГРАММЫ ===\n');

% Создаём фигуру с правильным размером
figure('Position', [50, 50, 1200, 700], 'Color', 'white');

% Основной 3D график
h_surf = surf(nu_dcwt, t_dcwt, Power_DW.', 'EdgeColor', 'none');
% Настройка цветов
colormap("parula");  % Или: parula, hot, turbo
cbar = colorbar;
%cbar.Label.String = '|V_{DW}|^2';
cbar.Label.FontSize = 11;
caxis([0, max(Power_DW(:))]);  % Фиксированный диапазон для сравнения

% Подписи осей
xlabel('Частота \nu, Гц', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('Время t, с', 'FontSize', 13, 'FontWeight', 'bold');
zlabel('|V_{DСWT}(\nu,t)|^2', 'FontSize', 13, 'FontWeight', 'bold');
title('Мощность DCWT. Запись 105', 'FontSize', 15, 'FontWeight', 'bold');

% Настройка вида
view([60, 75])
camlight('headlight');  % Источник света
lighting gouraud;  % Плавное затенение
shading interp;  % Интерполяция цветов

% Ограничения осей
xlim([nu_min, nu_max_dcwt]);
ylim([t_dcwt(1), t_dcwt(end)]);
zlim([0, max(Power_DW(:))]);

% Сетка для ориентации
grid on;
set(gca, 'GridAlpha', 0.3);

% Сохранение
saveas(gcf, 'dcwt_scalogram_105.png', 'png');
fprintf('3D скейлограмма сохранена: dcwt_scalogram_105.png\n');

%% 5. СТАТИСТИКА

fprintf('=== СТАТИСТИКА ===\n');
fprintf('Диапазон частот: [%.5f, %.3f] Гц\n', nu_min, nu_max_dcwt);
fprintf('Число частот: %d (логарифмическая сетка)\n', N_nu);
fprintf('Число временных точек: %d\n', length(t_dcwt));
fprintf('Шаг по времени: %.3f с\n', dt_dcwt);
fprintf('Граничные интервалы: %.2f с\n', t_L_border);
fprintf('Время вычисления: %.2f с\n', elapsed_time);

[max_val, max_idx] = max(Power_DW(:));
[row, col] = ind2sub(size(Power_DW), max_idx);
fprintf('Пик мощности: %.4e (t = %.2f с, ν = %.4f Гц)\n', ...
    max_val, t_dcwt(col), nu_dcwt(row));

%% 6. ОЧИСТКА ПАРАЛЛЕЛЬНОГО ПУЛА

delete(gcp('nocreate'));  % Закрытие пула после вычислений
fprintf('\nDCWT ЗАВЕРШЁН (4 ядра)\n');