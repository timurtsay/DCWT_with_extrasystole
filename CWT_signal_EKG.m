%% =========================================================================
%% модуль 3 - CWT 
%% =========================================================================
clear; clc; close all;

%% 0. НАСТРОЙКИ
fprintf([repmat('=', 1, 70), '\n']);
fprintf('ВЕЙВЛЕТ-АНАЛИЗ ЭКГ: ОПТИМИЗАЦИЯ 10σ\n');
%fprintf([repmat('=', 1, 70), '\n\n');

m = 7;
Delta_x = 0.707;
Delta_R = 1.5;
nu_break = 2.0;
nu_max = 10.0;
fs_low = 50;
fs_high = 10;
filename = 'rr_final_105_1280sec.csv';

%% 1. ЗАГРУЗКА ДАННЫХ 
fprintf('1. Загрузка данных...\n');
data = readtable(filename);
t_L = data{:, 1};
t_L = sort(t_L(:));
tau_L = 0.02 * ones(size(t_L));
T = t_L(end);
t_start = t_L(1);
N = length(t_L);
fprintf('   Ударов: %d, T = %.3f с\n\n', N, T);

%% 2. ГРАНИЧНЫЕ ИНТЕРВАЛЫ 
fprintf('2. Расчёт границ...\n');
nu_min = (m + 4 * Delta_R) / T;
t_L_border = 2 * Delta_x / nu_min;  
t_R_border = T - t_L_border;
t_work_start = t_start + t_L_border;
t_work_end = t_R_border;
fprintf('   Рабочий интервал: [%.2f, %.2f] с\n\n', t_work_start, t_work_end);

%% 3. ЧАСТОТНАЯ СЕТКА 
fprintf('3. Частотная сетка...\n');
dlog_critical = 0.005;
dlog_mid = 0.01;
dlog_high = 0.03;

nu_critical = logspace(log10(nu_min), log10(0.4), ...
    round((log10(0.4) - log10(nu_min)) / dlog_critical) + 1);
nu_mid = logspace(log10(0.4), log10(2.0), ...
    round((log10(2.0) - log10(0.4)) / dlog_mid) + 1);
nu_high_viz = logspace(log10(2.0), log10(10.0), ...
    round((log10(10.0) - log10(2.0)) / dlog_high));

nu_save = unique([nu_critical, nu_mid]);
nu_viz = unique([nu_save, nu_high_viz]);
n_save = length(nu_save);
n_viz = length(nu_viz);
n_high_viz = n_viz - n_save;
fprintf('   nu_save: %d точек, nu_viz: %d точек\n\n', n_save, n_viz);

%% 4. ВРЕМЕННАЯ СЕТКА 
fprintf('4. Временная сетка...\n');
T_duration = T - t_start;
n_t_save = round(T_duration * fs_low) + 1;
t_save_grid = linspace(t_start, T, n_t_save);
n_t_viz = round(T_duration * fs_high) + 1;
t_viz_grid = linspace(t_start, T, n_t_viz);
fprintf('   t_save: %d точек, t_viz: %d точек\n\n', n_t_save, n_t_viz);

%% 5. КОНСТАНТЫ
fprintf('5. Константы...\n');
pi2 = pi^2;
sqrt_pi = sqrt(pi);
Dm = 1/sqrt(sqrt_pi*(1 - 2*exp(-3*pi2) + exp(-4*pi2)));
fprintf('   Dm = %.6f\n\n', Dm);

%% 6. ВЫЧИСЛЕНИЕ CWT С ОКНОМ 8σ 
fprintf('6. Вычисление CWT (окно 8σ)...\n');

% 6.1. ЗОНА СОХРАНЕНИЯ (0-2 Гц) 
fprintf('6.1. Низкочастотная зона (0-2 Гц)...\n');
tic_save = tic;

V_save_complex = zeros(length(t_save_grid), n_save, 'like', 1i);

for f = 1:n_save
    nu = nu_save(f);
    
    for L = 1:N
        tau_current = tau_L(L);
        a_l = sqrt(1 + 2 * nu^2 * tau_current^2);
        sigma = a_l / nu;
        window_half = 5 * sigma;
        
        t_center = t_L(L);
        t_min = t_center - window_half;
        t_max = t_center + window_half;
        
        idx_start = find(t_save_grid >= t_min, 1, 'first');
        idx_end = find(t_save_grid <= t_max, 1, 'last');
        
        if isempty(idx_start) || isempty(idx_end) || idx_end < idx_start
            continue;
        end
        
        idx_start = max(1, idx_start);
        idx_end = min(length(t_save_grid), idx_end);
        
        t_window = t_save_grid(idx_start:idx_end).';
        
        b_l = 2 * tau_current * sqrt_pi;
        D = Dm * nu * b_l / a_l;
        E = exp( - nu^2 * (t_window - t_L(L)).^2 / (2 * a_l^2) );
        A = exp( -2 * pi2 * (a_l^2 - 1) / a_l^2 );
        C_L = D * E * A;
        phi_L = -2 * pi * (t_window - t_L(L)) * nu / (a_l^2);
        alpha_L = exp( -2 * pi2 / a_l^2 );
        V_L = C_L .* (exp(1i * phi_L) - alpha_L);
        
        V_save_complex(idx_start:idx_end, f) = V_save_complex(idx_start:idx_end, f) + V_L;
    end
    
    if mod(f, 10) == 0
        fprintf('       %d/%d частот\n', f, n_save);
    end
end

Power_save_full = abs(V_save_complex).^2;
clear V_save_complex;
fprintf('   Готово за %.2f с\n\n', toc(tic_save));

% 6.2. ЗОНА ВИЗУАЛИЗАЦИИ (2-10 Гц)
fprintf('6.2. Высокочастотная зона (2-10 Гц)...\n');

if n_high_viz > 0
    tic_viz = tic;
    
    nu_high_actual = nu_viz(n_save+1:end);
    V_viz_complex = zeros(length(t_viz_grid), length(nu_high_actual), 'like', single(1i));
    
    Dm_s = single(Dm);
    pi2_s = single(pi2);
    t_L_s = single(t_L);
    tau_L_s = single(tau_L);
    
    for f = 1:length(nu_high_actual)
        nu = single(nu_high_actual(f));
        
        for L = 1:N
            tau_current = tau_L_s(L);
            a_l = sqrt(single(1) + single(2) * tau_current^2 * nu^2);
            sigma = a_l / nu;
            window_half = single(4) * sigma;
            
            t_center = t_L_s(L);
            t_min = t_center - window_half;
            t_max = t_center + window_half;
            
            idx_start = find(t_viz_grid >= t_min, 1, 'first');
            idx_end = find(t_viz_grid <= t_max, 1, 'last');
            
            if isempty(idx_start) || isempty(idx_end) || idx_end < idx_start
                continue;
            end
            
            idx_start = max(1, idx_start);
            idx_end = min(length(t_viz_grid), idx_end);
            
            t_window = t_viz_grid(idx_start:idx_end).';
            
            b_l = single(2) * tau_current * sqrt(single(pi));
            D = Dm_s * nu * b_l / a_l;
            E = exp( - nu^2 * (t_window - t_L_s(L)).^2 / (single(2) * a_l.^2) );
            A = exp( -single(2) * pi2_s * (a_l.^2 - single(1)) / a_l.^2 );
            C_L = D .* E .* A;
            phi_L = -single(2) * single(pi) .* (t_window - t_L_s(L)) * nu / (a_l.^2);
            alpha_L = exp( -single(2) * pi2_s / a_l.^2 );
            V_L = C_L .* (exp(single(1i) * phi_L) - alpha_L);
            
            V_viz_complex(idx_start:idx_end, f) = V_viz_complex(idx_start:idx_end, f) + V_L;
        end
        
        if mod(f, 10) == 0
            fprintf('       %d/%d частот\n', f, length(nu_high_actual));
        end
    end
    
    Power_viz_high = abs(V_viz_complex).^2;
    clear V_viz_complex;
    fprintf('   Готово за %.2f с\n\n', toc(tic_viz));
else
    Power_viz_high = [];
end

%% 7. ОБРЕЗКА КРАЁВ 
fprintf('7. Обрезка по рабочему интервалу...\n');

idx_work = (t_save_grid >= t_work_start) & (t_save_grid <= t_work_end);
t_save_final = t_save_grid(idx_work);
Power_save_final = Power_save_full(idx_work, :);

fprintf('   До обрезки: %d точек\n', length(t_save_grid));
fprintf('   После обрезки: %d точек\n', length(t_save_final));
fprintf('   Границы: [%.2f, %.2f] с\n\n', t_work_start, t_work_end);

%% 8. СОХРАНЕНИЕ 
fprintf('8. Сохранение...\n');
save_filename = 'CWT_105.mat';
save(save_filename, 'Power_save_final', 't_save_final', 'nu_save', ...
    'nu_min', 'nu_break', 't_L_border', 't_R_border', ...
    'm', 'Delta_x', 'Delta_R', 'T', 't_start', 'fs_low', '-v7.3');
fprintf('   Сохранено: %s (%.1f МБ)\n\n', save_filename, dir(save_filename).bytes/1024^2);

%% 9. ВИЗУАЛИЗАЦИЯ
fprintf('9. Визуализация...\n');

if ~isempty(Power_viz_high)
    Power_viz_high_interp = zeros(length(t_save_grid), size(Power_viz_high, 2));
    for f = 1:size(Power_viz_high, 2)
        Power_viz_high_interp(:, f) = interp1(t_viz_grid, double(Power_viz_high(:, f)), ...
                                               t_save_grid, 'pchip');
    end
    
    Power_viz_full = zeros(length(t_save_grid), n_viz);
    Power_viz_full(:, 1:n_save) = Power_save_full;
    Power_viz_full(:, n_save+1:end) = Power_viz_high_interp;
    Power_viz_display = Power_viz_full(idx_work, :);
else
    Power_viz_display = Power_save_final;
end

figure('Position', [100, 100, 900, 600]);
surf(nu_viz, t_save_final, Power_viz_display, 'EdgeColor', 'none');
xlabel('\nu, Гц', 'FontSize', 12);
ylabel('t, с', 'FontSize', 12);
zlabel('|V(\nu,t)|^2', 'FontSize', 12);
title('Power CWT (визуализация: 0-10 Гц)', 'FontSize', 14);
colorbar; camlight; view([60, 85]);
xline(nu_break, '--w', 'LineWidth', 1, 'Label', '2 Гц', 'LabelVerticalAlignment', 'bottom');
saveas(gcf, 'CWT_polnoe_105.png');

figure('Position', [100, 100, 900, 600]);
surf(nu_save, t_save_final, Power_save_final, 'EdgeColor', 'none');
xlabel('\nu, Гц', 'FontSize', 12);
ylabel('t, с', 'FontSize', 12);
zlabel('|V(\nu,t)|^2', 'FontSize', 12);
title('Power CWT (сохранённые данные: 0-2 Гц)', 'FontSize', 14);
colorbar; camlight; view([60, 85]);
saveas(gcf, 'CWT_2_105.png');

fprintf('   Графики сохранены.\n\n');

%% 10. СТАТИСТИКА 
fprintf('%s\n', repmat('=', 1, 70));
fprintf('ИТОГИ\n');
fprintf('%s\n', repmat('=', 1, 70));
fprintf('Обрезка краёв: [%.2f, %.2f] с (как в оригинале)\n', t_work_start, t_work_end);
fprintf('Окно 10σ: для каждого удара и частоты\n');
fprintf('Комплексная сумма: V = ΣV_L, потом Power = |V|²\n');
fprintf('Файл: %s\n', save_filename);
fprintf('%s\n', repmat('=', 1, 70));