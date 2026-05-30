%% =========================================================================
%% модуль 6 - РАсчет спектральных интегралов 
%% =========================================================================
clear; clc; close all;

%% 1. ЗАГРУЗКА ДАННЫХ
load('dcwt_power_105.mat');  % Power_DW, nu_dcwt, t_dcwt

fprintf('=== ЗАГРУЗКА ДАННЫХ ===\n');
fprintf('Power_DW: [%d × %d]\n', size(Power_DW,1), size(Power_DW,2));
fprintf('nu_dcwt: %d точек, [%g; %g] Гц\n', ...
    length(nu_dcwt), nu_dcwt(1), nu_dcwt(end));
fprintf('t_dcwt: %d точек, [%g; %g] с\n', ...
    length(t_dcwt), t_dcwt(1), t_dcwt(end));

t_norm_start = 250;   % начало интервала [с] — измените под свои нужды
t_norm_end   = 1000;  % конец интервала [с] — измените под свои нужды
fprintf('Нормировка по интервалу: [%.1f; %.1f] с\n', t_norm_start, t_norm_end);
%% 2. ПАРАМЕТРЫ
C_psi = 1.0132;  % Константа нормировки вейвлета Морле (статья BSPC 2014)

% Границы частотных диапазонов (Гц) — ДОБАВЛЕН VLF
ranges = struct( ...
    'VLF', [0.015, 0.04], ...  % очень низкие частоты
    'LF',  [0.04,  0.15], ...  % низкие частоты
    'HF',  [0.15,  0.40]  ...  % высокие частоты
);

%% 3. ПЛОТНОСТЬ ЭНЕРГИИ ε(ν,t) — Ур. (11)
% ε(ν,t) = (2/C_ψ) * |V_DW|² / ν
epsilon = (2 / C_psi) * (Power_DW ./ nu_dcwt');  % [N_nu × N_t]

E = struct();
for name = {'VLF', 'LF', 'HF'}
    name = name{1};
    nu_low = ranges.(name)(1);
    nu_high = ranges.(name)(2);
    
    idx = (nu_dcwt >= nu_low) & (nu_dcwt <= nu_high);
    
    % Интегрирование по d(ln ν)
    % ∫ ε(ν) dν = ∫ [ε(ν)·ν] d(ln ν), т.к. d(ln ν) = dν/ν
    ln_nu = log(nu_dcwt(idx));
    integrand = epsilon(idx, :) .* nu_dcwt(idx)';  % |V|² · (2/C_ψ)
    
    % Трапеции по равномерной сетке ln(ν)
    E.(name) = trapz(ln_nu, integrand, 1) / (log(nu_high) - log(nu_low));
    
    fprintf('E_%s(t): mean = %.3e, std = %.3e\n', ...
        name, mean(E.(name)), std(E.(name)));
end
%% 5. ПАРАМЕТРЫ d(t) И ОТНОШЕНИЙ МЕЖДУ ДИАПАЗОНАМИ
fprintf('\n=== РАСЧЁТ ПАРАМЕТРОВ d(t) ===\n');

for name = {'VLF', 'LF', 'HF'}
    name = name{1};
    
    % Находим индексы времени в нормировочном окне
    idx_norm = (t_dcwt >= t_norm_start) & (t_dcwt <= t_norm_end);
    
    % Проверка: есть ли точки в интервале?
    if ~any(idx_norm)
        error('Нет точек времени в интервале нормировки [%.2f, %.2f] с!', ...
            t_norm_start, t_norm_end);
    end
    
    % Вычисляем среднее ТОЛЬКО на этом интервале
    E_mean_norm = mean(E.(name)(idx_norm));
    
    % Нормируем весь сигнал на это значение
    d.(name) = E.(name) ./ E_mean_norm;
    
    fprintf('d_%s(t): mean = %.3f, max = %.3f, <E_%s>_norm = %.3e (на интервале [%.1f; %.1f] с)\n', ...
        name, mean(d.(name)), max(d.(name)), name, E_mean_norm, ...
        t_norm_start, t_norm_end);
end

% d_μ(t) = E_μ(t) / <E_μ> для всех трёх диапазонов
% for name = {'VLF', 'LF', 'HF'}
%     name = name{1};
%     E_mean = mean(E.(name));
%     d.(name) = E.(name) ./ E_mean;
%     fprintf('d_%s(t): mean = %.3f, max = %.3f, <E_%s> = %.3e\n', ...
%         name, mean(d.(name)), max(d.(name)), name, E_mean);
% end

% Параметры отношений между диапазонами (обратный порядок)
fprintf('\n=== РАСЧЁТ ПАРАМЕТРОВ ОТНОШЕНИЙ ===\n');

% d_VLF/HF(t) = [E_VLF(t)/E_HF(t)] / <[E_VLF/E_HF]>
ratio_VLF_HF = E.VLF ./ E.HF;
ratio_VLF_HF_mean = mean(E.VLF(idx_norm)) ./ mean(E.HF(idx_norm));
d_VLF_HF = ratio_VLF_HF ./ ratio_VLF_HF_mean;
fprintf('d_VLF/HF(t): mean = %.3f, max = %.3f, <VLF/HF> = %.3f\n', ...
    mean(d_VLF_HF(idx_norm)), max(d_VLF_HF), ratio_VLF_HF_mean);

% d_VLF/LF(t) = [E_VLF(t)/E_LF(t)] / <[E_VLF/E_LF]>
ratio_VLF_LF = E.VLF ./ E.LF;
ratio_VLF_LF_mean = mean(E.VLF(idx_norm))./mean(E.LF(idx_norm));
d_VLF_LF = ratio_VLF_LF ./ ratio_VLF_LF_mean;
fprintf('d_VLF/LF(t): mean = %.3f, max = %.3f, <VLF/LF> = %.3f\n', ...
    mean(d_VLF_LF(idx_norm)), max(d_VLF_LF), ratio_VLF_LF_mean);

% d_LF/HF(t) = [E_LF(t)/E_HF(t)] / <[E_LF/E_HF]>  (классический вагосимпатический индекс)
ratio_LF_HF = E.LF./E.HF  ;
ratio_LF_HF_mean = mean(E.LF(idx_norm))./mean(E.HF(idx_norm));
d_LF_HF = ratio_LF_HF ./ ratio_LF_HF_mean;
fprintf('d_LF/HF(t): mean = %.3f, max = %.3f, <LF/HF> = %.3f\n', ...
    mean(d_LF_HF(idx_norm)), max(d_LF_HF), ratio_LF_HF_mean);

% Отдельные переменные
d_VLF = d.VLF;
d_LF  = d.LF;
d_HF  = d.HF;
E_VLF_mean = mean(E.VLF(idx_norm));
E_LF_mean  = mean(E.LF(idx_norm));
E_HF_mean  = mean(E.HF(idx_norm));

%% 6. СОХРАНЕНИЕ РЕЗУЛЬТАТОВ
fprintf('\n=== СОХРАНЕНИЕ ===\n');

% Приводим все векторы к столбцам [N_t × 1]
t_col = t_dcwt(:);
E_VLF_col = E.VLF(:);
E_LF_col  = E.LF(:);
E_HF_col  = E.HF(:);
d_VLF_col = d_VLF(:);
d_LF_col  = d_LF(:);
d_HF_col  = d_HF(:);
d_VLF_HF_col = d_VLF_HF(:);  
d_VLF_LF_col = d_VLF_LF(:);  
d_LF_HF_col  = d_LF_HF(:);   

% Проверка согласованности размеров
assert(length(t_col) == length(E_VLF_col), 'Размер t и E_VLF не совпадает!');
assert(length(E_VLF_col) == length(E_LF_col), 'Размер E_VLF и E_LF не совпадает!');
assert(length(E_LF_col) == length(E_HF_col), 'Размер E_LF и E_HF не совпадает!');
assert(length(E_HF_col) == length(d_VLF_col), 'Размер E_HF и d_VLF не совпадает!');
assert(length(d_VLF_col) == length(d_LF_col), 'Размер d_VLF и d_LF не совпадает!');
assert(length(d_LF_col) == length(d_HF_col), 'Размер d_LF и d_HF не совпадает!');
assert(length(d_HF_col) == length(d_VLF_HF_col), 'Размер d_HF и d_VLF_HF не совпадает!');
assert(length(d_VLF_HF_col) == length(d_VLF_LF_col), 'Размер d_VLF_HF и d_VLF_LF не совпадает!');
assert(length(d_VLF_LF_col) == length(d_LF_HF_col), 'Размер d_VLF_LF и d_LF_HF не совпадает!');

% Сохранение в .mat (добавляем новые параметры)
save('dcwt_E_d_results_105.mat', ...
    'E', 'd_VLF', 'd_LF', 'd_HF', ...
    'd_VLF_HF', 'd_VLF_LF', 'd_LF_HF', ...  
    't_dcwt', 'nu_dcwt', ...
    'E_VLF_mean', 'E_LF_mean', 'E_HF_mean', ...
    'ratio_VLF_HF_mean', 'ratio_VLF_LF_mean', 'ratio_LF_HF_mean', ...  
    'ranges', 'C_psi');
fprintf('Сохранено: dcwt_E_d_results_105.mat\n');

% Экспорт в CSV (добавляем новые столбцы)
export_table = table(t_col, E_VLF_col, E_LF_col, E_HF_col, ...
                     d_VLF_col, d_LF_col, d_HF_col, ...
                     d_VLF_HF_col, d_VLF_LF_col, d_LF_HF_col, ...  
    'VariableNames', {'time_s', 'E_VLF', 'E_LF', 'E_HF', ...
                      'd_VLF', 'd_LF', 'd_HF', ...
                      'd_VLF_HF', 'd_VLF_LF', 'd_LF_HF'});  
writetable(export_table, 'dcwt_HF_LF_VLF_d_results_105.csv');
fprintf('Экспортировано: dcwt_HF_LF_VLF_d_results_105.csv\n');

%% 7. ВИЗУАЛИЗАЦИЯ
% 7.4 Сводный график: d_VLF, d_LF, d_HF на одном полотне
figure('Position', [100, 100, 900, 400], 'Color', 'white');

plot(t_dcwt, d_VLF, 'g-', 'LineWidth', 1.5, 'DisplayName', 'd_{VLF}(t)');
hold on;
plot(t_dcwt, d_LF,  'b-', 'LineWidth', 1.5, 'DisplayName', 'd_{LF}(t)');
plot(t_dcwt, d_HF,  'r-', 'LineWidth', 1.5, 'DisplayName', 'd_{HF}(t)');

% Базовая линия d = 1
yline(1, '--k', 'LineWidth', 1, 'DisplayName', 'Средний уровень');

grid on;
xlabel('Время t, с');
ylabel('d_{\mu}(t)');
%title('Параметры d_{VLF}, d_{LF}, d_{HF} (нормированные на среднее)');
legend('Location', 'best', 'FontSize', 9);
xlim([t_dcwt(1), t_dcwt(end)]);

sgtitle('Параметры d_{HF, LF, VLF}. Запись 105', ...
    'FontSize', 10, 'FontWeight', 'bold');
saveas(gcf, 'd_param_105.png');
fprintf('Сохранено: d_param_105.png\n');

fprintf('\nРАСЧЁТ ЗАВЕРШЁН\n');