#модуль 2 - расчет границ локальной частоты

import numpy as np
import matplotlib.pyplot as plt

filename = 'rr_final_105_1280sec.csv'

# ------------------------------------------------------------
# 1. ЗАГРУЗКА RR-ИНТЕРВАЛОВ
# ------------------------------------------------------------
try:
    data = np.loadtxt(filename, delimiter=',', skiprows=1)
    t_load = data[:, 0]  # Время (с)
    RR_ms = data[:, 1]  # RR-интервалы (мс)

    print(f"✓ Файл загружён: {filename}")
    print(f"✓ Количество RR-интервалов: {len(RR_ms)}")
    print(f"✓ Диапазон RR: [{np.min(RR_ms):.0f}; {np.max(RR_ms):.0f}] мс")
except Exception as e:
    print(f"❌ Ошибка: {e}")
    raise SystemExit

# ------------------------------------------------------------
# 2. ПАРАМЕТРЫ ИЗ СТАТЬИ (как в MATLAB-коде)
# ------------------------------------------------------------
m = 7
Delta_x = 0.707
Delta_R = 1.5
fs_low = 50  # ← КРИТИЧЕСКИ ВАЖНО: шаг 1/50 = 0.02 с

RR_s = RR_ms / 1000.0
T_total = t_load[-1] - t_load[0]  # Общая длительность записи
t_start = t_load[0]

# ------------------------------------------------------------
# 3. РАСЧЁТ РАБОЧЕГО ИНТЕРВАЛА (как в MATLAB!)
# ------------------------------------------------------------
nu_min_global = (m + 4 * Delta_R) / T_total
t_L_border = 2 * Delta_x / nu_min_global
t_R_border = T_total - t_L_border
t_work_start = t_start + t_L_border
t_work_end = t_start + t_R_border  # ← Важно: + t_start!

print(f"✓ Параметры обрезки (как в MATLAB):")
print(f"  T_total = {T_total:.3f} с")
print(f"  nu_min_global = {nu_min_global:.5f} Гц")
print(f"  t_L_border = {t_L_border:.3f} с")
print(f"  Рабочий интервал: [{t_work_start:.2f}, {t_work_end:.2f}] с")

# ------------------------------------------------------------
# 4. СОЗДАНИЕ СЕТКИ ГРАНИЦ = СЕТКЕ СКЕЙЛОГРАММЫ
# ------------------------------------------------------------
# Шаг точно такой же, как в MATLAB: 1/fs_low = 0.02 с
dt_save = 1.0 / fs_low
t_boundary = np.arange(t_work_start, t_work_end + 1e-9, dt_save)

print(f"\n✓ Сетка границ СОВПАДАЕТ со скейлограммой:")
print(f"  Шаг: {dt_save:.4f} с (1/{fs_low} Гц)")
print(f"  Точек: {len(t_boundary)}")
print(f"  Диапазон: [{t_boundary[0]:.2f}, {t_boundary[-1]:.2f}] с")

# ------------------------------------------------------------
# 5. ВЫЧИСЛЕНИЕ ДИСКРЕТНЫХ ГРАНИЦ (ИСПРАВЛЕНО)
# ------------------------------------------------------------
B_min, B_max = 0.001, 0.21
RR_cr, tau_RR = 0.84, 0.12
A_normal, A_comp = 0.25, 0.5
ECTOPIC_THRESH = 0.5

f_n = 1.0 / RR_s


t_n = t_load.copy()


def compute_B(RR):
    return B_min + (B_max - B_min) / (1.0 + np.exp(-(RR - RR_cr) / tau_RR))


B_vals = compute_B(RR_s)
is_ecopic = RR_s < ECTOPIC_THRESH

f_max = np.zeros_like(f_n)
f_min = np.zeros_like(f_n)
A_vals = np.full_like(f_n, A_normal)

for i in range(len(RR_s)):
    # Относительный критерий паузы (работает с вашими данными)
    if i > 0 and is_ecopic[i - 1] and RR_s[i] > RR_s[i - 1] * 1.3:
        A_vals[i] = A_comp
    f_max[i] = (1.0 + B_vals[i]) * f_n[i]
    f_min[i] = (1.0 - A_vals[i]) * f_n[i]


# ------------------------------------------------------------
# 6. СИГМОИДАЛЬНАЯ АППРОКСИМАЦИЯ (ИСПРАВЛЕНО)
# ------------------------------------------------------------
def sigmoid_approx(f_boundary, t_n, t_eval, RR_s, is_ecopic):
    N = len(f_boundary)
    F_cont = np.full_like(t_eval, f_boundary[0], dtype=float)

    for n in range(1, N):
        # Для экстрасистолы центр перехода = точное время R-зубца из CSV
        if is_ecopic[n - 1]:
            t_c = t_n[n]
        else:
            t_c = (t_n[n] + t_n[n - 1]) / 2.0

        tau_n = RR_s[n - 1] / 6
        delta_f = f_boundary[n] - f_boundary[n - 1]
        sigmoid = 1.0 / (1.0 + np.exp(-(t_eval - t_c) / tau_n))
        F_cont += delta_f * sigmoid
    return F_cont


# Передаём is_ecopic в функцию
F_upper = sigmoid_approx(f_max, t_n, t_boundary, RR_s, is_ecopic)
F_lower = sigmoid_approx(f_min, t_n, t_boundary, RR_s, is_ecopic)
# ------------------------------------------------------------
# 7. ЖЁСТКОЕ ОГРАНИЧЕНИЕ ПО ЧАСТОТЕ (как в статье)
# ------------------------------------------------------------
NU_MAX_CWT1 = 2.0
nu_min_cwt = max(13.0 / T_total, 0.04)  # Формула из статьи

F_upper = np.clip(F_upper, nu_min_cwt, NU_MAX_CWT1)
F_lower = np.clip(F_lower, nu_min_cwt, NU_MAX_CWT1)

print(f"\n✓ Частотный диапазон: [{nu_min_cwt:.4f}; {NU_MAX_CWT1:.2f}] Гц")
print(f"✓ Верхняя граница: [{F_upper.min():.4f}, {F_upper.max():.4f}] Гц")
print(f"✓ Нижняя граница: [{F_lower.min():.4f}, {F_lower.max():.4f}] Гц")

# ------------------------------------------------------------
# 8. СОХРАНЕНИЕ (БЕЗ ИНТЕРПОЛЯЦИИ!)
# ------------------------------------------------------------
# Верхняя граница
output_file = 'upper_boundary_105.txt'
data_to_save = np.column_stack((t_boundary, F_upper))
np.savetxt(
    output_file,
    data_to_save,
    fmt='%.6f',
    delimiter='\t',
    header='time(s)\tF_upper(Hz)',
    comments='',
    encoding='utf-8'
)
print(f"\n✓ Верхняя граница: {len(t_boundary)} точек → {output_file}")

# Нижняя граница
output_file = 'lower_boundary_105.txt'
data_to_save = np.column_stack((t_boundary, F_lower))
np.savetxt(
    output_file,
    data_to_save,
    fmt='%.6f',
    delimiter='\t',
    header='time(s)\tF_lower(Hz)',
    comments='',
    encoding='utf-8'
)
print(f"✓ Нижняя граница: {len(t_boundary)} точек → {output_file}")

# ------------------------------------------------------------
# 9. ВИЗУАЛИЗАЦИЯ
# ------------------------------------------------------------
plt.figure(figsize=(14, 5))

# Частотный коридор
plt.fill_between(t_boundary, F_lower, F_upper, color='lightblue', alpha=0.4, label='Частотный коридор')
plt.plot(t_boundary, F_upper, 'b--', linewidth=2.0, label='$F^{>}(t)$ (верхняя)')
plt.plot(t_boundary, F_lower, 'r--', linewidth=2.0, label='$F^{<}(t)$ (нижняя)')

# Дискретные точки (для контекста)
plt.scatter(t_n, f_n, c='gray', s=15, alpha=0.5, zorder=2, label='$f_n = 1/RR_n$')

plt.xlabel('Время, с', fontsize=12, fontweight='bold')
plt.ylabel('Частота, Гц', fontsize=12, fontweight='bold')
plt.title('Границы локальной частоты (сетка = скейлограмме, шаг 0.02 с)',
          fontsize=13, fontweight='bold', pad=10)
plt.grid(True, alpha=0.3, linestyle='--')
plt.legend(fontsize=10, loc='upper right')
plt.xlim(t_boundary[0], t_boundary[-1])
plt.tight_layout()

plt.show()

# ------------------------------------------------------------
# 10. ДИАГНОСТИКА СОВПАДЕНИЯ СЕТОК
# ------------------------------------------------------------
print("\n" + "=" * 70)
print("ДИАГНОСТИКА: СОВПАДЕНИЕ СЕТОК ГРАНИЦ И СКЕЙЛОГРАММЫ")
print("=" * 70)
print(f"Шаг границ:       {np.mean(np.diff(t_boundary)):.6f} с")
print(f"Ожидаемый шаг:    {1 / fs_low:.6f} с (1/50 Гц)")
print(f"Разница:          {abs(np.mean(np.diff(t_boundary)) - 1 / fs_low):.2e} с")
print("   Просто используйте:")
print("      for i in range(len(t_boundary)):")
print("          nu_min = F_lower[i]")
print("          nu_max = F_upper[i]")
print("          p_vec = Power[:, i]  # Power[частота, время]")
print("=" * 70)