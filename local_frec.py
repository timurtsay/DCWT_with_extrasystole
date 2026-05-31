#модуль 4 - расчет локальной частоты

import numpy as np
import matplotlib.pyplot as plt
import h5py

# ------------------------------------------------------------
# 1. НАСТРОЙКИ И ЗАГРУЗКА ДАННЫХ
# ------------------------------------------------------------
print("=" * 60)
print("АНАЛИЗ ТРАЕКТОРИИ ЧАСТОТЫ (ПИК НА ГРАНИЦЕ ТОЛЬКО ПРИ ЭКСТРАСИСТОЛЕ)")
print("=" * 60)

mat_filename = 'CWT_105.mat'
boundary_upper_file = 'upper_boundary_105.txt'
boundary_lower_file = 'lower_boundary_105.txt'
rr_filename = 'rr_final_105_1280sec.csv'  # ← НУЖЕН ФАЙЛ С РЯДОМ RR!

print(f"\n1. Загрузка скейлограммы из {mat_filename}...")
try:
    with h5py.File(mat_filename, 'r') as f:
        Power = f['Power_save_final'][:]  # [частота × время]
        t = f['t_save_final'][:].flatten()
        nu = f['nu_save'][:].flatten()

    print(f"   ✓ Power.shape = {Power.shape}")
    print(f"   Время: t∈[{t.min():.2f}, {t.max():.2f}] с ({len(t)} точек)")
    print(f"   Частота: ν∈[{nu.min():.4f}, {nu.max():.2f}] Гц ({len(nu)} точек)")

except Exception as e:
    print(f"Ошибка загрузки: {e}")
    raise

# ------------------------------------------------------------
# 2. ИДЕНТИФИКАЦИЯ ЭКТРАСИСТОЛ И ПАУЗ ПО РЯДУ RR
# ------------------------------------------------------------
print(f"\n2. Идентификация экстрасистол и компенсаторных пауз...")

try:
    rr_data = np.loadtxt(rr_filename, delimiter=',', skiprows=1)
    t_rr = rr_data[:, 0]
    RR_ms = rr_data[:, 1]
    RR_s = RR_ms / 1000.0

    print(f"   ✓ Загружено {len(RR_s)} RR-интервалов")

    # === КРИТЕРИИ ===
    ECTOPIC_THRESH = 0.5  #экстрасистола
    PAUSE_THRESH = 0.8  # пауза Для разных ритмов можно подобрать разные значения

    is_ectopic = np.zeros(len(RR_s), dtype=bool)
    is_pause = np.zeros(len(RR_s), dtype=bool)

    for i in range(len(RR_s)):

        if RR_s[i] < ECTOPIC_THRESH:
            is_ectopic[i] = True

        if RR_s[i] > PAUSE_THRESH:
            is_pause[i] = True

    print(f"Экстрасистол: {np.sum(is_ectopic)}, Пауз: {np.sum(is_pause)}")

    # === СОЗДАНИЕ МАСОК ДЛЯ СКЕЙЛОГРАММЫ ===
    ectopic_mask = np.zeros_like(t, dtype=bool)
    pause_mask = np.zeros_like(t, dtype=bool)
    WINDOW_HALF = 0.35

    for i in range(len(RR_s)):
        t_curr = t_rr[i]

        if is_ectopic[i]:
            mask = (t >= t_curr - WINDOW_HALF) & (t <= t_curr + WINDOW_HALF)
            ectopic_mask[mask] = True
        if is_pause[i]:
            t_center_pause = t_curr - (RR_s[i] / 2.0)
            mask = (t >= t_center_pause - WINDOW_HALF) & (t <= t_center_pause + WINDOW_HALF)
            pause_mask[mask] = True

    print(f"Точек с экстрасистолой: {np.sum(ectopic_mask)}")
    print(f"Точек с паузой: {np.sum(pause_mask)}")

except Exception as e:
    print(f"Ошибка RR: {e}")
    ectopic_mask = None
    pause_mask = None

# ------------------------------------------------------------
# 3. ЗАГРУЗКА ГРАНИЦ
# ------------------------------------------------------------
print(f"\n3. Загрузка границ...")


def load_boundary(filename):
    try:
        data = np.loadtxt(filename, delimiter='\t', skiprows=1, encoding='utf-8')
    except:
        data = np.loadtxt(filename, skiprows=1, encoding='utf-8')
    return data[:, 0], data[:, 1]


t_upper, nu_max_nodes = load_boundary(boundary_upper_file)
t_lower, nu_min_nodes = load_boundary(boundary_lower_file)

print(f"   Верхняя граница: {len(t_upper)} точек")
print(f"   Нижняя граница: {len(t_lower)} точек")

# ------------------------------------------------------------
# 4. ОБРЕЗКА ГРАНИЦ ПОД СКЕЙЛОГРАММУ
# ------------------------------------------------------------
print(f"\n4. Обрезка границ под скейлограмму...")

mask_upper = (t_upper >= t.min()) & (t_upper <= t.max())
mask_lower = (t_lower >= t.min()) & (t_lower <= t.max())

t_upper = t_upper[mask_upper]
nu_max_nodes = nu_max_nodes[mask_upper]
t_lower = t_lower[mask_lower]
nu_min_nodes = nu_min_nodes[mask_lower]

print(f"   Верхняя граница (обрезана): {len(t_upper)} точек")
print(f"   Нижняя граница (обрезана): {len(t_lower)} точек")

# ------------------------------------------------------------
# 5. СИНХРОНИЗАЦИЯ ВРЕМЕННЫХ СЕТОК И МАСОК
# ------------------------------------------------------------
t_samples = t_upper
nu_max_interp = nu_max_nodes
nu_min_interp = nu_min_nodes

print(f"   Сетка анализа: {len(t_samples)} точек")

# 1. Маска экстрасистолы
if ectopic_mask is not None:
    ectopic_mask_samples = np.zeros_like(t_samples, dtype=bool)
    for i, t_curr in enumerate(t_samples):
        idx = np.argmin(np.abs(t - t_curr))
        ectopic_mask_samples[i] = ectopic_mask[idx]
else:
    ectopic_mask_samples = None

# 2. Маска паузы
if pause_mask is not None:
    pause_mask_samples = np.zeros_like(t_samples, dtype=bool)
    for i, t_curr in enumerate(t_samples):
        idx = np.argmin(np.abs(t - t_curr))
        pause_mask_samples[i] = pause_mask[idx]
    print(f"Пауз на сетке анализа: {np.sum(pause_mask_samples)}")
else:
    pause_mask_samples = None


# ============================================================================
# ФУНКЦИЯ ПОИСКА ПИКА
# ============================================================================
def find_peak_with_conditions(Power, nu, t, t_curr, nu_min, nu_max,
                              is_ectopic=False, is_pause=False):
    """
    - Экстрасистола: пик на ВЕРХНЕЙ границе
    - Пауза: пик на НИЖНЕЙ границе
    - Норма: поиск максимума внутри
    """
    # 1. Интерполяция спектра
    if t_curr <= t[0]:
        t_idx_left, t_idx_right, alpha = 0, 1, 0.0
    elif t_curr >= t[-1]:
        t_idx_left, t_idx_right, alpha = len(t) - 2, len(t) - 1, 1.0
    else:
        t_idx_right = np.searchsorted(t, t_curr)
        t_idx_left = t_idx_right - 1
        alpha = (t_curr - t[t_idx_left]) / (t[t_idx_right] - t[t_idx_left])

    p_interp = (1 - alpha) * Power[:, t_idx_left] + alpha * Power[:, t_idx_right]

    if is_ectopic:
        return nu_max

    if is_pause:
        return nu_min

    nu_mask = (nu >= nu_min) & (nu <= nu_max)
    nu_valid = nu[nu_mask]
    p_valid = p_interp[nu_mask]

    if len(nu_valid) < 2: return (nu_min + nu_max) / 2.0

    p_max, p_min = p_valid.max(), p_valid.min()
    p_range = p_max - p_min

    if p_range < 1e-12:
        return np.clip(np.sum(nu_valid * p_valid) / (p_valid.sum() or 1), nu_min, nu_max)

    pk_idx = np.argmax(p_valid)
    nu_peak = nu_valid[pk_idx]

    # Параболическая интерполяция
    if 0 < pk_idx < len(p_valid) - 1:
        try:
            a, b, _ = np.polyfit(nu_valid[pk_idx - 1:pk_idx + 2], p_valid[pk_idx - 1:pk_idx + 2], 2)
            if a < -1e-6:
                nu_opt = -b / (2 * a)
                if nu_min <= nu_opt <= nu_max: nu_peak = nu_opt
        except:
            pass

    return np.clip(nu_peak, nu_min, nu_max)


# ------------------------------------------------------------
# 6. ЗАПУСК ПОИСКА
# ------------------------------------------------------------
nu_peaks_raw = np.zeros_like(t_samples)

for i, t_curr in enumerate(t_samples):
    is_e = ectopic_mask_samples[i] if ectopic_mask_samples is not None else False
    is_p = pause_mask_samples[i] if pause_mask_samples is not None else False

    if is_e: is_p = False

    nu_peaks_raw[i] = find_peak_with_conditions(
        Power, nu, t, t_curr,
        nu_min_interp[i], nu_max_interp[i],
        is_ectopic=is_e,
        is_pause=is_p
    )

print(f" Диапазон: [{nu_peaks_raw.min():.3f}, {nu_peaks_raw.max():.3f}] Гц")
on_upper = np.sum(np.isclose(nu_peaks_raw, nu_max_interp, atol=1e-3))
on_lower = np.sum(np.isclose(nu_peaks_raw, nu_min_interp, atol=1e-3))
print(f"   Пики на верхней границе: {on_upper} ({on_upper / len(nu_peaks_raw) * 100:.1f}%)")
print(f"   Пики на нижней границе: {on_lower} ({on_lower / len(nu_peaks_raw) * 100:.1f}%)")

# ------------------------------------------------------------
# 7. СГЛАЖИВАНИЕ ТРАЕКТОРИИ
# ------------------------------------------------------------
nu_peaks = np.full_like(t_samples, nu_peaks_raw[0], dtype=float)

#tau_values = []
tau = 0.1
for n in range(1, len(t_samples)):
    delta = nu_peaks_raw[n] - nu_peaks_raw[n - 1]
    tc = (t_samples[n] + t_samples[n - 1]) / 2

    dt_local = t_samples[n] - t_samples[n - 1]
    nu_peaks += delta / (1 + np.exp(-(t_samples - tc) / tau))

nu_peaks = np.clip(nu_peaks, nu_min_interp, nu_max_interp)

# Диагностика
dt_mean = np.mean(np.diff(t_samples))
#tau_mean = np.mean(tau_values)
#tau_min = np.min(tau_values)
#tau_max = np.max(tau_values)

print(f"Средний шаг времени: {dt_mean:.6f} с")
#print(f"   Характеристическое время τ: среднее={tau_mean:.6f} с, диапазон=[{tau_min:.6f}, {tau_max:.6f}] с")
print(f"Диапазон сырой: [{nu_peaks_raw.min():.4f}, {nu_peaks_raw.max():.4f}] Гц")
print(f"Диапазон сглаженной: [{nu_peaks.min():.4f}, {nu_peaks.max():.4f}] Гц")

# ------------------------------------------------------------
# 8. ВИЗУАЛИЗАЦИЯ
# ------------------------------------------------------------
print(f"\n8. Построение графика локальной частоты...")

fig, ax = plt.subplots(1, 1, figsize=(14, 6))

# Заливка частотного коридора
ax.fill_between(t_samples, nu_min_interp, nu_max_interp,
                color='lightblue', alpha=0.3, label='Частотный коридор')

# Границы коридора
ax.plot(t_upper, nu_max_nodes, 'b--', linewidth=2.0, label='Верхняя граница $F^{>}(t)$')
ax.plot(t_lower, nu_min_nodes, 'r--', linewidth=2.0, label='Нижняя граница $F^{<}(t)$')

# Сглаженная локальная частота
ax.plot(t_samples, nu_peaks, 'k-', linewidth=3.0, label='Локальная частота $\\nu_{peak}(t)$')

ax.set_xlabel('Время, с', fontsize=13, fontweight='bold')
ax.set_ylabel('Частота, Гц', fontsize=13, fontweight='bold')
ax.set_title('Локальная частота (пик на границе ТОЛЬКО при экстрасистоле)',
             fontsize=14, fontweight='bold', pad=15)
ax.legend(loc='upper right', fontsize=11, framealpha=0.95)
ax.grid(True, alpha=0.4, linestyle='--', linewidth=0.8)
ax.set_xlim(t_samples.min(), t_samples.max())
ax.set_ylim(nu_min_interp.min() * 0.95, nu_max_interp.max() * 1.05)

plt.tight_layout()
plt.savefig('local_frequency_105.png', dpi=300, bbox_inches='tight')
plt.show()

# ------------------------------------------------------------
# 9. СОХРАНЕНИЕ РЕЗУЛЬТАТОВ
# ------------------------------------------------------------
print(f"\n9. Сохранение результатов...")

dt_save = 0.1
t_save = np.arange(t_samples.min(), t_samples.max() + dt_save, dt_save)

nu_peaks_save = np.interp(t_save, t_samples, nu_peaks)

output_file = 'peaks_105.txt'
np.savetxt(output_file, np.column_stack([t_save, nu_peaks_save]),
           fmt='%.6f', delimiter='\t', header='Time(s)\tNu_peak(Hz)', comments='')

print(f"Траектория сохранена: {output_file}")
print(f"Количество точек: {len(nu_peaks_save)} (было {len(nu_peaks)})")
print(f"Шаг времени: {dt_save:.3f} с (было {np.mean(np.diff(t_samples)):.4f} с)")
print(f"Диапазон времени: [{t_save.min():.2f}, {t_save.max():.2f}] с")
print(f"Диапазон частот: [{nu_peaks_save.min():.4f}, {nu_peaks_save.max():.4f}] Гц")