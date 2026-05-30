#модуль 1 - фильтрация ЭКГ

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import wfdb
from scipy.interpolate import interp1d
from datetime import datetime

# ============================================================================
# КОНФИГУРАЦИЯ
# ============================================================================

RECORD_NUM = 105
DATABASE_DIR = 'mit-bih-arrhythmia-database-1.0.0'

# окно обратботки
T_START = 0
T_END = 1280

# окно визуализации
T_START_TAH = 0
T_END_TAH = 1280

# Артефакты, шум, неопределённые или нефизиологические маркеры
OUTLIER_BEAT_SYMBOLS = ['!', '"', '+', '|', '~', '*', 'Q', '?', '[', ']', '#']
# Валидные сердечные сокращения
KEEP_BEAT_SYMBOLS = ['N', 'L', 'R', 'B', 'A', 'a', 'J', 'S', 'V', 'E', 'F', 'f', 'e', 'j', 'n', '/']


# ============================================================================
# ФУНКЦИИ
# ============================================================================

def load_and_process(record_num, database_dir, t_start, t_end):
    print("=" * 70)
    print("ЗАГРУЗКА И ОБРАБОТКА ДАННЫХ")
    print("=" * 70)
    print(f"Запись: {record_num} | Окно: [{t_start}, {t_end}] сек")
    print(f"Дата: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)

    # 1. Загрузка данных (современный API wfdb)
    print("\n[1/4] Загрузка данных...")
    record = wfdb.rdrecord(f'{database_dir}/{record_num}')
    ecg_signal = record.p_signal[:, 0]
    fs = record.fs

    annotation = wfdb.rdann(f'{database_dir}/{record_num}', 'atr')
    samples = np.array(annotation.sample)
    symbols = np.array(annotation.symbol)
    print(f"fs={fs} Гц | Аннотаций: {len(symbols)}")

    # 2. Фильтрация по временному окну
    print("\n[2/4] Фильтрация по окну...")
    sample_start = int(t_start * fs)
    sample_end = int(t_end * fs)

    time_mask = (samples >= sample_start) & (samples <= sample_end)
    seg_samples = samples[time_mask] - sample_start
    seg_symbols = symbols[time_mask]
    seg_times = seg_samples / fs

    unique_symbols, counts = np.unique(seg_symbols, return_counts=True)
    symbol_stats = dict(zip(unique_symbols, counts))
    print(f"Статистика аннотаций: {symbol_stats}")

    # 3. RR-интервалы
    print("\n[3/4] Вычисление RR-интервалов...")
    rr_intervals = np.diff(seg_samples) / fs * 1000  # мс
    rr_times = seg_times[1:]
    rr_symbols_pre = seg_symbols[:-1]
    rr_symbols_post = seg_symbols[1:]

    valid_mask = rr_intervals > 0
    rr_intervals = rr_intervals[valid_mask]
    rr_times = rr_times[valid_mask]
    rr_symbols_pre = rr_symbols_pre[valid_mask]
    rr_symbols_post = rr_symbols_post[valid_mask]

    print(f"Всего вычислено интервалов: {len(rr_intervals)}")

    # 4. Фильтрация выбросов
    print("\n[4/4] Детекция выбросов по меткам аннотаций...")

    # Интервал считается артефактом, если хотя бы один из двух ударов имеет метку из OUTLIER
    final_outlier_mask = np.isin(rr_symbols_pre, OUTLIER_BEAT_SYMBOLS) | \
                         np.isin(rr_symbols_post, OUTLIER_BEAT_SYMBOLS)

    print(f"Статистика выбросов:")
    print(f"Интервалы с артефактами: {np.sum(final_outlier_mask)}")

    # интерполяция
    rr_final = rr_intervals.copy()
    if np.any(final_outlier_mask):
        good_mask = ~final_outlier_mask
        rr_good = rr_intervals[good_mask]

        if np.sum(good_mask) >= 4:
            t_good = rr_times[good_mask]
            f = interp1d(t_good, rr_good, kind='linear', fill_value='extrapolate', bounds_error=False)
            rr_final[final_outlier_mask] = f(rr_times[final_outlier_mask])
            print(f" Интерполировано выбросов: {np.sum(final_outlier_mask)}")
        else:
            median_val = np.median(rr_good)
            rr_final[final_outlier_mask] = median_val
            print(f"Использована медиана для {np.sum(final_outlier_mask)} точек")
    else:
        print(f"Выбросы не обнаружены")

    # Защита от экстраполяционных артефактов (физиологически допустимый диапазон 300–2000 мс)
    rr_final = np.clip(rr_final, 300, 2000)

    pvc_mask = (rr_symbols_pre == 'V') | (rr_symbols_post == 'V')
    pvc_count = np.sum(pvc_mask)

    print(f"Обработано интервалов: {len(rr_final)}")
    print(f"Сохранено желудочковых экстрасистол (PVC): {pvc_count}")

    return {
        'rr_final': rr_final, 'rr_times': rr_times,
        'rr_symbols_pre': rr_symbols_pre, 'rr_symbols_post': rr_symbols_post,
        'pvc_mask': pvc_mask, 'outlier_mask': final_outlier_mask,
        'pvc_count': pvc_count, 't_start': t_start, 't_end': t_end,
        'record_num': record_num, 'seg_symbols': seg_symbols, 'seg_times': seg_times,
        'outlier_count': np.sum(final_outlier_mask),
        'outlier_indices': np.where(final_outlier_mask)[0]
    }


def export_for_analysis(data, output_file=None):
    if output_file is None:
        output_file = f'rr_final_{data["record_num"]}_{int(data["t_end"])}sec.csv'
    df = pd.DataFrame({'t_sec': data['rr_times'], 'RR_ms': data['rr_final']})
    df.to_csv(output_file, index=False, float_format='%.3f')
    print(f"Экспорт: {output_file}")
    return df


def plot_tachogram(data, t_start_tah=None, t_end_tah=None, output_path=None):
    print("\n" + "=" * 70)
    print("ВИЗУАЛИЗАЦИЯ ТАХОГРАММЫ")
    print("=" * 70)

    if t_start_tah is None: t_start_tah = data['t_start']
    if t_end_tah is None: t_end_tah = data['t_end']
    print(f"Окно визуализации: [{t_start_tah}, {t_end_tah}] секунд")

    viz_mask = (data['rr_times'] >= t_start_tah) & (data['rr_times'] <= t_end_tah)
    rr_viz = data['rr_final'][viz_mask]
    rr_times_viz = data['rr_times'][viz_mask]

    pvc_indices = np.where(data['seg_symbols'] == 'V')[0]
    pvc_times = data['seg_times'][pvc_indices]
    pvc_mask_viz = (pvc_times >= t_start_tah) & (pvc_times <= t_end_tah)
    pvc_times_viz = pvc_times[pvc_mask_viz]
    pvc_count_viz = len(pvc_times_viz)

    fig, ax = plt.subplots(figsize=(18, 7))
    ax.plot(rr_times_viz, rr_viz, 'k-', linewidth=1.5, label='RR-интервалы', alpha=0.95)
    ax.fill_between(rr_times_viz, rr_viz, 0, color='black', alpha=0.3, label='Заливка')

    if pvc_count_viz > 0:
        for t in pvc_times_viz:
            ax.axvline(x=t, color='red', linestyle='--', linewidth=1.5, alpha=0.85, zorder=5)

    ax.set_xlim(t_start_tah, t_end_tah)
    ax.set_xlabel('Время, секунды', fontsize=14, fontweight='bold')
    ax.set_ylabel('RR, мс', fontsize=14, fontweight='bold')

    legend_lines = [plt.Line2D([0], [0], color='black', linewidth=1.5)]
    legend_labels = ['RR-интервалы']
    legend_lines.append(plt.Line2D([0], [0], color='none'))
    legend_labels.append(f'Обработка: [{data["t_start"]}, {data["t_end"]}] сек')

    if pvc_count_viz > 0:
        pvc_str = ', '.join([f'{t:.1f}' for t in pvc_times_viz[:5]])
        if pvc_count_viz > 5: pvc_str += '...'
        legend_lines.append(plt.Line2D([0], [0], color='red', linestyle='--', linewidth=1.5))
        legend_labels.append(f'PVC ({pvc_count_viz}): {pvc_str} сек')

    ax.legend(legend_lines, legend_labels, loc='upper right', fontsize=10, framealpha=0.95)
    ax.set_title(f'Тахограмма RR-интервалов — Запись {data["record_num"]}', fontsize=15, fontweight='bold', pad=15)
    ax.grid(True, alpha=0.45, linestyle='--', linewidth=0.8)
    ax.set_axisbelow(True)

    # 🛠 ИСПРАВЛЕНИЕ: сохранение графика перед показом
    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"График сохранен: {output_path}")
    plt.show()


def export_detailed_report(data, output_file=None):
    if output_file is None:
        output_file = f'report_{data["record_num"]}_{int(data["t_end"])}sec.txt'
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write("=" * 70 + "\nДЕТАЛЬНЫЙ ОТЧЕТ ОБ ОБРАБОТКЕ ДАННЫХ\n" + "=" * 70 + "\n\n")
        f.write(f"Запись: {data['record_num']}\nВременное окно: [{data['t_start']}, {data['t_end']}] сек\n")
        f.write(f"Дата: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write(f"Всего интервалов: {len(data['rr_final'])}\nИнтерполировано артефактов: {data['outlier_count']}\n")
        f.write(f"Сохранено PVC: {data['pvc_count']}\n\n")
        f.write(f"Средний RR: {np.mean(data['rr_final']):.2f} ± {np.std(data['rr_final']):.2f} мс\n")
        f.write(
            f"Медиана: {np.median(data['rr_final']):.2f} мс | ЧСС: {60000 / np.mean(data['rr_final']):.2f} уд/мин\n\n")

        pvc_indices = np.where(data['seg_symbols'] == 'V')[0]
        pvc_times = data['seg_times'][pvc_indices]
        if len(pvc_times) > 0:
            for i, t in enumerate(pvc_times):
                f.write(f"PVC #{i + 1}: {t:.3f} с от начала записи ({t - data['t_start']:.3f} с от окна)\n")
        else:
            f.write("PVC не обнаружены\n")
    print(f"Отчет сохранен: {output_file}")
    return output_file


# ============================================================================
# ОСНОВНОЙ ПЛАЙПЛАЙН
# ============================================================================
if __name__ == '__main__':
    try:
        print("\n" + "=" * 70 + "\nЗАПУСК ПЛАЙПЛАЙНА ОБРАБОТКИ ЭКГ\n" + "=" * 70 + "\n")
        data = load_and_process(RECORD_NUM, DATABASE_DIR, T_START, T_END)

        stats = {
            'mean_rr': np.mean(data['rr_final']), 'std_rr': np.std(data['rr_final']),
            'min_rr': np.min(data['rr_final']), 'max_rr': np.max(data['rr_final']),
            'median_rr': np.median(data['rr_final']),
            'heart_rate': 60000 / np.mean(data['rr_final']),
            'total_beats': len(data['rr_final']), 'pvc_count': data['pvc_count'],
            'outlier_count': data['outlier_count']
        }
        print("\n" + "=" * 70 + "\nФИНАЛЬНАЯ СТАТИСТИКА\n" + "=" * 70)
        for k, v in stats.items():
            print(f"{k.replace('_', ' ').title()}: {v:.2f}" if isinstance(v,
                                                                          float) else f"{k.replace('_', ' ').title()}: {v}")

        print("\n" + "=" * 70 + "\nЭКСПОРТ ДАННЫХ\n" + "=" * 70)
        csv_file = export_for_analysis(data)
        report_file = export_detailed_report(data)

        print("\n" + "=" * 70 + "\nВИЗУАЛИЗАЦИЯ ТАХОГРАММЫ\n" + "=" * 70)
        plot_tachogram(data, t_start_tah=T_START_TAH, t_end_tah=T_END_TAH,
                       output_path=f'tachogram_{RECORD_NUM}_{int(T_START_TAH)}-{int(T_END_TAH)}sec.png')

        print("\n" + "=" * 70 + "\n✅ ПЛАЙПЛАЙН ЗАВЕРШЕН УСПЕШНО!\n" + "=" * 70)
    except Exception as e:
        print(f"\nКРИТИЧЕСКАЯ ОШИБКА: {e}")
        import traceback

        traceback.print_exc()