"""Gene-set AUROCs over published PsychAD differential-expression statistics."""
from pathlib import Path
import re
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

DATA_PATH = Path(__file__).parent / 'data/Donghoon_Lee...Panos_Roussos/PsychAD_SupplementaryTable6.csv'
PHENOTYPES = ['dx_AD', 'CERAD', 'Braak', 'Dementia']
LABELS = {'dx_AD': 'AD diagnosis', 'CERAD': 'CERAD', 'Braak': 'Braak', 'Dementia': 'Dementia'}
SCORES = {'Signed t-statistic': ('statistic', False), 'Absolute t-statistic': ('statistic', True), 'Signed effect estimate': ('estimate', False)}


def parse_genes(text):
    return list(dict.fromkeys(x for x in re.split(r'[\s,;]+', text.strip()) if x))


def load_data(path=DATA_PATH):
    data = pd.read_csv(path, usecols=['coef', 'assay', 'ID', 'statistic', 'estimate'])
    if data.duplicated(['coef', 'assay', 'ID']).any():
        raise ValueError('Duplicate genes within a phenotype/cell type in the source table.')
    return data


def adjust_bh(p):
    p = np.asarray(p, dtype=float)
    if not len(p):
        return p
    order = np.argsort(p)
    q = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    result = np.empty(len(p))
    result[order] = np.minimum(q, 1)
    return result


def stars(q):
    return '***' if q < .001 else '**' if q < .01 else '*' if q < .05 else ''


def run_analysis(data, targets, background=None, score='Signed t-statistic'):
    targets = set(targets)
    if not targets:
        raise ValueError('Enter at least one target gene.')
    universe = set(data.ID) if background is None else set(background)
    if not universe:
        raise ValueError('The custom background is empty.')
    if not targets.intersection(universe).intersection(data.ID):
        raise ValueError('No target genes overlap both the source table and the background.')
    column, absolute = SCORES[score]
    rows = []
    for (cell, phenotype), group in data.groupby(['assay', 'coef'], sort=True):
        available = group[np.isfinite(group[column])]
        used = available[available.ID.isin(universe)]
        mask = used.ID.isin(targets)
        x, y = used.loc[mask, column].to_numpy(), used.loc[~mask, column].to_numpy()
        if absolute:
            x, y = abs(x), abs(y)
        n, m = len(x), len(y)
        auc = p = np.nan
        status = 'OK'
        if not n or not m:
            status = 'No target genes' if not n else 'No comparison genes'
        else:
            # Average ranks handle ties; asymptotic test applies tie and continuity corrections.
            test = mannwhitneyu(x, y, alternative='two-sided', method='asymptotic', use_continuity=True)
            auc, p = test.statistic / (n * m), test.pvalue
            if np.all(x == x[0]) and np.all(y == x[0]):
                auc, p = 0.5, 1.0
        rows.append({'Cell type': cell, 'Phenotype': phenotype, 'AUROC': auc,
                     'Target genes (n)': n, 'Comparison genes (n)': m,
                     'Universe genes (n)': n + m, 'Available genes (n)': len(available),
                     'p-value': p, 'Status': status,
                     'Target genes used': ';'.join(sorted(used.loc[mask, 'ID'])),
                     'Target genes not used': ';'.join(sorted(targets - set(used.loc[mask, 'ID'])))})
    result = pd.DataFrame(rows)
    result['FDR (BH)'] = np.nan
    valid = result['p-value'].notna()
    result.loc[valid, 'FDR (BH)'] = adjust_bh(result.loc[valid, 'p-value'])
    result['Significance'] = result['FDR (BH)'].map(stars)
    return result
