"""Streamlit interface for PsychAD Supplementary Table 6 gene-set testing."""
from io import BytesIO
import os
from pathlib import Path
import tempfile
os.environ.setdefault('MPLCONFIGDIR', str(Path(tempfile.gettempdir()) / 'psychad_matplotlib'))
import matplotlib.pyplot as plt
import numpy as np
import streamlit as st
from analysis import LIFESPAN_PATH, LIFESPAN_SCORES, AGING_PATH, AGING_GROUPS, AGING_LABELS, DATA_PATH, PHENOTYPES, LABELS, SCORES, load_data, parse_genes, run_analysis

st.set_page_config(page_title='PsychAD · Gene-set AUROCs', page_icon='🧠', layout='wide')
dataset = st.sidebar.selectbox('Dataset', ['AD phenotypes', 'Aging — whole lifespan', 'Aging — within age bins'])
is_aging = dataset != 'AD phenotypes'
is_lifespan = dataset == 'Aging — whole lifespan'
source_path = LIFESPAN_PATH if is_lifespan else AGING_PATH if is_aging else DATA_PATH
phenotypes = [f'AveExpr_z_{g}' for g in AGING_GROUPS] + ['Coarse_AveExpr_slope', 'Whole_lifespan'] if is_lifespan else AGING_GROUPS if is_aging else PHENOTYPES
labels = {**dict(zip([f'AveExpr_z_{g}' for g in AGING_GROUPS], ['Development\nAveExpr z-score', 'Young adulthood\nAveExpr z-score', 'Middle adulthood\nAveExpr z-score', 'Late adulthood\nAveExpr z-score'])), 'Whole_lifespan': 'Whole lifespan\nassociation (F)', 'Coarse_AveExpr_slope': 'Coarse AveExpr\nslope (4 bins)'} if is_lifespan else AGING_LABELS if is_aging else LABELS
st.title('PsychAD · Gene-set AUROCs · Cell type × ' + ('whole lifespan' if is_lifespan else 'age period' if is_aging else 'AD phenotype'))
st.write('Test a gene set against cell-type-specific differential expression across ' + ('the whole lifespan.' if is_lifespan else 'age periods.' if is_aging else 'Alzheimer’s disease phenotypes.'))
if is_lifespan:
    st.markdown('<p style="font-size:18px;color:#000000;line-height:1.6;">Differential expression statistics from '
                '<a href="https://www.nature.com/articles/s41586-026-10271-7">Lifespan single-cell transcriptomic atlas of the human prefrontal cortex</a> '
                '(Hui Yang et al., <em>Nature</em>) · Supplementary Data 7. '
                'Donor-level pseudobulk Dreamlet model across 284 neurotypical donors aged 0–97 years and 26 cell types: '
                'expression ~ poly(log2(age + 1), degree = 2) + covariates. '
                ''
                'Both polynomial terms are fitted across the whole lifespan, without age bins. ENSG IDs are excluded.</p>', unsafe_allow_html=True)
elif is_aging:
    st.markdown(
        '<p style="font-size: 18px; color: #000000; line-height: 1.6;">'
        'Differential expression statistics from '
        '<a href="https://www.nature.com/articles/s41586-026-10271-7" style="color: #000000; text-decoration: underline;">'
        'Lifespan single-cell transcriptomic atlas of the human prefrontal cortex</a> '
        '(Hui Yang et al., <em>Nature</em>) · Supplementary Data 4 (supplied filtered file). '
        'Donor-level pseudobulk age-association analysis using Dreamlet in human dorsolateral prefrontal cortex. '
        ''
        'The cohort includes 284 neurotypical donors: development n = 53, young adulthood n = 54, '
        'middle adulthood n = 95 and late adulthood n = 82. The supplied data cover 26 cell types; '
        'EN_L5_ET was excluded in the study because of low nucleus counts.</p>', unsafe_allow_html=True,
    )
else:
    st.markdown(
        '<p style="font-size: 18px; color: #000000; line-height: 1.6;">'
        'Differential expression statistics from '
        '<a href="https://doi.org/10.1038/s41586-025-09573-z" target="_blank" rel="noopener noreferrer" '
        'style="color: #000000; text-decoration: underline;">'
        'Single-cell atlas of transcriptomic vulnerability across brain disorders</a>'
        ' (Donghoon Lee et al., <em>Nature</em>) · Supplementary Table 6. '
        'Donor-level pseudobulk differential expression using Dreamlet in human dorsolateral '
        'prefrontal cortex, across 27 cell types (<strong>n = 696 donors</strong>).</p>',
        unsafe_allow_html=True,
    )

@st.cache_data(show_spinner='Loading differential-expression data…')
def cached_data(path, mtime):
    return load_data(Path(path))


try:
    data = cached_data(str(source_path), source_path.stat().st_mtime_ns)
except (OSError, ValueError) as exc:
    st.error(str(exc))
    st.stop()

with st.sidebar:
    st.header('Input genes')
    st.caption('Paste human symbols separated by spaces, commas, semicolons or newlines. Results update automatically when inputs change.')
    target_text = st.text_area('Input gene set', value='TREM2\nC4B\nLILRA2\nNPNT\nMSR1\nSLC11A1\nRP11-203I2\nMYL3\nGYPC\nCD44\nSERPINA3\nCHIT1\nTNFAIP2\nPIRT\nPADI2', height=180)
    background_text = st.text_area('Background gene set (optional)', height=140, help='Leave blank to use all genes in the selected dataset. A custom background is the eligible gene universe, including target genes.')
    if is_lifespan:
        score = 'Overall lifespan association (F-statistic)'
        st.caption('Whole lifespan: Four age-bin AveExpr z-scores, coarse slope and F-statistic are tested side by side.')
    else:
        score = st.selectbox('Gene ranking score', list(SCORES), key='standard_score')
    st.caption('Default example: 15 positively tau-correlated genes from Figure 2b of [Huang W, Bartosch AM, Xiao H, et al. (2021), Nature Communications 12:5659](https://doi.org/10.1038/s41467-021-25902-y), “An immune response characterizes early Alzheimer’s disease pathology and subjective cognitive impairment in hydrocephalus biopsies.”')

targets = parse_genes(target_text)
background = parse_genes(background_text) or None
if not targets:
    st.info('Enter a target gene set to see the results.')
    st.stop()
try:
    with st.spinner('Testing each cell type and phenotype…'):
        results = run_analysis(data, targets, background, score)
except ValueError as exc:
    st.error(str(exc))
    st.stop()
# Positive scores indicate AUROC > 0.5; negative scores indicate AUROC < 0.5.
# Clip numerical zero p-values to the smallest positive float for finite log scores.
results['Signed -log10(p)'] = np.sign(results['AUROC'] - 0.5) * -np.log10(
    results['p-value'].clip(lower=np.nextafter(0.0, 1.0))
)
results = results.sort_values(
    ['Signed -log10(p)', 'AUROC'], ascending=[False, False],
    kind='stable', na_position='last',
).reset_index(drop=True)
# Display alias only; source data and analysis identifiers remain unchanged.
results['Cell type'] = results['Cell type'].replace({'Adaptive': 'Adaptive Immune'})
columns = [column for column in results.columns if column not in ('FDR (BH)', 'Significance')]
status_index = columns.index('Status')
columns[status_index:status_index] = ['FDR (BH)', 'Significance']
results = results[columns]
if is_lifespan:
    results['Score'] = results['Phenotype'].map({**{f'AveExpr_z_{g}': f'AveExpr z-score: {g}' for g in AGING_GROUPS}, 'Whole_lifespan': 'F-statistic', 'Coarse_AveExpr_slope': 'AveExpr OLS slope across bins 1–4'})
run_score = 'Age-bin AveExpr z-scores + coarse slope + F-statistic' if is_lifespan else score
st.caption(f'{run_score} · {len(targets):,} unique target genes · background: ' + ('all source genes' if background is None else f'{len(background):,} submitted genes'))
missing = sorted(set(targets) - set(data.ID))
outside = sorted(set(targets) - set(background)) if background is not None else []
if missing:
    st.warning('Target genes absent from the selected dataset: ' + ', '.join(missing))
if outside:
    st.warning('Target genes excluded by the custom background: ' + ', '.join(outside))
if background is not None:
    absent_background = set(background) - set(data.ID)
    if absent_background:
        st.caption(f'{len(absent_background):,} background genes are absent from the selected dataset.')


view = results.copy()
view['Significance'] = view['Significance'].replace('', 'NS')
# Full subclass order shown top-to-bottom in the paper's Figure 6a.
paper_order = [
    'IN_LAMP5_RELN', 'IN_LAMP5_LHX6', 'IN_ADARB2', 'IN_VIP',
    'IN_PVALB', 'IN_PVALB_CHC', 'IN_SST',
    'EN_L2_3_IT', 'EN_L6_IT_1', 'EN_L3_5_IT_2', 'EN_L3_5_IT_3',
    'EN_L3_5_IT_1', 'EN_L6_IT_2', 'EN_L6_CT', 'EN_L6B', 'EN_L5_ET',
    'EN_L5_6_NP', 'Astro', 'OPC', 'Oligo', 'Micro', 'PVM',
    'Adaptive Immune', 'VLMC', 'SMC', 'PC', 'Endo',
]
available_cells = set(results['Cell type'])
cells = [cell for cell in paper_order if cell in available_cells]
cells += sorted(available_cells - set(paper_order))
indexed = results.set_index(['Cell type', 'Phenotype'])
matrix = np.array([[indexed.loc[(c, p), 'AUROC'] for p in phenotypes] for c in cells])
fig, ax = plt.subplots(figsize=(16 if is_lifespan else 11, max(3, len(cells) * .53 + 1.3)))
cmap = plt.get_cmap('RdBu_r').copy()
cmap.set_bad('#e5e7eb')
im = ax.imshow(matrix, vmin=0, vmax=1, cmap=cmap, aspect='auto')
ax.set_xticks(range(len(phenotypes)), [labels[p] for p in phenotypes])
ax.xaxis.tick_top()
ax.set_yticks(range(len(cells)), cells)
ax.tick_params(length=0, pad=8)
for i, cell in enumerate(cells):
    for j, phenotype in enumerate(phenotypes):
        row = indexed.loc[(cell, phenotype)]
        auc = row['AUROC']
        counts = f"n={row['Target genes (n)']:,} / {row['Comparison genes (n)']:,}"
        label = f"{auc:.3f}\n{counts}\np={row['p-value']:.2g}" if np.isfinite(auc) else f"Not testable\n{counts}"
        ax.text(j, i, label, ha='center', va='center', fontsize=8, color='white' if np.isfinite(auc) and abs(auc-.5) > .32 else '#111827')
        if np.isfinite(auc) and row['Significance']:
            ax.text(j + .34, i - .24, row['Significance'], ha='center', va='center',
                    fontsize=16, fontweight='bold',
                    color='white' if abs(auc-.5) > .32 else '#111827')
ax.set_xticks(np.arange(-.5, len(phenotypes), 1), minor=True)
ax.set_yticks(np.arange(-.5, len(cells), 1), minor=True)
ax.grid(which='minor', color='white', linewidth=1.5)
ax.tick_params(which='minor', length=0)
for spine in ax.spines.values():
    spine.set_visible(False)
fig.colorbar(im, ax=ax, fraction=.025, pad=.025, label='AUROC')
fig.tight_layout()
st.pyplot(fig, use_container_width=True)
interpretation = (
    'the input gene set is shifted toward <strong>higher expression with increasing age within the indicated age period</strong> '
    if is_aging else
    'the input gene set is shifted toward <strong>higher expression in Alzheimer’s disease cases</strong> (versus controls), or higher expression with <strong>increasing CERAD, Braak or dementia severity</strong> '
)
if is_lifespan:
    st.markdown('<p style="font-size:18px;color:#000000;line-height:1.6;"><strong>How to read the heatmap:</strong> The first four columns rank each gene–cell-type pair’s AveExpr z-scores across the four age bins: (bin AveExpr − four-bin mean) / four-bin sample standard deviation. Red indicates that target genes have relatively higher standardized expression in that age bin than background genes; it does not indicate an age slope within the bin. The F-statistic column tests overall lifespan association strength, regardless of direction. The coarse slope column fits a straight line through each gene–cell-type pair’s four Supplementary Data 4 AveExpr values, against equally spaced age-bin positions 1, 2, 3, 4 (development, young, middle, late adulthood). Positive slopes indicate increasing average expression; negative slopes indicate decreasing expression. In each column, red means target genes rank higher than background genes. A red slope AUROC indicates relatively more positive slopes, not necessarily that all target genes increase.</p>'
                '<p style="font-size:18px;color:#000000;line-height:1.6;">The slope is a descriptive change per bin step, not per year or a donor-level covariate-adjusted age effect. All four finite AveExpr values are required; missing bins are not imputed. Constant four-bin profiles have undefined z-scores and are excluded from z-score tests, but retain a zero slope. Counts can therefore differ between columns. No ENSG IDs or significance filtering are used.</p>'
                '<p style="font-size:18px;color:#000000;line-height:1.6;">Color: AUROC (midpoint 0.5). Text: target n / comparison n and raw p-value. Stars use Benjamini–Hochberg correction across all valid tests in all six columns: <strong>***</strong> q &lt; 0.001; <strong>**</strong> q &lt; 0.01; <strong>*</strong> q &lt; 0.05. Rows retain the AD paper’s Figure 6a cell-subclass order.</p>', unsafe_allow_html=True)
else:
    st.markdown(
        f'<p style="font-size: 18px; color: #000000; line-height: 1.6;"><strong>How to read the heatmap:</strong> With signed scores, <strong>AUROC &gt; 0.5 (red)</strong> means {interpretation}relative to the background genes. <strong>AUROC &lt; 0.5 (blue)</strong> indicates the opposite direction; <strong>0.5</strong> indicates no directional shift. This describes the gene set overall, not necessarily every gene.</p>'
        '<p style="font-size: 18px; color: #000000; line-height: 1.6;">Color: AUROC (0–1; midpoint 0.5). Text: target n / comparison n, raw p-value, and FDR stars. FDR includes all valid tests in the completed run.</p>'
        '<p style="font-size: 18px; color: #000000; line-height: 1.6;">Significance (Benjamini–Hochberg adjusted p-value, q): <strong>***</strong> q &lt; 0.001; <strong>**</strong> q &lt; 0.01; <strong>*</strong> q &lt; 0.05; no stars: q ≥ 0.05 or not testable.</p>'
        '<p style="font-size: 18px; color: #000000; line-height: 1.6;">Rows use the AD paper’s Figure 6a ordering for cell types present in the selected dataset: inhibitory neurons, excitatory neurons, glia, immune cells, then vascular cells. Adaptive Immune is labelled Adaptive in the source table.</p>'
        '<p style="font-size: 18px; color: #000000; line-height: 1.6;">If you select Absolute t-statistic, red instead means stronger associations regardless of whether expression increases or decreases.</p>'
        , unsafe_allow_html=True,
    )
image = BytesIO()
fig.savefig(image, format='png', dpi=180, bbox_inches='tight')
st.download_button('Download heatmap PNG', image.getvalue(), f'psychad_{dataset.lower().replace(" ", "_")}_auroc_heatmap.png', 'image/png')
plt.close(fig)

st.subheader('AUROC results')
st.caption('Sorted by signed −log10(p), descending: sign = sign(AUROC − 0.5). Strong positive shifts appear first; strong negative shifts appear last. Untestable results appear at the bottom.')
st.dataframe(view, hide_index=True, use_container_width=True, column_config={'AUROC': st.column_config.NumberColumn(format='%.4f'), 'p-value': st.column_config.NumberColumn(format='%.3g'), 'FDR (BH)': st.column_config.NumberColumn(format='%.3g')})
st.download_button('Download all results CSV', view.assign(Dataset=dataset, **({} if is_lifespan else {'Score': run_score})).to_csv(index=False), f'psychad_{dataset.lower().replace(" ", "_")}_auroc_results.csv', 'text/csv')
with st.expander('Gene coverage by cell type and phenotype'):
    st.dataframe(results[['Cell type', 'Phenotype', 'Target genes (n)', 'Comparison genes (n)', 'Universe genes (n)', 'Available genes (n)', 'Target genes used', 'Target genes not used']], hide_index=True, use_container_width=True)
