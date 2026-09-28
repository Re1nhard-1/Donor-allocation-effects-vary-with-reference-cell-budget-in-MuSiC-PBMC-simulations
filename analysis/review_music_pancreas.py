import argparse
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import io

from scientific_style import configure_style

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data/processed/music_pancreas'
OUT = ROOT / 'results/music_pancreas'


def digest(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def save(frame, name):
    frame.to_csv(OUT / name, index=False, lineterminator='\n')


def read_validate(mode):
    design = json.loads((DATA / 'design.json').read_text())
    donors, types = design['donors'], design['cell_types']
    for name, expected in design['files'].items():
        assert digest(DATA / name) == expected
    targets = pd.read_csv(DATA / 'targets.tsv', sep='\t')
    refs = pd.read_csv(DATA / 'references.tsv', sep='\t')
    cells = pd.read_csv(DATA / 'cells.tsv', sep='\t')
    counts = io.mmread(DATA / 'counts.mtx.gz').tocsc()
    bulk = io.mmread(DATA / 'targets.mtx.gz').tocsc()
    for index, target in targets.iterrows():
        ids = np.array(list(map(int, target.columns_R.split('|')))) - 1
        assert len(np.unique(ids)) == 100
        assert cells.iloc[ids].donor.eq(target.held_out).all()
        np.testing.assert_array_equal(counts[:, ids].sum(axis=1), bulk[:, index].toarray())
        for cell_type in types:
            assert np.isclose(cells.iloc[ids].cell_type.eq(cell_type).mean(), target[f'true_{cell_type}'])
    for ref in refs.itertuples():
        ids = np.array(list(map(int, ref.columns_R.split('|')))) - 1
        assert len(ids) == len(np.unique(ids)) == ref.budget * len(types)
        selected = cells.iloc[ids]
        assert not selected.donor.eq(ref.held_out).any()
        for donor in set(donors) - {ref.held_out}:
            quota = ref.budget // 3 if ref.dominant == 'balanced' else (ref.budget * 10 // 12 if donor == ref.dominant else ref.budget // 12)
            for cell_type in types:
                assert ((selected.donor == donor) & (selected.cell_type == cell_type)).sum() == quota
    paths = [OUT / donor / f'{mode}_predictions.csv' for donor in (donors[:1] if mode == 'pilot' else donors)]
    frame = pd.concat([pd.read_csv(path) for path in paths], ignore_index=True)
    expected_refs = refs[(refs.held_out == donors[0]) & (refs.block == 0)] if mode == 'pilot' else refs
    assert len(frame) == len(expected_refs) * 60 * 2
    assert set(frame.reference_id) == set(expected_refs.reference_id)
    assert not frame.duplicated(['reference_id', 'target_name', 'method']).any()
    joined = frame.merge(targets[['target_name', 'held_out'] + [f'true_{t}' for t in types]], on='target_name', suffixes=('', '_expected'), validate='many_to_one')
    assert joined.held_out.eq(joined.held_out_expected).all()
    expected_metadata = refs[['reference_id', 'held_out', 'block', 'budget', 'dominant']]
    metadata = frame[['reference_id', 'held_out', 'block', 'budget', 'dominant']].drop_duplicates()
    actual_metadata = metadata.merge(expected_metadata, on='reference_id', suffixes=('', '_expected'), validate='one_to_one')
    for field in ['held_out', 'block', 'budget', 'dominant']:
        assert actual_metadata[field].eq(actual_metadata[f'{field}_expected']).all()
    assert frame.groupby(['reference_id', 'method']).size().eq(60).all()
    true = frame[[f'true_{t}' for t in types]].to_numpy()
    pred = frame[[f'pred_{t}' for t in types]].to_numpy()
    np.testing.assert_allclose(true, joined[[f'true_{t}_expected' for t in types]], atol=1e-14, rtol=0)
    assert np.isfinite(pred).all() and pred.min() >= -1e-10
    np.testing.assert_allclose(pred.sum(axis=1), 1, atol=1e-10, rtol=0)
    np.testing.assert_allclose(frame.mae_pp, np.abs(pred - true).mean(axis=1) * 100, atol=1e-10, rtol=0)
    np.testing.assert_allclose(frame.rmse_pp, np.sqrt(np.square(pred - true).mean(axis=1)) * 100, atol=1e-10, rtol=0)
    checks = {'mode': mode, 'rows': len(frame), 'references': len(expected_refs), 'weighted_fits': int(frame.method.eq('weighted').sum()),
              'maxiter_weighted': int((frame.method.eq('weighted') & ~frame.convergence.str.startswith('Converge')).sum()),
              'nonfinite_variance_weighted': int((frame.method.eq('weighted') & ~frame.variance_finite).sum()),
              'effective_gene_range': [int(frame.n_features.min()), int(frame.n_features.max())],
              'weighted_prediction_maximum_by_type': {t: float(frame.loc[frame.method.eq('weighted'), f'pred_{t}'].max()) for t in types},
              'source_contract_and_score_checks': 'passed', 'design_sha256': digest(DATA / 'design.json'),
              'prediction_sha256': {p.relative_to(ROOT).as_posix(): digest(p) for p in paths},
              'review_script_sha256': digest(__file__)}
    (OUT / f'{mode}_validation.json').write_text(json.dumps(checks, indent=2) + '\n', encoding='utf-8', newline='\n')
    print(json.dumps(checks, indent=2))
    return frame, design


def summarize(frame, design):
    types = design['cell_types']
    for t in types:
        frame[f'ae_{t}_pp'] = 100 * abs(frame[f'pred_{t}'] - frame[f'true_{t}'])
    metrics = ['mae_pp', 'rmse_pp'] + [f'ae_{t}_pp' for t in types]
    keys = ['method', 'held_out', 'block', 'budget', 'dominant']
    references = frame.groupby(keys, as_index=False)[metrics].mean()
    paired = references.loc[references.dominant.ne('balanced')].merge(
        references.loc[references.dominant.eq('balanced')].drop(columns='dominant'),
        on=keys[:-1], suffixes=('_unequal', '_balanced'), validate='many_to_one')
    for metric in metrics:
        paired[f'delta_{metric}'] = paired[f'{metric}_unequal'] - paired[f'{metric}_balanced']
    delta = [f'delta_{metric}' for metric in metrics]
    donor = paired.groupby(['method', 'held_out', 'budget'], as_index=False)[delta].mean()
    arrangements = paired.groupby(['method', 'held_out', 'budget', 'dominant'], as_index=False)[delta].mean()
    block = paired.groupby(['method', 'block', 'budget'], as_index=False)[delta].mean()
    contrasts = []
    for method in ['weighted', 'nnls']:
        wide = block.loc[block.method.eq(method)].pivot(index='block', columns='budget', values='delta_mae_pp')
        assert wide.shape == (12, 3) and wide.notna().all().all()
        for label, values in [('I24', wide[24]), ('I60', wide[60]), ('I120', wide[120]), ('C120_60', wide[120] - wide[60]), ('C120_24', wide[120] - wide[24])]:
            contrasts.append({'method': method, 'contrast': label, 'mean_pp': values.mean(), 'mcse_pp': values.std(ddof=1) / np.sqrt(len(values)),
                              'blocks': len(values), 'negative_blocks': int(values.lt(0).sum())})
    precision = pd.DataFrame(contrasts)
    references['allocation'] = np.where(references.dominant.eq('balanced'), 'balanced', '10:1:1')
    absolute = references.groupby(['method', 'budget', 'allocation'], as_index=False)[metrics].mean()
    donor_c = donor.pivot(index=['method', 'held_out'], columns='budget', values='delta_mae_pp').reset_index()
    donor_c.columns = ['method', 'held_out', 'I24_pp', 'I60_pp', 'I120_pp']
    donor_c['C120_60_pp'] = donor_c.I120_pp - donor_c.I60_pp
    donor_c['C120_24_pp'] = donor_c.I120_pp - donor_c.I24_pp
    distribution = arrangements.groupby(['method', 'budget']).delta_mae_pp.agg(
        mean_pp='mean', mean_absolute_pp=lambda x: x.abs().mean(), positive_arrangements=lambda x: x.gt(0).sum(),
        arrangements='size', q10_pp=lambda x: x.quantile(.1), q90_pp=lambda x: x.quantile(.9)).reset_index()
    for tab, name in [(references, 'reference_scores.csv'), (paired, 'paired_scores.csv'), (donor, 'donor_allocation_effects.csv'),
                      (arrangements, 'arrangement_effects.csv'), (block, 'block_effects.csv'), (precision, 'conditional_precision.csv'),
                      (absolute, 'absolute_errors.csv'), (donor_c, 'donor_budget_contrasts.csv'), (distribution, 'arrangement_distribution.csv')]:
        save(tab, name)
    return donor, arrangements, precision, absolute, donor_c, distribution


def plot(donor, arrangements, absolute):
    configure_style(7)
    fig, axes = plt.subplots(2, 2, figsize=(183/25.4, 140/25.4), layout='constrained')
    budgets = [24, 60, 120]
    colors = ['#0072B2', '#D55E00', '#009E73', '#CC79A7']
    donors = sorted(donor.held_out.unique())
    for allocation, marker, color in [('balanced', 'o', '#0072B2'), ('10:1:1', 's', '#D55E00')]:
        values = absolute[(absolute.method == 'weighted') & (absolute.allocation == allocation)].set_index('budget').loc[budgets, 'mae_pp']
        axes[0, 0].plot(budgets, values, marker=marker, color=color, label=allocation)
    axes[0, 0].set(ylabel='Mean absolute error (pp)', title='a  Weighted MuSiC absolute error')
    axes[0, 0].legend()
    for ax, method, title in [(axes[0, 1], 'weighted', 'b  Weighted MuSiC allocation effect'), (axes[1, 0], 'nnls', 'c  Internal NNLS allocation effect')]:
        data = donor[donor.method == method]
        for d, color in zip(donors, colors):
            values = data[data.held_out == d].set_index('budget').loc[budgets, 'delta_mae_pp']
            ax.plot(budgets, values, marker='o', color=color, alpha=.8, label=d)
        ax.plot(budgets, data.groupby('budget').delta_mae_pp.mean().loc[budgets], 'k--s', linewidth=1.3, label='Equal-donor mean')
        ax.axhline(0, color='.5', linewidth=.6, linestyle=':')
        ax.set(ylabel='I(B) = unequal − balanced MAE (pp)', title=title)
    axes[0, 1].legend(ncol=2, fontsize=6)
    ax = axes[1, 1]
    high = arrangements[(arrangements.method == 'weighted') & (arrangements.budget == 120)]
    for i, (d, color) in enumerate(zip(donors, colors)):
        values = high[high.held_out == d].sort_values('dominant').delta_mae_pp.to_numpy()
        ax.scatter(values, np.arange(3)*.14 + i - .14, color=color, s=18)
    ax.axvline(0, color='.5', linewidth=.6, linestyle=':')
    ax.set(yticks=range(4), yticklabels=donors, xlabel='I(120) by dominant-donor arrangement (pp)', title='d  Weighted effects at the largest budget')
    for ax in [axes[0, 0], axes[0, 1], axes[1, 0]]:
        ax.set(xlabel='Reference cells per type across three donors', xticks=budgets)
    for extension in ['png', 'pdf', 'svg']:
        fig.savefig(OUT / f'pancreas_budget_effects.{extension}', dpi=300)
    plt.close(fig)


def write_results(design, precision, absolute, donor_c, distribution):
    lines = ['# Raw-count MuSiC pancreas extension', '',
             'This exploratory cross-tissue extension uses the original annotated UMI counts from Baron et al., GSE84133, DOI 10.1016/j.cels.2016.08.011. The source was previously explored with different normalized-cell workflows; this is not an untouched confirmation cohort.', '',
             'All four source donors and every author cell type with at least 100 cells in each donor were retained: alpha, beta, delta and ductal cells. Each fold holds out one donor and uses the other three as the reference. Budgets are 24, 60 and 120 cells per type. Balanced and 10:1:1 allocations rotate all three dominant donors. Twelve independently seeded global reference blocks use nested cell prefixes within allocation. Each donor supplies 60 fixed 100-cell raw-count targets: one exactly balanced composition and 59 symmetric Dirichlet-multinomial compositions, shared across donors. Cells are sampled without replacement within a target or donor-type reference; reuse across targets, folds and blocks is allowed.', '',
             'MuSiC 1.0.0 at commit f21fe67f5670d5e9fca0ad7550abaae3423eb59c is unchanged. Full-gene counts, markers=NULL, cell_size=NULL, ct.cov=FALSE, iter.max=1000, nu=0.0001, eps=0.01, centered=FALSE and normalize=FALSE match the manuscript workflow. Both weighted and internal ordinary NNLS outputs are retained.', '',
             'The primary endpoint, fixed in data/processed/music_pancreas/design.json before fitting, is weighted I(120)−I(60), where I(B) is 10:1:1-minus-balanced MAE in percentage points. Means equally weight targets, dominant donor choices and held-out donors. MCSE is the SD of twelve complete global-block estimates divided by sqrt(12), conditional on the fixed donor inventories, targets and grouping. Overlapping donor summaries are descriptive; no population P values or confidence intervals are estimated. All twelve blocks were required without effect-direction or precision stopping.', '',
             '## Results', '', '| Output | Budget | Balanced MAE (pp) | 10:1:1 MAE (pp) | I(B) (pp) | Conditional MCSE (pp) |', '| --- | --- | --- | --- | --- | --- |']
    for method in ['weighted', 'nnls']:
        for budget in [24, 60, 120]:
            data = absolute[(absolute.method == method) & (absolute.budget == budget)].set_index('allocation')
            p = precision[(precision.method == method) & (precision.contrast == f'I{budget}')].iloc[0]
            lines.append(f'| {method} | {budget} | {data.loc["balanced", "mae_pp"]:.3f} | {data.loc["10:1:1", "mae_pp"]:.3f} | {p.mean_pp:+.3f} | {p.mcse_pp:.4f} |')
    for method in ['weighted', 'nnls']:
        p = precision[(precision.method == method) & (precision.contrast == 'C120_60')].iloc[0]
        dc = donor_c[donor_c.method == method]
        dist = distribution[(distribution.method == method) & (distribution.budget == 120)].iloc[0]
        lines += ['', f'{method}: I(120)−I(60) = {p.mean_pp:+.3f} pp (conditional MCSE {p.mcse_pp:.4f}); negative in {dc.C120_60_pp.lt(0).sum()}/4 donor means and {p.negative_blocks}/12 blocks. At B=120, {int(dist.positive_arrangements)}/12 arrangements have positive I; mean absolute arrangement effect is {dist.mean_absolute_pp:.3f} pp.']
    primary = precision[(precision.method == 'weighted') & (precision.contrast == 'C120_60')].iloc[0]
    secondary = precision[(precision.method == 'weighted') & (precision.contrast == 'C120_24')].iloc[0]
    wp = precision[precision.method == 'weighted'].set_index('contrast')
    wd = donor_c[donor_c.method == 'weighted']
    dist = distribution[(distribution.method == 'weighted') & (distribution.budget == 120)].iloc[0]
    lines += ['', '## Primary and secondary findings', '',
              f'The prespecified weighted 60-to-120 contrast was {primary.mean_pp:+.3f} pp with conditional MCSE {primary.mcse_pp:.4f}; only {int(primary.negative_blocks)}/12 blocks and {wd.C120_60_pp.lt(0).sum()}/4 donor means were negative. This extension did not reproduce a mean decline over that interval. The small mean relative to reference-randomization variability is not evidence of equivalence or a universal plateau.', '',
              f'The separately specified secondary 24-to-120 contrast was {secondary.mean_pp:+.3f} pp with conditional MCSE {secondary.mcse_pp:.4f}, negative in {int(secondary.negative_blocks)}/12 blocks and {wd.C120_24_pp.lt(0).sum()}/4 donor means. It shows a low-budget penalty in this design, but cannot replace the primary endpoint after seeing the results. The PBMC 60-to-300 endpoint remains untested in pancreas.', '',
              '## Manuscript-ready result', '',
              f'An exploratory pancreas extension showed a lower-budget allocation penalty but did not reproduce the decline over the prespecified 60-to-120-cell interval. Using four held-out donors, four cell types and twelve reference-draw blocks, weighted I(B) was {wp.loc["I24", "mean_pp"]:+.3f}, {wp.loc["I60", "mean_pp"]:+.3f} and {wp.loc["I120", "mean_pp"]:+.3f} percentage points at budgets of 24, 60 and 120 cells per type. The primary contrast, I(120)−I(60), was {primary.mean_pp:+.3f} points (conditional MCSE {primary.mcse_pp:.4f}). The secondary I(120)−I(24) contrast was {secondary.mean_pp:+.3f} points (MCSE {secondary.mcse_pp:.4f}) and was negative for all four donor means. At B=120, {int(dist.positive_arrangements)}/12 arrangements still had higher error under unequal allocation, with a mean absolute contrast of {dist.mean_absolute_pp:.3f} points. These findings extend the observation of budget-dependent allocation sensitivity beyond PBMC, while limiting any claim that the same decline occurs over a common budget interval. The pancreas source had been explored previously with other workflows; its four-type, 100-cell target design and feasible budget range differed from the PBMC analyses.', '',
              '## Interpretation boundaries', '',
              'The budgets, number and identity of cell types, target size, donor population, tissue and platform differ from the PBMC experiments. A common B=300 design is infeasible under the retained four-type inventories. This extension tests allocation sensitivity and its budget dependence within pancreas, not the exact PBMC 60-to-300 endpoint or a causal tissue difference. Four donors and potentially high sampling fractions limit generalization. Synthetic targets do not establish measured-bulk performance. The source is biologically distinct from the PBMC cohorts but was used in earlier project exploration.', '',
              '## Reproduction', '',
              'Use the existing project-local Python environment and the fixed R library described in environment/MUSIC_RUNTIME.md. Input preparation validates the retained original count and annotation hashes against the earlier source manifest; source URLs and download times remain in design.json. It refuses to overwrite prepared inputs. Run analysis/prepare_music_pancreas.py once; from the project root run Rscript --vanilla analysis/run_music_pancreas.R human1 pilot. Review the input and estimate checks with analysis/review_music_pancreas.py --mode pilot before full runs. Run the R script with each human1–human4 and full, then run analysis/review_music_pancreas.py --mode full. The full run reuses completed pilot reference fits. No old PBMC outputs are recomputed or changed.', '',
              'full_validation.json records checks of all new raw-count target sums, reference quotas and donor exclusion, prediction completeness, unit-sum proportions and independent recomputation of error metrics. conditional_precision.csv, donor_budget_contrasts.csv and arrangement_distribution.csv separate block precision, donor summaries and arrangement distributions; paired_scores.csv retains all cell-type errors. Runtime warnings about dependencies built under R 4.5.3 are preserved in command output; actual runtime is R 4.5.2.', '',
              'Figure: Panels a–c use the three numeric reference budgets. Panel a gives equal-donor weighted MAE; panels b and c show four overlapping held-out donor means and their equal-weight mean. Panel d shows each of twelve dominant-donor arrangements at B=120, averaged over targets and twelve blocks. Vertical offsets only separate points. All errors are percentage points; ordinary NNLS is internal to MuSiC. No population uncertainty intervals are shown.']
    (OUT / 'RESULTS.md').write_text('\n'.join(lines) + '\n', encoding='utf-8', newline='\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', choices=['pilot', 'full'], required=True)
    args = parser.parse_args()
    frame, design = read_validate(args.mode)
    if args.mode == 'full':
        donor, arrangements, precision, absolute, donor_c, distribution = summarize(frame, design)
        plot(donor, arrangements, absolute)
        write_results(design, precision, absolute, donor_c, distribution)
        print(precision.to_string(index=False))


if __name__ == '__main__':
    main()
