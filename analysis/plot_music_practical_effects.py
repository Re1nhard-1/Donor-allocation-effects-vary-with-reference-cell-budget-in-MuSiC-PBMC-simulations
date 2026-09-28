import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from scientific_style import configure_style

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'results/music_mc_extension_review/20260918T163322469416Z'
OUT = ROOT / 'results/reference_design_interpretation'


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    metrics = pd.read_csv(SOURCE / 'aggregate_metrics_by_block_set.csv')
    penalties = pd.read_csv(SOURCE / 'aggregate_penalties_by_block_set.csv')
    key = pd.read_csv(SOURCE / 'cell_type_metric_key.csv')
    records = []
    for block_set in ['original_3', 'new_12', 'combined_15']:
        for budget in [60, 300]:
            m = metrics[(metrics.block_set == block_set) & (metrics.method == 'music_weighted') & (metrics.budget == budget)].set_index('level')
            p = penalties[(penalties.block_set == block_set) & (penalties.method == 'music_weighted') & (penalties.budget == budget) & (penalties.level == 'ratio10')].iloc[0]
            fields = [('Overall', 'mae_pp', 'mae_penalty_pp')] + [(r.cell_type, f'absolute_error_{r.type_index}_pp', f'absolute_error_{r.type_index}_penalty_pp') for r in key.itertuples()]
            for label, metric, contrast in fields:
                baseline, unequal, delta = m.loc['balanced', metric], m.loc['ratio10', metric], p[contrast]
                records.append({'block_set': block_set, 'budget': budget, 'cell_type': label, 'balanced_mae_pp': baseline,
                                'unequal_mae_pp': unequal, 'delta_mae_pp': delta, 'relative_error_change_percent': 100 * delta / baseline})
    table = pd.DataFrame(records)
    table.to_csv(OUT / 'pbmc_absolute_and_relative_effects.csv', index=False, lineterminator='\n')
    configure_style(7)
    fig, axes = plt.subplots(1, 2, figsize=(183/25.4, 95/25.4), sharey=True, layout='constrained')
    for ax, budget in zip(axes, [60, 300]):
        data = table[(table.block_set == 'new_12') & (table.budget == budget)].reset_index(drop=True)
        y = np.arange(len(data))
        ax.hlines(y, data.balanced_mae_pp, data.unequal_mae_pp, color='.55', linewidth=1.2)
        ax.scatter(data.balanced_mae_pp, y, color='#0072B2', marker='o', s=23, label='Balanced')
        ax.scatter(data.unequal_mae_pp, y, color='#D55E00', marker='s', s=23, label='10:1:1')
        for i, row in data.iterrows():
            ax.text(10.2, i, f'{row.delta_mae_pp:+.2f} pp ({row.relative_error_change_percent:+.1f}%)', va='center', fontsize=6.3)
        ax.set(xlim=(0, 14.2), xticks=[0, 2, 4, 6, 8, 10], yticks=y, yticklabels=data.cell_type,
               title=f'B = {budget} reference cells per type', xlabel='Mean absolute error (pp)')
        ax.text(10.2, -.75, 'Unequal − balanced', fontsize=6.3)
        ax.set_ylim(6.5, -1)
    axes[0].legend(loc='lower left', bbox_to_anchor=(0, -.25), ncol=2)
    for extension in ['png', 'pdf', 'svg']:
        fig.savefig(OUT / f'pbmc_cell_type_effects.{extension}', dpi=300)
    plt.close(fig)
    caption = ('Weighted MuSiC absolute and relative allocation effects in the twelve additional PBMC reference-draw blocks. '
               'All six cell types and overall MAE are shown in source order. Points give saved equal-donor mean absolute errors under balanced and 10:1:1 allocations. '
               'Annotations show the saved signed allocation contrast in percentage points (pp) and that contrast divided by the balanced error, expressed as a percent. '
               'Relative changes describe error, not changes in cell proportions or application importance. Both panels use the same absolute-error scale; no population intervals are implied. '
               'This is a presentation of existing summaries, with newly calculated relative error ratios, not a new MuSiC fit. The CSV also includes the original three and combined fifteen blocks, kept separate.\n')
    (OUT / 'pbmc_cell_type_effects.caption.md').write_text(caption, encoding='utf-8', newline='\n')
    provenance = {'sources': {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in [SOURCE / 'aggregate_metrics_by_block_set.csv', SOURCE / 'aggregate_penalties_by_block_set.csv', SOURCE / 'cell_type_metric_key.csv']},
                  'operation': 'Reuse saved absolute errors and signed penalties; calculate relative ratios only; plot all types in new_12',
                  'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (OUT / 'source.json').write_text(json.dumps(provenance, indent=2) + '\n', encoding='utf-8', newline='\n')
    print(table[(table.block_set == 'new_12')].to_string(index=False))


if __name__ == '__main__':
    main()
