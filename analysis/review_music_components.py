from pathlib import Path
import hashlib
import itertools
import json
import numpy as np
import pandas as pd

root = Path('results/music_components')
d = pd.read_csv(root / 'predictions.csv')
keys = ['triple', 'block', 'budget', 'dominant', 'target_name', 'held_out']
assert len(d) == 64512 and not d.duplicated(keys + ['mask']).any()
assert d.groupby(keys).size().eq(8).all()
p = d[[f'pred_{k}' for k in range(6)]].to_numpy()
t = d[[f'true_{k}' for k in range(6)]].to_numpy()
assert np.isfinite(p).all() and p.min() >= -1e-10
assert np.max(np.abs(p.sum(axis=1)-1)) < 1e-10
mae = np.abs(p-t).mean(axis=1)*100
mae_diff = float(np.max(np.abs(mae-d.mae)))
assert mae_diff < 1e-10
assert d.convergence.str.fullmatch(r'Converge at [0-9]+').all()
f = d.pivot(index=keys, columns='mask', values='mae').sort_index()
phi = np.zeros((len(f), 3))
for order in itertools.permutations(range(3)):
    mask = 0
    for j in order:
        next_mask = mask | (1 << j)
        phi[:, j] += (f[next_mask]-f[mask]).to_numpy()/6
        mask = next_mask
assert np.max(np.abs(phi.sum(axis=1)-(f[7]-f[0]).to_numpy())) < 1e-10
stored = d[d['mask'].eq(0)].set_index(keys).sort_index()
phi_diff = float(np.max(np.abs(phi-stored[['phi_Theta','phi_S','phi_Sigma']].to_numpy())))
assert phi_diff < 1e-10
a = f.index.to_frame(index=False)
for j, name in enumerate(['Theta', 'S', 'Sigma']):
    a[name] = phi[:, j]
a['common_effect'] = (f[7]-f[0]).to_numpy()
a['original_effect'] = (stored.original_unequal-stored.original_balanced).to_numpy()
a['support_effect'] = a.common_effect-a.original_effect
a.to_csv(root/'paired_contributions.csv',index=False)
metrics = ['Theta','S','Sigma','common_effect','original_effect','support_effect']
donor = a.groupby(['held_out','block','budget'])[metrics].mean().reset_index()
donor.to_csv(root/'donor_block_contributions.csv',index=False)
block = donor.groupby(['block','budget'])[metrics].mean()
contrast = block.xs(300,level='budget')-block.xs(60,level='budget')
contrast.to_csv(root/'block_budget_contrasts.csv')
summary = pd.DataFrame({'mean':contrast.mean(),'block_sd':contrast.std(ddof=1),
                        'min_block':contrast.min(),'max_block':contrast.max()})
summary.to_csv(root/'budget_contribution_summary.csv')
absolute = d.groupby(['held_out','block','budget','mask']).mae.mean().groupby(['budget','mask']).mean()
absolute.to_csv(root/'absolute_errors_by_combination.csv')
checks = {'complete':True,'logical_predictions':len(d),'paired_cases':len(a),
          'mae_max_difference_pp':mae_diff,'shapley_permutation_max_difference_pp':phi_diff,
          'all_converged':True,'common_gene_range':[int(d.n_common.min()),int(d.n_common.max())],
          'review_script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(root/'validation.json').write_text(json.dumps(checks,indent=2)+'\n')
print(summary.to_string())
print(json.dumps(checks))
