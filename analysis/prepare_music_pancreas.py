import gzip
import hashlib
import json
import platform
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy import io, sparse

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data/processed/music_pancreas'
SOURCE = ROOT / 'results/deconv_input/20260916T210631005068Z/input_manifest.json'


def digest(path):
    return hashlib.file_digest(Path(path).open('rb'), 'sha256').hexdigest()


def rng(key):
    seed = int(hashlib.sha256(('music-pancreas-v1|' + key).encode()).hexdigest()[:16], 16)
    return np.random.Generator(np.random.PCG64(seed))


def matrix_write(path, matrix):
    with gzip.open(path, 'wb') as handle:
        io.mmwrite(handle, matrix, field='integer')


def main():
    source = json.loads(SOURCE.read_text())
    OUT.mkdir(parents=True, exist_ok=False)
    donors = [record['donor'] for record in source['records']]
    inventory = pd.DataFrame({r['donor']: r['cell_type_counts'] for r in source['records']}).fillna(0).astype(int)
    types = sorted(inventory.index[(inventory >= 100).all(axis=1)])
    assert types == ['alpha', 'beta', 'delta', 'ductal']
    inventory.to_csv(OUT / 'inventory.csv', lineterminator='\n')
    genes = json.loads((ROOT / source['genes_path']).read_text())
    assert digest(ROOT / source['genes_path']) == source['genes_sha256']
    (OUT / 'genes.tsv').write_text('\n'.join(genes) + '\n', encoding='utf-8', newline='\n')
    matrices, cells, inputs = [], [], []
    for record in source['records']:
        for path_key, hash_key in [('raw_path', 'sha256'), ('matrix_path', 'matrix_sha256'), ('cells_path', 'cells_sha256')]:
            assert digest(ROOT / record[path_key]) == record[hash_key]
        obs = pd.read_csv(ROOT / record['cells_path'])
        counts = sparse.load_npz(ROOT / record['matrix_path']).tocsr()
        assert counts.shape == tuple(record['shape'])
        np.testing.assert_array_equal(np.asarray(counts.sum(axis=1)).ravel(), obs.umi_total)
        keep = obs.cell_type.isin(types).to_numpy()
        counts, obs = counts[keep], obs.loc[keep].copy()
        obs['source_cell_id'] = obs.cell_id
        obs['cell_id'] = record['donor'] + ':' + obs.cell_id.astype(str)
        cells.append(obs)
        matrices.append(counts)
        inputs.append({k: record[k] for k in ['donor', 'accession', 'url', 'retrieved_at_utc', 'raw_path', 'sha256', 'matrix_path', 'matrix_sha256', 'cells_path', 'cells_sha256']})
    cells = pd.concat(cells, ignore_index=True)
    cells['column_R'] = np.arange(1, len(cells) + 1)
    assert cells.cell_id.is_unique
    counts = sparse.vstack(matrices, format='csr')
    assert np.isfinite(counts.data).all() and (counts.data >= 0).all()
    assert np.equal(counts.data, np.floor(counts.data)).all()
    cells.to_csv(OUT / 'cells.tsv', sep='\t', index=False, lineterminator='\n')
    matrix_write(OUT / 'counts.mtx.gz', counts.T.tocsc())
    pools = {(d, t): cells.index[(cells.donor == d) & (cells.cell_type == t)].to_numpy() for d in donors for t in types}
    mixtures = [np.full(len(types), 25)]
    composition_rng = rng('compositions')
    for composition in composition_rng.dirichlet(np.ones(len(types)), 59):
        mixtures.append(composition_rng.multinomial(100, composition))
    mixtures = np.asarray(mixtures)
    target_rows, target_matrices = [], []
    for donor in donors:
        for mixture_id, n_cells in enumerate(mixtures):
            selected = np.concatenate([rng(f'target|{donor}|{mixture_id}|{t}').choice(pools[donor, t], n, replace=False) for t, n in zip(types, n_cells)])
            assert len(selected) == len(np.unique(selected)) == 100
            assert set(cells.iloc[selected].donor) == {donor}
            vector = sparse.csc_matrix(counts[selected].sum(axis=0).T)
            target_matrices.append(vector)
            target_rows.append({'target_name': f'{donor}_m{mixture_id:02d}', 'held_out': donor, 'mixture_id': mixture_id,
                                'columns_R': '|'.join(map(str, selected + 1)), **{f'true_{t}': n / 100 for t, n in zip(types, n_cells)}})
    pd.DataFrame(target_rows).to_csv(OUT / 'targets.tsv', sep='\t', index=False, lineterminator='\n')
    matrix_write(OUT / 'targets.mtx.gz', sparse.hstack(target_matrices, format='csc'))
    references = []
    for block in range(12):
        lists = {(d, t): rng(f'reference|{block}|{d}|{t}').permutation(pools[d, t]) for d in donors for t in types}
        for held_out in donors:
            reference_donors = [d for d in donors if d != held_out]
            for budget in [24, 60, 120]:
                for dominant in ['balanced'] + reference_donors:
                    quotas = {d: budget // 3 if dominant == 'balanced' else (budget * 10 // 12 if d == dominant else budget // 12) for d in reference_donors}
                    selected = np.concatenate([lists[d, t][:quotas[d]] for d in reference_donors for t in types])
                    assert len(selected) == len(np.unique(selected)) == budget * len(types)
                    assert held_out not in set(cells.iloc[selected].donor)
                    references.append({'reference_id': f'{held_out}_b{block:02d}_n{budget}_{dominant}', 'held_out': held_out,
                                       'block': block, 'budget': budget, 'dominant': dominant,
                                       'reference_donors': '|'.join(reference_donors), 'columns_R': '|'.join(map(str, selected + 1))})
    pd.DataFrame(references).to_csv(OUT / 'references.tsv', sep='\t', index=False, lineterminator='\n')
    config = {'version': 'music-pancreas-v1', 'scope': 'Exploratory cross-tissue raw-count extension; previously explored source, not untouched validation',
              'source_accession': 'GSE84133', 'paper_doi': '10.1016/j.cels.2016.08.011', 'source_manifest': SOURCE.relative_to(ROOT).as_posix(),
              'source_manifest_sha256': digest(SOURCE), 'source_inputs': inputs, 'donors': donors, 'cell_types': types,
              'eligibility': 'All four source donors; all author cell types with at least 100 cells in every donor',
              'budgets': [24, 60, 120], 'allocations': ['1:1:1', '10:1:1'], 'reference_blocks': 12,
              'target_cells': 100, 'targets_per_donor': 60, 'target_distribution': 'One exactly equal composition; 59 symmetric Dirichlet(1) then multinomial(100); same compositions across donors',
              'sampling': 'Without replacement within target and donor-type reference; target donor excluded in each fold; shared cells across targets, folds and blocks allowed; global independently seeded block permutations with nested prefixes within allocation',
              'primary_contrast': 'Weighted I(120)-I(60), I(B)=MAE(10:1:1)-MAE(balanced), percentage points; equal means over targets, dominant choices and held-out donors',
              'secondary': ['I(24)', 'I(120)-I(24)', 'absolute MAE and RMSE', 'all four type-specific errors', '12 arrangement effects', 'internal ordinary NNLS'],
              'inference': 'Descriptive donor and arrangement results; conditional MCSE across 12 whole blocks; no population p-values or CI; no precision-based stopping',
              'pilot_gate': 'First held-out donor, block 0, all 12 references; complete finite nonnegative unit-sum estimates, gene identity and quota checks, convergence diagnostics and examination for whole-lineage collapse; retain pilot predictions in full analysis; do not select using allocation-effect sign',
              'music_version': '1.0.0', 'music_commit': 'f21fe67f5670d5e9fca0ad7550abaae3423eb59c', 'music_license': 'GPL >= 3',
              'parameters': {'markers': None, 'cell_size': None, 'ct.cov': False, 'iter.max': 1000, 'nu': 0.0001, 'eps': 0.01, 'centered': False, 'normalize': False},
              'source_terms': 'Public GEO annotated count tables; retain original accession, authors and citation; no additional dataset license inferred from public access',
              'software': {'python': platform.python_version(), 'numpy': np.__version__, 'pandas': pd.__version__, 'scipy': scipy.__version__},
              'seed': 'PCG64 initialized from first 16 SHA256 hex digits of music-pancreas-v1|purpose|identifiers',
              'limitations': ['Four donors with overlapping reference pools', 'One donor recorded as type 2 diabetes; no disease contrast', 'Finite inventories and high sampling fractions', 'Different tissue, platform, retained types and 100-cell targets versus PBMC; no between-source causal comparison', 'B=300 not feasible for the common four types; no test of the PBMC upper endpoint', 'Synthetic targets only; no measured-bulk validation', 'No causal mechanism identification'],
              'scripts': {name: digest(ROOT / 'analysis' / name) for name in ['prepare_music_pancreas.py', 'run_music_pancreas.R']},
              'files': {p.name: digest(p) for p in OUT.iterdir() if p.is_file()}}
    (OUT / 'design.json').write_text(json.dumps(config, indent=2) + '\n', encoding='utf-8', newline='\n')
    print(json.dumps({'cells': len(cells), 'genes': len(genes), 'references': len(references), 'targets': len(target_rows), 'types': types}))


if __name__ == '__main__':
    main()
