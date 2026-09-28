"""Download original Baron human count tables and validate sparse inputs."""
from __future__ import annotations
import argparse
from collections import Counter
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import platform
from urllib.request import Request, urlopen
import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]

def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False)+"\n", encoding="utf-8", newline="\n")

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("metadata_audit",type=Path)
    args=parser.parse_args()
    metadata_path=args.metadata_audit.resolve()
    metadata=json.loads(metadata_path.read_bytes())
    run=datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    raw=ROOT/'data/raw/GSE84133'/run
    processed=ROOT/'data/processed/baron'/run
    output=ROOT/'results/deconv_input'/run
    for p in (raw,processed,output): p.mkdir(parents=True,exist_ok=False)
    records=[]
    expected_genes=None
    for i,acc in enumerate(['GSM2230757','GSM2230758','GSM2230759','GSM2230760'],1):
        donor=f'human{i}'
        meta=next(r for r in metadata['records'] if r['accession']==acc)
        source=next(u for u in meta['parsed']['supplementary_files'] if u.endswith('_counts.csv.gz'))
        url=source.replace('ftp://ftp.ncbi.nlm.nih.gov/','https://ftp.ncbi.nlm.nih.gov/',1)
        if not url.startswith('https://ftp.ncbi.nlm.nih.gov/geo/'): raise ValueError(url)
        target=raw/url.rsplit('/',1)[1]
        request=Request(url,headers={'User-Agent':'bioinformatics-paper-pilot/0.1'})
        total=0
        with urlopen(request,timeout=45) as response, target.open('xb') as handle:
            final_url=response.geturl()
            while chunk:=response.read(1024*1024):
                total+=len(chunk)
                if total>128*1024*1024: raise ValueError('Unexpected download size')
                handle.write(chunk)
        print(f'{donor}: downloaded {total} bytes',flush=True)
        chunks=[]; observation=[]; counts_by_type=Counter()
        with gzip.open(target,'rt') as handle:
            header=pd.read_csv(handle,nrows=0).columns.tolist()
        if header[1:3]!=['barcode','assigned_cluster']:
            raise ValueError(f'Unexpected annotation columns: {header[:3]}')
        genes=header[3:]
        if len(set(genes))!=len(genes): raise ValueError('Duplicate features')
        if expected_genes is None: expected_genes=genes
        if genes!=expected_genes: raise ValueError('Gene order differs among donors')
        for frame in pd.read_csv(target,chunksize=128):
            values=frame.iloc[:,3:].to_numpy()
            if not np.isfinite(values).all() or (values<0).any() or not np.equal(values,np.floor(values)).all():
                raise ValueError('Counts are not finite nonnegative integers')
            if values.max()>np.iinfo(np.int32).max: raise ValueError('Count exceeds int32')
            counts=sparse.csr_matrix(values.astype(np.int32))
            if (np.asarray(counts.sum(axis=1)).ravel()<=0).any(): raise ValueError('Zero UMI cell')
            obs=frame.iloc[:,:3].copy()
            obs.columns=['cell_id','barcode','cell_type']
            if obs.isna().any().any(): raise ValueError('Missing cell annotation')
            obs.insert(0,'donor',donor)
            obs['umi_total']=np.asarray(counts.sum(axis=1)).ravel()
            observation.append(obs); chunks.append(counts)
            counts_by_type.update(obs['cell_type'])
        matrix=sparse.vstack(chunks,format='csr')
        obs=pd.concat(observation,ignore_index=True)
        if obs['cell_id'].duplicated().any(): raise ValueError('Repeated cell ID within donor')
        matrix_path=processed/f'{donor}_counts.npz'
        obs_path=processed/f'{donor}_cells.csv'
        sparse.save_npz(matrix_path,matrix)
        obs.to_csv(obs_path,index=False,lineterminator='\n')
        records.append({'donor':donor,'accession':acc,'url':url,'deposited_url':source,
                        'final_url':final_url,'retrieved_at_utc':datetime.now(timezone.utc).isoformat(),
                        'raw_path':target.relative_to(ROOT).as_posix(),'bytes':total,'sha256':digest(target),
                        'matrix_path':matrix_path.relative_to(ROOT).as_posix(),'matrix_sha256':digest(matrix_path),
                        'cells_path':obs_path.relative_to(ROOT).as_posix(),'cells_sha256':digest(obs_path),
                        'shape':list(matrix.shape),'nonzero_entries':int(matrix.nnz),
                        'cell_type_counts':dict(sorted(counts_by_type.items())),
                        'median_umi':float(obs['umi_total'].median()),
                        'characteristics':meta['parsed'].get('characteristics_ch1',[])})
        print(f'{donor}: {matrix.shape}; {dict(sorted(counts_by_type.items()))}',flush=True)
    genes_path=processed/'genes.json'
    write_json(genes_path,expected_genes)
    report={'run_id':run,'dataset':'GSE84133 original annotated human UMI-filtered counts',
            'source_paper':'https://pmc.ncbi.nlm.nih.gov/articles/PMC5228327/',
            'metadata_audit':metadata_path.relative_to(ROOT).as_posix(),'metadata_sha256':digest(metadata_path),
            'script':'analysis/prepare_baron.py','script_sha256':digest(__file__),
            'python':platform.python_version(),'genes_path':genes_path.relative_to(ROOT).as_posix(),
            'genes_sha256':digest(genes_path),'records':records,
            'limitations':['Original author cell labels retained; annotation uncertainty not validated.',
                'Four source donors, one recorded type 2 diabetes; no causal disease comparison.',
                'Processed count tables, not new wet-lab measurements or new biological samples.',
                'Public deposited inputs; cite original creators; publication-specific terms still to review.']}
    write_json(output/'input_manifest.json',report)
    pd.DataFrame([{'donor':r['donor'],**r['cell_type_counts']} for r in records]).fillna(0).to_csv(
        output/'cell_type_inventory.csv',index=False,lineterminator='\n')
    print('MANIFEST='+str(output/'input_manifest.json'),flush=True)

if __name__=='__main__': main()
