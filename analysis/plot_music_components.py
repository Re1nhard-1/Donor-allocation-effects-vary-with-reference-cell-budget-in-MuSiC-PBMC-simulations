"""Complete prespecified secondary summaries and plot saved component predictions."""
from pathlib import Path
import itertools
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

root = Path('results/music_components')
d = pd.read_csv(root/'predictions.csv')
keys = ['triple','block','budget','dominant','target_name','held_out']
f = d.pivot(index=keys, columns='mask', values='mae')
# Factorial interaction contrasts: inclusion-exclusion against the 000 baseline.
inter = f.index.to_frame(index=False)
for mask, name in [(3,'Theta_S'),(5,'Theta_Sigma'),(6,'S_Sigma'),(7,'Theta_S_Sigma')]:
    subsets = [s for s in range(8) if s & mask == s]
    inter[name] = sum((-1)**(mask.bit_count()-s.bit_count())*f[s] for s in subsets).to_numpy()
inter.to_csv(root/'factorial_interactions.csv',index=False)
types = ['B','CD14 Mono','CD4 T','CD8 T','FCGR3A Mono','NK']
out = []
for k, label in enumerate(types):
    temp = d[keys+['mask']].copy()
    temp['error'] = np.abs(d[f'pred_{k}']-d[f'true_{k}'])*100
    pf = temp.pivot(index=keys, columns='mask', values='error')
    a = pf.index.to_frame(index=False)
    for j, name in enumerate(['Theta','S','Sigma']):
        phi = np.zeros(len(pf))
        for order in itertools.permutations(range(3)):
            mask = 0
            for factor in order:
                nxt = mask | (1 << factor)
                if factor == j:
                    phi += (pf[nxt]-pf[mask]).to_numpy()/6
                mask = nxt
        a[name] = phi
    a['cell_type'] = label
    out.append(a.groupby(['held_out','block','budget','cell_type'])[['Theta','S','Sigma']].mean().reset_index())
pd.concat(out).to_csv(root/'cell_type_contributions.csv',index=False)
a = pd.read_csv(root/'paired_contributions.csv')
a.groupby(['held_out','triple','dominant','budget'])[['Theta','S','Sigma','common_effect']].mean().to_csv(root/'arrangement_contributions.csv')
checks = [json.loads(p.read_text()) for p in root.glob('*/component_full.json')]
assert len(checks)==14 and all(c['completed'] for c in checks)
control=json.loads((root/'d00_d01_d02/component_control.json').read_text())
counts={'completed_triples':14,'logical_predictions':len(d),'adapter_calls_including_controls':sum(c['actual_adapter_fits'] for c in checks)+control['actual_adapter_fits'], 'additional_official_endpoint_control_calls':control['standard_endpoint_control_fits']}
(root/'execution_counts.json').write_text(json.dumps(counts,indent=2)+'\n')

plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.titlesize':9,'axes.labelsize':8,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
blocks=pd.read_csv(root/'block_budget_contrasts.csv',index_col=0)
fig,axes=plt.subplots(1,2,figsize=(7.2,3.35),gridspec_kw={'width_ratios':[1.2,1]})
labels=['Mean profile Θ','Cell size S','Variance Σ','Total']
for j,metric in enumerate(['Theta','S','Sigma','common_effect']):
    vals=blocks[metric].to_numpy()
    axes[0].scatter(vals,j+np.array([-.10,0,.10]),s=20,facecolors='white',edgecolors='#0072B2',zorder=3)
    axes[0].scatter(vals.mean(),j,s=27,color='black',marker='D',zorder=4)
axes[0].set_yticks(range(4),labels); axes[0].invert_yaxis()
axes[0].axvline(0,color='0.55',lw=.7,ls='--')
axes[0].set_xlabel('Contribution to I(300) − I(60) (pp)')
axes[0].set_title('a  Component contributions',loc='left',weight='bold')
axes[0].scatter([],[],s=20,facecolors='white',edgecolors='#0072B2',label='Three draw blocks')
axes[0].scatter([],[],s=27,color='black',marker='D',label='Mean')
axes[0].legend(frameon=False,loc='lower left',bbox_to_anchor=(0,-.30),ncol=2,fontsize=7)
absolute=pd.read_csv(root/'absolute_errors_by_combination.csv')
for budget,color,marker in [(60,'#D55E00','o'),(300,'#0072B2','s')]:
    v=absolute[absolute.budget.eq(budget)].set_index('mask').loc[range(8),'mae']
    axes[1].plot(v.to_numpy(),range(8),ls='',marker=marker,color=color,ms=4,label=f'B = {budget}')
axes[1].set_yticks(range(8),['000','100','010','110','001','101','011','111'])
axes[1].invert_yaxis(); axes[1].set_ylabel('Unequal component flags (Θ, S, Σ)')
axes[1].set_xlabel('Absolute MAE (pp)')
axes[1].set_title('b  All eight combinations',loc='left',weight='bold')
axes[1].legend(frameon=False,loc='lower left',bbox_to_anchor=(0,-.30),ncol=2,fontsize=7)
fig.subplots_adjust(left=.18,right=.99,top=.89,bottom=.24,wspace=.78)
for ext in ['png','pdf','svg']:
    fig.savefig(root/f'component_interventions.{ext}',dpi=300)
print(json.dumps(counts))
