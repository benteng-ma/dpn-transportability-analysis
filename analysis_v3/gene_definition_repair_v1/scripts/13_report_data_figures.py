from common import *
import pandas as pd,numpy as np,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

LABELS={'late_shared_concordant_neuronal_core':'Late core','late_neuron_residual':'Neuron residual','late_allcell_residual':'All-cell residual','severity_neuron_shared_concordant_core':'Severity core','severity_neuron_residual':'Severity residual','original_early_allcell':'Early parent','original_late_allcell':'All-cell parent','original_late_neuron':'Neuron parent','original_severity':'Severity parent','original_xenium':'Xenium parent'}
ENDS={'original_hDRG_DPN_control':'Human hDRG: DPN - control','original_JCI_sural_DPN_control':'Sural nerve: DPN - control','original_JCI_sural_severe_moderate':'Sural nerve: severe - moderate','GSE24290_C01':'GSE24290: progressive - non-progressive','GSE148059_C02':'GSE148059: degenerator - regenerator'}
def read(f):
 override=B/'strict_background_addendum/consolidated_results'/Path(f).name
 return pd.read_csv(override if str(f).startswith('results/') and override.exists() else B/f,sep='\t')
def save(fig,name):
 fig.savefig(B/'figures'/f'{name}.png',dpi=300,bbox_inches='tight',facecolor='white')
 fig.savefig(B/'figures'/f'{name}.tiff',dpi=400,bbox_inches='tight',facecolor='white',pil_kwargs={'compression':'tiff_lzw'})
 plt.close(fig)
def main():
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'text.color':'black','axes.labelcolor':'black','xtick.color':'black','ytick.color':'black','axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
 mem=read('inputs/REPAIRED_MEMBERS.tsv');old=read('results/OLD_NEW_MEMBER_ROW_COMPARISON.tsv');cv=read('results/COVERAGE.tsv');ass=read('results/ASSOCIATION_IMPACT_ALL.tsv');hist=read('results/HISTORICAL_DONOR_TEST_IMPACT.tsv');be=read('results/REAL_BENCHMARK_ALL.tsv');rnd=read('results/MATCHED_RANDOM_REFERENCE.tsv')
 counts=mem.groupby(['module_id','direction']).size().unstack(fill_value=0).reset_index();counts['total']=counts.up+counts.down;put(counts,'results/REPAIRED_PROGRAM_COUNTS.tsv')
 # Flat machine-readable comparison: each endpoint/program has three methods explicitly labeled.
 keys=['endpoint','module_id'];wide=None
 for ver in ['legacy_inherited','repaired_inherited','repaired_statfree']:
  f=ass[ass.version==ver].drop(columns='version').merge(cv[cv.version==ver].drop(columns='version'),on=keys,validate='one_to_one')
  f=f.rename(columns={c:ver+'__'+c for c in f.columns if c not in keys});wide=f if wide is None else wide.merge(f,on=keys,validate='one_to_one')
 put(wide,'results/MAIN_EFFECT_COVERAGE_COMPARISON.tsv')
 scores=read('results/DONOR_SCORES.tsv.gz');sw=scores.pivot(index=['endpoint','sample_id','patient_id','group','module_id'],columns='version',values='score').reset_index();put(sw,'results/DONOR_SCORE_COMPARISON.tsv.gz')
 # All selected display rows are entire fixed families, not significance-selected figures.
 part=read('results/PARTITION_COUNTS.tsv');pt=part.pivot(index='family',columns='bin',values='genes').fillna(0).reindex(['late_allcell','late_neuron','severity'])
 fig,ax=plt.subplots(figsize=(7,3.5));left=np.zeros(3)
 for col,label,color in [('core','Core','#31688e'),('residual','Residual','#9ab5a5'),('excluded_cross_source_opposed','Opposed source directions','#b58b78'),('excluded_partner_unresolved','Unresolved partner','#d4d4d4')]:
  vals=pt[col].to_numpy();ax.barh(range(3),vals,left=left,color=color,label=label,height=.58)
  for j,(l,v) in enumerate(zip(left,vals)):
   if v>=75:ax.text(l+v/2,j,str(int(v)),ha='center',va='center',fontsize=8)
  left+=vals
 ax.set_yticks(range(3),['Late all-cell','Late neuron','Severity']);ax.invert_yaxis();ax.set_xlabel('Unique GeneIDs in repaired signed parent');ax.legend(loc='upper center',bbox_to_anchor=(.5,-.22),ncol=2,frameon=False);fig.tight_layout();save(fig,'Figure_R1_partition')
 comp=list(LABELS)[:5];allmods=comp+list(LABELS)[5:]
 neural=[e for e in ass.endpoint.unique() if e.startswith('original_')]
 fig,axes=plt.subplots(1,3,figsize=(10.6,4.7),sharey=True)
 for k,(ax,eid) in enumerate(zip(axes,neural)):
  for j,mod in enumerate(comp):
   for ver,shift,color,mark in [('legacy_inherited',-.13,'#777777','o'),('repaired_statfree',.13,'#31688e','s')]:
    z=ass[(ass.endpoint==eid)&(ass.module_id==mod)&(ass.version==ver)].iloc[0]
    if z.status=='COMPLETED':ax.errorbar(z.mean_difference,j+shift,xerr=[[z.mean_difference-z.ci_low],[z.ci_high-z.mean_difference]],fmt=mark,color=color,ms=4,lw=1,capsize=2)
    else:ax.text(.98,j+shift,'NE',ha='right',va='center',transform=ax.get_yaxis_transform(),fontsize=8)
  ax.axvline(0,color='black',lw=.6);ax.set_xlabel('Mean score difference (95% CI)');ax.text(-.12,1.02,chr(65+k),transform=ax.transAxes,fontweight='bold',fontsize=12)
 axes[0].set_yticks(range(5),[LABELS[m] for m in comp]);axes[0].invert_yaxis();fig.tight_layout();save(fig,'Figure_R2_neural_effects')
 alc=[e for e in ass.endpoint.unique() if e.startswith('GSE')]
 fig,axes=plt.subplots(1,2,figsize=(9,5.5),sharey=True)
 for k,(ax,eid) in enumerate(zip(axes,alc)):
  for j,mod in enumerate(allmods):
   for ver,shift,color,mark in [('legacy_inherited',-.13,'#777777','o'),('repaired_statfree',.13,'#31688e','s')]:
    z=ass[(ass.endpoint==eid)&(ass.module_id==mod)&(ass.version==ver)].iloc[0]
    if z.status=='COMPLETED':ax.errorbar(z.mean_difference,j+shift,xerr=[[z.mean_difference-z.ci_low],[z.ci_high-z.mean_difference]],fmt=mark,color=color,ms=3.5,lw=.8,capsize=2)
    else:ax.text(.98,j+shift,'NE',ha='right',va='center',transform=ax.get_yaxis_transform(),fontsize=8)
  ax.axvline(0,color='black',lw=.6);ax.set_xlabel('Mean score difference (95% CI)');ax.text(-.12,1.02,chr(65+k),transform=ax.transAxes,fontweight='bold',fontsize=12)
 axes[0].set_yticks(range(10),[LABELS[m] for m in allmods]);axes[0].invert_yaxis();fig.tight_layout();save(fig,'Figure_R3_ALC_effects')
 rr=rnd[rnd.background=='statfree'].copy();put(rr,'figures/Figure_R4_source.tsv')
 fig,axes=plt.subplots(1,2,figsize=(10.5,5),sharey=True)
 names=[]
 for j,z in enumerate(rr.itertuples()):
  short={'original_hDRG_DPN_control':'hDRG DPN','original_JCI_sural_DPN_control':'Sural DPN','original_JCI_sural_severe_moderate':'Sural axonal loss','GSE24290_C01':'GSE24290','GSE148059_C02':'GSE148059'}
  names.append(short[z.endpoint]+' / '+z.family.replace('late_',''))
  row=be[(be.endpoint==z.endpoint)&(be.background=='statfree')&(be.family==z.family)&(be.test=='core_minus_residual')].iloc[0]
  axes[0].errorbar(row.effect,j,xerr=[[row.effect-row.ci_low],[row.ci_high-row.effect]],fmt='o',color='#31688e',ms=4,capsize=2)
  axes[1].plot(z.conditional_midrank_percentile,j,'o',color='#31688e')
 axes[0].axvline(0,color='black',lw=.6);axes[0].set_xlabel('Core - residual group contrast (95% CI)');axes[1].set_xlim(0,1);axes[1].set_xlabel('Conditional matched-split percentile');axes[1].axvline(.5,color='#999999',lw=.6,ls=':')
 axes[0].set_yticks(range(len(rr)),names);axes[0].invert_yaxis()
 for ax,label in zip(axes,'AB'):ax.text(-.12,1.02,label,transform=ax.transAxes,fontweight='bold',fontsize=12)
 fig.tight_layout();save(fig,'Figure_R4_benchmarks')
 for name,f in [('Figure_R1',part),('Figure_R2',ass[ass.endpoint.isin(neural)]),('Figure_R3',ass[ass.endpoint.isin(alc)])]:put(f,'figures/'+name+'_source.tsv')
 audit=json.loads((B/'reports/FINAL_AUDIT.json').read_text());summary={'created_utc':now(),'main_association_rows':len(ass),'donor_score_rows':len(scores),'corrected_donor_scores':int(scores.version.str.startswith('repaired').sum()),'benchmark_rows':len(be),'benchmark_complete':int(be.status.eq('COMPLETED').sum()),'benchmark_holm_below_05':int((be.holm_p<.05).sum()),'minimum_benchmark_holm':float(be.holm_p.min()),'matched_partitions':int(rnd.draws.sum()),'historical_tests':len(hist),'tests_passed':7234+audit['counts']['PASS'],'repaired_memberships':len(mem),'secondary':audit['secondary_result_families']}
 summary['tests_passed']=json.loads((B/'reports/HANDOFF_CHECKS.json').read_text())['all_passes']
 for name in summary['secondary']:
  f=read('results/'+name);summary['secondary'][name]=dict(rows=len(f),finite_q=int(f.q_new.notna().sum()),q_below_05=int((f.q_new<.05).sum()),q_below_10=int((f.q_new<.1).sum()),minimum_q=float(f.q_new.min()))
 dump(summary,'reports/REPORT_NUMBERS.json');print(json.dumps(summary,indent=2),flush=True)
if __name__=='__main__':run(main)
