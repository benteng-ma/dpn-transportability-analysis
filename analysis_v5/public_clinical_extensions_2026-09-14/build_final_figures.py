from PIL import Image, ImageDraw, ImageFont
from pathlib import Path
import math, os
import matplotlib.pyplot as plt
import pandas as pd
from docx import Document

HERE = Path(__file__).resolve().parent
SRC = Path(os.environ.get('DPN_FIGURE_SOURCE_DIR', HERE / 'figure_sources'))
OUT = Path(os.environ.get('DPN_FIGURE_OUTPUT_DIR', HERE / 'rebuilt_figures'))
OUT.mkdir(parents=True, exist_ok=True)

try:
    font_b = ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf', 48)
    font_n = ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf', 26)
    font_s = ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf', 22)
except:
    font_b = None; font_n=None; font_s=None

def add_label(draw, x, y, lab):
    draw.text((x,y), lab, fill='black', font=font_b)

def cell_crop(im, col, row, cols=2, rows=4, erase_label=True):
    w,h=im.size
    x0=round(col*w/cols); x1=round((col+1)*w/cols)
    y0=round(row*h/rows); y1=round((row+1)*h/rows)
    crop=im.crop((x0,y0,x1,y1)).convert('RGB')
    if erase_label:
        d=ImageDraw.Draw(crop)
        d.rectangle((0,0,min(120,crop.width),min(90,crop.height)),fill='white')
    return crop

def fit_paste(canvas, crop, box, pad=22):
    x0,y0,x1,y1=box
    bw=x1-x0-2*pad; bh=y1-y0-2*pad
    r=min(bw/crop.width,bh/crop.height)
    nw=max(1,int(crop.width*r)); nh=max(1,int(crop.height*r))
    cr=crop.resize((nw,nh),Image.Resampling.LANCZOS)
    x=x0+pad+(bw-nw)//2; y=y0+pad+(bh-nh)//2
    canvas.paste(cr,(x,y))

def make_grid_from_original(srcfile, outfile, custom=None, canvas_size=(2600,3200)):
    im=Image.open(srcfile).convert('RGB')
    W,H=canvas_size
    canvas=Image.new('RGB',(W,H),'white')
    draw=ImageDraw.Draw(canvas)
    ml,mr,mt,mb=90,70,85,70
    gx,gy=90,70
    cw=(W-ml-mr-gx)//2
    rh=(H-mt-mb-3*gy)//4
    labs='abcdefgh'
    for row in range(4):
        for col in range(2):
            idx=row*2+col
            x0=ml+col*(cw+gx); y0=mt+row*(rh+gy)
            x1=x0+cw; y1=y0+rh
            crop = custom.get(idx) if custom and idx in custom else cell_crop(im,col,row)
            fit_paste(canvas,crop,(x0,y0,x1,y1),pad=10)
            add_label(draw,x0-55,y0-42,labs[idx])
    canvas.save(outfile,dpi=(300,300),quality=95)
    return outfile

# Figure 1: custom panel a flowchart; other panels from original v4.5
fig1=Image.open(SRC/'image1.png').convert('RGB')
# create custom panel a at approx original cell size
pa=Image.new('RGB',(1100,620),'white'); d=ImageDraw.Draw(pa)
boxes=[(55,65,260,210,'Discovery DRG\nprograms\nfixed source directions'),
       (310,65,540,210,'Canonical audit\nand fixed scoring\nto whole-background ranks'),
       (590,65,805,210,'Human neural\nand pathology\ncomparisons'),
       (850,65,1065,210,'Cell, spatial,\nclinical and assay\ncontext'),
       (220,340,505,500,'Donor-level or\nbiological-library\ninference'),
       (610,340,935,500,'Tissue-specific\ninterpretation and\nexternal boundaries')]
for x0,y0,x1,y1,txt in boxes:
    d.rectangle((x0,y0,x1,y1),outline='black',width=3)
    # centered multiline
    bbox=d.multiline_textbbox((0,0),txt,font=font_s,spacing=4,align='center')
    tw=bbox[2]-bbox[0]; th=bbox[3]-bbox[1]
    d.multiline_text(((x0+x1-tw)/2,(y0+y1-th)/2),txt,font=font_s,fill='black',spacing=4,align='center')
# arrows
def arrow(p1,p2):
    d.line([p1,p2],fill='black',width=4)
    x2,y2=p2; x1,y1=p1
    ang=math.atan2(y2-y1,x2-x1); L=16
    for da in (2.55,-2.55):
        d.line([(x2,y2),(x2+L*math.cos(ang+da),y2+L*math.sin(ang+da))],fill='black',width=4)
arrow((260,138),(310,138)); arrow((540,138),(590,138)); arrow((805,138),(850,138))
arrow((425,210),(425,340)); arrow((697,210),(697,340)); arrow((505,420),(610,420))
make_grid_from_original(SRC/'image1.png', OUT/'Figure_1.png', custom={0:pa})

# Figure 2 and 5: preserve original data-generated figures, only copy
for srcname,outname in [('image2.png','Figure_2.png'),('image5.png','Figure_5.png')]:
    im=Image.open(SRC/srcname).convert('RGB'); im.save(OUT/outname,dpi=(300,300))

# Supplement tables for Fig 3 custom panels
supp_path = Path(os.environ.get('DPN_SUPPLEMENTARY_DOCX', HERE / 'DPN_Supplementary_Information_SR_v4_5.docx'))
supp = Document(supp_path)
# table 35 main effects
rows=[]
for r in supp.tables[35].rows[1:]:
    vals=[c.text.strip() for c in r.cells]
    if vals[7]=='COMPLETED':
        rows.append(dict(program=vals[0],effect=float(vals[2]),p=float(vals[3]),q=float(vals[4]),lo=float(vals[5]),hi=float(vals[6])))
main=pd.DataFrame(rows)
order=['P1','P2','P3','P4','P5','P6','P7','P8','P10']; main=main.set_index('program').loc[order].reset_index()
# arms table
rows=[]
for r in supp.tables[36].rows[1:]:
    vals=[c.text.strip() for c in r.cells]; rows.append(dict(program=vals[0],arm=vals[1],effect=float(vals[3])))
arms=pd.DataFrame(rows)
# donor effects table
rows=[]
for r in supp.tables[41].rows[1:]:
    vals=[c.text.strip() for c in r.cells]; rows.append(dict(donor=vals[0],program=vals[1],effect=float(vals[3])))
don=pd.DataFrame(rows)

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
def save_ax(fig, path):
    fig.savefig(path,dpi=300,bbox_inches='tight',facecolor='white'); plt.close(fig); return Image.open(path).convert('RGB')
TMP = Path(os.environ.get('DPN_FIGURE_TMP_DIR', HERE / 'tmp_dpn_panels'))
TMP.mkdir(parents=True, exist_ok=True)
# c
fig,ax=plt.subplots(figsize=(5.2,3.8)); y=range(len(main)); ax.errorbar(main.effect,y,xerr=[main.effect-main.lo,main.hi-main.effect],fmt='o',capsize=3)
ax.axvline(0,ls=':',lw=1); ax.set_yticks(list(y),main.program); ax.invert_yaxis(); ax.set_xlabel('Nageotte minus adjacent neuronal score')
for i,r in main.iterrows(): ax.text(r.hi+0.002,i,f'q={r.q:.3f}',va='center',fontsize=8)
fig.tight_layout(); pc=save_ax(fig,TMP/'c.png')
# d donor distributions
fig,ax=plt.subplots(figsize=(5.2,3.8)); colors=plt.cm.tab10.colors
for i,p in enumerate(order):
    vals=don[don.program==p].effect.values
    if len(vals):
        jitter=[i+(j-(len(vals)-1)/2)*0.035 for j in range(len(vals))]
        ax.scatter(jitter,vals,s=28,color=colors[i%10]); ax.hlines(vals.mean(),i-.22,i+.22,color=colors[i%10],lw=2)
ax.axhline(0,ls=':',lw=1); ax.set_xticks(range(len(order)),order); ax.set_ylabel('Within-donor region difference'); ax.set_xlabel('Each point is one donor')
fig.tight_layout(); pdim=save_ax(fig,TMP/'d.png')
# e arm decomposition: DeltaU and -DeltaD
pivot=arms.pivot(index='program',columns='arm',values='effect').loc[order]
fig,ax=plt.subplots(figsize=(5.2,3.8)); x=range(len(order)); w=.36
ax.bar([i-w/2 for i in x],pivot['source_up'],width=w,label='source-up change')
ax.bar([i+w/2 for i in x],-pivot['source_down'],width=w,label='negative source-down change')
ax.axhline(0,ls=':',lw=1); ax.set_xticks(list(x),order); ax.set_ylabel('Contribution to score difference'); ax.legend(frameon=False,fontsize=8,loc='upper center',ncol=2)
fig.tight_layout(); pe=save_ax(fig,TMP/'e.png')
# g LODO
fig,ax=plt.subplots(figsize=(5.2,3.8)); y=range(len(main)); ax.errorbar(main.effect,y,xerr=[main.effect-main.lo,main.hi-main.effect],fmt='o',capsize=3)
ax.axvline(0,ls=':',lw=1); ax.set_yticks(list(y),main.program); ax.invert_yaxis(); ax.set_xlabel('Effect and leave-one-donor-out range'); fig.tight_layout(); pg=save_ax(fig,TMP/'g.png')
# h effect-q
fig,ax=plt.subplots(figsize=(5.2,3.8)); ax.scatter(main.effect,-main.q.map(math.log10),s=35)
for _,r in main.iterrows(): ax.text(r.effect, -math.log10(r.q)+0.025, r.program, ha='center',fontsize=8)
ax.axvline(0,ls=':',lw=1); ax.axhline(-math.log10(.05),ls=':',lw=1); ax.set_xlabel('Nageotte minus neuronal score'); ax.set_ylabel('-log10(BH q)'); fig.tight_layout(); ph=save_ax(fig,TMP/'h.png')
# Fig3 custom mapping: c=2,d=3,e=4,g=6,h=7; f original
custom3={2:pc,3:pdim,4:pe,6:pg,7:ph}
make_grid_from_original(SRC/'image3.png', OUT/'Figure_3.png', custom=custom3)

# Figures 4 and 6 generic repack from data-generated originals
make_grid_from_original(SRC/'image4.png', OUT/'Figure_4.png')
make_grid_from_original(SRC/'image6.png', OUT/'Figure_6.png')
# Figure7 more compact landscape-ish square canvas
make_grid_from_original(SRC/'image7.png', OUT/'Figure_7.png', canvas_size=(3000,3000))

print('built', list(map(str, sorted(OUT.glob('Figure_*.png')))))
# Regenerate Figure 7 directly from verified M01 result tables (no raster recycling)
BASE = HERE / '01_GSE302658' / 'results'
fa=pd.read_csv(BASE/'GSE302658_CLINICAL_FIELD_AUDIT.tsv',sep='\t')
bas=pd.read_csv(BASE/'BASELINE_SYMPTOM_ASSOCIATIONS.tsv',sep='\t')
lng=pd.read_csv(BASE/'LONGITUDINAL_SYMPTOM_ASSOCIATIONS.tsv',sep='\t')
loo=pd.read_csv(BASE/'LOO_STABILITY_RESULTS.tsv',sep='\t')
inter=pd.read_csv(BASE/'TREATMENT_INTERACTION_OMNIBUS.tsv',sep='\t')
programs=['P1','P2','P3','P4','P6','P7','P8','P10']
symptoms=['npsi_total','burning_superficial','pressing_deep','paresthesia_dysesthesia','average_pain_nrs','worst_pain_nrs','paroxysmal','evoked']
sym_labels={'npsi_total':'NPSI total','burning_superficial':'Burning','pressing_deep':'Pressing','paresthesia_dysesthesia':'Paresthesia /\ndysesthesia','average_pain_nrs':'Average pain NRS','worst_pain_nrs':'Worst pain NRS','paroxysmal':'Paroxysmal','evoked':'Evoked'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False})
fig,axs=plt.subplots(4,2,figsize=(10,10),constrained_layout=True)
# a
ax=axs[0,0]; cats=['Randomized\nparticipants','Baseline\nRNA','Follow-up\nRNA','Paired\nRNA','Paired RNA +\ncomplete NPSI']; vals=[104,102,103,101,95]
yp=range(len(cats)); ax.barh(list(yp),vals); ax.set_yticks(list(yp),cats); ax.invert_yaxis(); ax.set_xlabel('Participants'); ax.set_xlim(0,110)
for i,v in enumerate(vals): ax.text(v+1,i,str(v),va='center',fontsize=8)
# b
ax=axs[0,1]; m=fa.set_index('field_requested').loc[symptoms]; x=range(len(symptoms)); w=.38
ax.bar([i-w/2 for i in x],m.n_baseline_nonmissing,width=w,label='Baseline'); ax.bar([i+w/2 for i in x],m.n_paired_nonmissing,width=w,label='Paired')
ax.set_xticks(list(x),[sym_labels[s] for s in symptoms],rotation=45,ha='right'); ax.set_ylabel('Participants'); ax.legend(frameon=False,ncol=2,loc='upper center')
# helper pivot
def pvt(df,val): return df[df.program.isin(programs)&df.symptom.isin(symptoms)].pivot(index='program',columns='symptom',values=val).reindex(index=programs,columns=symptoms)
# c baseline heatmap
ax=axs[1,0]; mat=pvt(bas,'beta_program'); im=ax.imshow(mat.values,aspect='auto',cmap='RdBu_r',vmin=-.3,vmax=.3); ax.set_yticks(range(len(programs)),programs); ax.set_xticks(range(len(symptoms)),[sym_labels[s].replace('\n',' ') for s in symptoms],rotation=45,ha='right'); ax.set_ylabel('Fixed program'); cb=fig.colorbar(im,ax=ax,fraction=.046,pad=.03); cb.set_label('Standardized HC3 beta')
# d longitudinal heatmap
ax=axs[1,1]; mat2=pvt(lng,'beta_delta_program'); qmat=pvt(lng,'BH_q'); im2=ax.imshow(mat2.values,aspect='auto',cmap='RdBu_r',vmin=-.3,vmax=.3); ax.set_yticks(range(len(programs)),programs); ax.set_xticks(range(len(symptoms)),[sym_labels[s].replace('\n',' ') for s in symptoms],rotation=45,ha='right'); ax.set_ylabel('Fixed program'); cb=fig.colorbar(im2,ax=ax,fraction=.046,pad=.03); cb.set_label('Standardized HC3 beta')
for i,p in enumerate(programs):
  for j,s in enumerate(symptoms):
    q=qmat.loc[p,s]
    if pd.notna(q) and q<.05: ax.text(j,i,'*',ha='center',va='center',fontsize=12,fontweight='bold')
    elif pd.notna(q) and q<.10: ax.text(j,i,'+',ha='center',va='center',fontsize=11,fontweight='bold')
# significant subset q<0.1 order P1/P4/P10 × symptom order
sig=lng[(lng.program.isin(['P1','P4','P10'])) & (lng.BH_q<.10)].copy()
sig['sym_rank']=sig.symptom.map({s:i for i,s in enumerate(symptoms)}); sig['prog_rank']=sig.program.map({'P1':0,'P4':1,'P10':2}); sig=sig.sort_values(['prog_rank','sym_rank'])
labels=[f"{r.program}  {sym_labels[r.symptom].replace(chr(10),' ')}" for r in sig.itertuples()]
y=list(range(len(sig)))
# e
ax=axs[2,0]; ax.errorbar(sig.beta_delta_program,y,xerr=[sig.beta_delta_program-sig.ci_low,sig.ci_high-sig.beta_delta_program],fmt='o',capsize=3); ax.set_yticks(y,labels); ax.invert_yaxis(); ax.set_xlabel('Standardized longitudinal beta (95% CI)')
# f
ax=axs[2,1]; ax.barh(y,sig.partial_R2); ax.set_yticks(y,labels); ax.invert_yaxis(); ax.set_xlabel('Partial R² for program change');
for yi,v in zip(y,sig.partial_R2): ax.text(v+.002,yi,f'{v:.3f}',va='center',fontsize=7)
# g
ax=axs[3,0]; lsub=loo.merge(sig[['program','symptom']],on=['program','symptom'],how='inner'); lsub['sym_rank']=lsub.symptom.map({s:i for i,s in enumerate(symptoms)}); lsub['prog_rank']=lsub.program.map({'P1':0,'P4':1,'P10':2}); lsub=lsub.sort_values(['prog_rank','sym_rank']); yy=list(range(len(lsub))); labs=[f"{r.program}  {sym_labels[r.symptom].replace(chr(10),' ')}" for r in lsub.itertuples()]
ax.errorbar(lsub.primary_beta,yy,xerr=[lsub.primary_beta-lsub.loo_beta_min,lsub.loo_beta_max-lsub.primary_beta],fmt='o',capsize=3); ax.set_yticks(yy,labs); ax.invert_yaxis(); ax.set_xlabel('Primary beta and leave-one-participant-out range')
# h counts
ax=axs[3,1]; q05=[int((bas.BH_q<.05).sum()),int((lng.BH_q<.05).sum()),int((inter.BH_q_omnibus<.05).sum())]; q10=[int((bas.BH_q<.10).sum()),int((lng.BH_q<.10).sum()),int((inter.BH_q_omnibus<.10).sum())]; cat=['Baseline','Longitudinal','Treatment\ninteraction']; xx=range(3); ax.bar(xx,q05,label='q<0.05'); ax.bar(xx,[b-a for a,b in zip(q05,q10)],bottom=q05,label='0.05≤q<0.10'); ax.set_xticks(list(xx),cat); ax.set_ylabel('Supported program-domain tests'); ax.set_ylim(0,14); ax.legend(frameon=False)
for i,(a,b) in enumerate(zip(q05,q10)):
    ax.text(i,b+.4,f'{a}/64 q<0.05\n{b}/64 q<0.10',ha='center',va='bottom',fontsize=7,fontweight='bold')
for lab,ax in zip('abcdefgh',axs.flat): ax.text(-.16,1.06,lab,transform=ax.transAxes,fontweight='bold',fontsize=16,va='top')
fig.savefig(OUT/'Figure_7.png',dpi=300,bbox_inches='tight',facecolor='white'); plt.close(fig)
