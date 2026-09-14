from common import *
import pandas as pd,numpy as np
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,PageBreak,Table,TableStyle,Image
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet,ParagraphStyle
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import A4
from PIL import Image as PILImage

LABELS={'late_shared_concordant_neuronal_core':'Late core','late_neuron_residual':'Neuron residual','late_allcell_residual':'All-cell residual','severity_neuron_shared_concordant_core':'Severity core','severity_neuron_residual':'Severity residual','original_early_allcell':'Early parent','original_late_allcell':'All-cell parent','original_late_neuron':'Neuron parent','original_severity':'Severity parent','original_xenium':'Xenium parent'}
ENDS={'original_hDRG_DPN_control':'Human hDRG: DPN vs control','original_JCI_sural_DPN_control':'Sural nerve: DPN vs control','original_JCI_sural_severe_moderate':'Sural nerve: severe vs moderate','GSE24290_C01':'GSE24290: progression grouping','GSE148059_C02':'GSE148059: pathology grouping'}
def read(f):
 override=B/'strict_background_addendum/consolidated_results'/Path(f).name
 return pd.read_csv(override if str(f).startswith('results/') and override.exists() else B/f,sep='\t')
def fmt(x):
 if pd.isna(x):return 'NE'
 if isinstance(x,str):return x
 if isinstance(x,(int,np.integer)):return str(x)
 return f'{x:.4g}'
def main():
 styles=getSampleStyleSheet()
 for s in styles.byName.values():s.textColor=colors.black
 styles.add(ParagraphStyle(name='Text',fontName='Helvetica',fontSize=10,leading=14,spaceAfter=9))
 styles.add(ParagraphStyle(name='SmallTable',fontName='Helvetica',fontSize=8,leading=10))
 story=[];W=A4[0]-88
 def p(t):story.append(Paragraph(t,styles['Text']))
 def heading(t):story.append(Paragraph(t,styles['Heading1']));story.append(Spacer(1,5))
 def page(t):
  if story:story.append(PageBreak())
  heading(t)
 def table(headers,rows,widths=None):
  data=[[Paragraph(str(c),styles['SmallTable']) for c in headers]]+[[Paragraph(fmt(c),styles['SmallTable']) for c in row] for row in rows]
  t=Table(data,colWidths=widths or [W/len(headers)]*len(headers),repeatRows=1,hAlign='LEFT')
  t.setStyle(TableStyle([('LINEABOVE',(0,0),(-1,0),1,colors.black),('LINEBELOW',(0,0),(-1,0),.6,colors.black),('LINEBELOW',(0,-1),(-1,-1),1,colors.black),('VALIGN',(0,0),(-1,-1),'TOP'),('LEFTPADDING',(0,0),(-1,-1),4),('RIGHTPADDING',(0,0),(-1,-1),4),('TOPPADDING',(0,0),(-1,-1),5),('BOTTOMPADDING',(0,0),(-1,-1),5)]));story.append(t);story.append(Spacer(1,12))
 def fig(name,caption):
  path=B/'figures'/name;im=PILImage.open(path);w,h=im.size;scale=min(W/w,410/h);story.append(Image(str(path),width=w*scale,height=h*scale));story.append(Spacer(1,12));p(caption)
 page('DPN gene-definition repair and impact assessment')
 p('Version 1, strict-background addendum | 11 September 2026 | Local review only')
 p('This is retrospective correction after known v2 and Batch01-04 results, not preregistration, a new discovery cohort, or a submission-ready manuscript. Original source files and frozen outputs remain unchanged. No repository or DOI was published.')
 heading('Decisions supported by the actual results')
 p('The repaired severity core has 16 upregulated and 7 downregulated GeneIDs. It fails the unchanged ten-per-arm gate. Its former axonal-loss association cannot be carried forward as a corrected result; this is not evidence that the biological association is absent.')
 p('Some late-program associations remain at the historical exploratory q&lt;0.10 threshold. They do not establish strong external validation. The complete corrected empirical decomposition benchmark gives no Holm&lt;0.05 advantage; the smallest adjusted P is 0.08658.')
 table(['Historical donor test','Old P / q','Repaired stat-free P / q'],[
 ['hDRG late core','0.02667 / 0.08000','0.03333 / 0.08000'],
 ['Sural DPN neuron residual','0.01407 / 0.02814','0.02814 / 0.05628'],
 ['Sural severity core','0.002521 / 0.005042','NE: seven down GeneIDs'],
 ['Sural severity residual','0.15546 / 0.15546','0.14580 / 0.14580']],[205,140,W-345])
 p('Values above use the original one-sided tests and original eligible correction families. The new common-method two-sided results are displayed separately. Loss of an eligible program also changes the finite-test BH denominator; this is not independent evidence.')
 p('ALC pathology, fixed-marker-adjusted models and corrected TF exclusion sensitivity remain unsupported at q&lt;0.10. Three CIAP reference contrasts have q=0.07033 at each of two cell-count thresholds; none has q&lt;0.05 and covariate sensitivity is unsupported. CIAP findings are not DPN validation.')
 page('1. Source repair and locked definitions')
 p('All 6,038 original rows from five source sheets were reread, retaining source comparison, Excel row, raw symbol, signed effect and adjusted P. The archived human NCBI snapshot was reused. No target effect or P determined membership. Exact duplicated source rows: zero.')
 p('Resolution priority was unique official symbol, unique case-folded official symbol, then unique synonym. Ambiguous/unmapped rows are not guessed. All deposited rows for a GeneID within a contrast must have concordant nonzero signs and satisfy unchanged q&lt;0.05 and |log2FC|&gt;0.585. Median source effect is annotation, not a new test.')
 p('There are 4,700 eligible gene-contrast records, 892 threshold/finiteness disagreements and 156 within-source direction conflicts. These are records across contrasts, not independent genes or patients. The deposited source lists are not a complete tested transcriptome.')
 counts=read('results/REPAIRED_PROGRAM_COUNTS.tsv');table(['Repaired program','Up','Down','Total'],[[LABELS[r.module_id],r.up,r.down,r.total] for r in counts.itertuples()],[260,75,75,W-410])
 p('The actual rules/members/input lock preceded corrected disease scoring: 2026-09-10 12:53:50 UTC. Later historical adapters and reporting are post-lock implementation, not retrospectively preregistered work. Canonical members SHA256: 4a1a4a604366bab914ec04b98a06ec78863417e943b23e39fea913288308ce9c.')
 page('2. Parent/core/residual/excluded identities')
 fig('Figure_R1_partition.png','Figure R1. Unique GeneIDs in the three repaired signed parents. Core and residual are disjoint; opposed source directions and unresolved partners remain explicit excluded bins. Parent directions are retained, not chosen to fit target results.')
 p('Late all-cell: 3,102 = 128 core + 2,906 residual + 27 opposed + 41 unresolved. Late neuron: 443 = 128 + 285 + 27 + 3. Severity: 455 = 23 + 337 + 22 + 73. Within-source conflicts do not enter signed parents.')
 p('All 30 measured parent identities pass. Score reconstruction uses separate up/down weights for each bin, including excluded members. A simple unweighted sum of component scores is not the full-program score. Severity structural identity passes, but the separate score-coverage gate fails.')
 page('3. Outcome-statistic-free background sensitivity')
 p('The three complete branches are legacy members/inherited background, repaired members/inherited background, and repaired members/stat-free background. The second isolates membership on existing matrices; the third also changes outcome-dependent eligibility and expression duplicate handling.')
 inv=read('inputs/MATRIX_INVENTORY.tsv');strict=pd.read_csv(B/'strict_background_addendum/inputs/MATRIX_INVENTORY.tsv',sep='\t');inv=pd.concat([inv[~((inv.background=='statfree')&inv.endpoint.isin(strict.endpoint))],strict]);rows=[]
 for eid,f in inv.groupby('endpoint',sort=False):
  z=f.set_index('background');rows.append([ENDS[eid],int(z.loc['inherited','n_genes']),int(z.loc['statfree','n_genes']),int(z.iloc[0].n_samples)])
 table(['Endpoint','Old genes','Stat-free genes','Samples'],rows,[250,85,85,W-420])
 p('Historical hDRG/JCI DE eligibility required finite target statistics and positive baseMean. DE alias rows were prioritized by finite statistic, smallest adjusted P and largest absolute statistic. hDRG expression duplicates themselves used completeness and source-row order, not target P. JCI expression and DE information were coupled within the deposited sheets.')
 p('The new background uses only deposited expression columns, resolvable GeneIDs and finite values, with per-sample median aggregation. All-zero nonnegative abundance rows are removed; rlog is not converted to counts. No target label, DE statistic, P or q selects a row. The same gene count does not imply the same expression matrix.')
 p('Final strict-background counts additionally exclude three ambiguously assigned GeneIDs in four affected datasets. GSE24290 matrices remain identical. This is a joint background/aggregation sensitivity, not a pure causal estimate of statistical leakage. The deposited expression tables may already have upstream selection.')
 p('The original 1,530 scores and original historical P/q were reproduced. New scores were independently reconstructed with pandas average percentile ranks. All score rows and old/new differences are in DONOR_SCORE_COMPARISON.tsv.gz; backgrounds and row maps are retained.')
 page('3b. Implementation discrepancy found at handoff')
 p('A final resolver-parity check failed for four target matrices. The inherited helper used a last-row dictionary for official symbols HBD, MMD2 and TEC, although the archived GeneInfo contains nonunique official-symbol mappings. Source program members already excluded these ambiguities. The first stat-free background implementation did not fully implement that rule.')
 p('The first results and failed checks are retained. A separate strict_background_addendum was locked at 2026-09-10 15:59:02 UTC before its own associations, with no change to source members, threshold or direction. Three ambiguous assignments were removed from each of the four affected expression backgrounds. No favorable replacement GeneID was selected.')
 p('Only affected main stat-free scores/tests/benchmarks, the GSE148059 pathology/proxy response and program-only cell references were recalculated. Unchanged GSE24290 fits were reused. Final benchmark Holm adjustment still includes all 40 eligible tests per background, not a reduced family from the partial rerun.')
 p('All figures and primary numeric tables in this final PDF read the consolidated strict branch. Original results/ files remain the first implementation; use strict_background_addendum/consolidated_results where a replacement exists. Explicit first-implementation versus strict score/effect/P comparisons are retained in the addendum reports.')
 page('4. Main neural effect impact')
 fig('Figure_R2_neural_effects.png','Figure R2. A, hDRG DPN vs control (5/7); B, sural DPN vs control (6/6); C, sural severe vs moderate axonal loss (13/4). Grey circles: legacy; blue squares: repaired stat-free. Error bars: stratified bootstrap 95% intervals for mean differences. NE means the signed-arm gate fails. All five components are displayed, including negative directions.')
 p('Historical and new common-method tests have different statistics and correction families. Do not treat a changed test as a membership effect. Bootstrap percentile intervals are not inversions of the permutation test; interval exclusion and permutation P need not correspond exactly.')
 p('In the original hDRG parent dual-layer framework, the late-neuron parent changes from failing to meeting the historical exploratory gate (gene q 0.24798 to 0.00910; donor q 0.09333 to 0.06333). This is retained, not suppressed. It remains a retrospective sensitivity on the same small cohort and deposited gene statistics, not a new validation cohort.')
 page('5. ALC impact')
 fig('Figure_R3_ALC_effects.png','Figure R3. A, GSE24290, 35 patients (18 progressive/17 non-progressive). B, GSE148059, 77 participants, including 22 degenerators, 27 regenerators and 28 intermediate records retained in the three-level model. Grey circles: legacy; blue squares: repaired stat-free. All ten programs are shown. Error bars: bootstrap 95% intervals.')
 p('The two ALC data products are not independent validation cohorts. Effects concern pathology classifications and their documented tissue timing, not prospective clinical prediction. Batch01 original 100,000-permutation/5,000-bootstrap pathology and ordinal analyses are separately rerun in BATCH01_PATHOLOGY_OLD_NEW.tsv.')
 page('6. Empirical decomposition benchmark')
 fig('Figure_R4_benchmarks.png','Figure R4. Entire stat-free late-family benchmark display. A, core-minus-residual group effect and 95% interval. B, observed component max-|t| percentile among the 999 expression/MAD/direction-matched random partitions for that endpoint/family. The percentile is conditional on these data and is not a validation P value.')
 p('Eighty late-family tests completed across two backgrounds; 40 severity-family tasks remain NE because the core arm gate fails. No Holm-adjusted test is below 0.05. Joint maxT accounts for choosing either component; all eligible family/endpoint/test rows enter Holm within background. This does not establish split superiority or equivalence to the whole program.')
 p('All 19,980 actual random partitions, GeneID memberships and checksums are retained. Matching ignores outcome labels when constructing strata. Every stored partition was checked for signed stratum counts and identity hashes. No favorable split was selected. No new simulations were run.')
 ass=read('results/ASSOCIATION_IMPACT_ALL.tsv');cv=read('results/COVERAGE.tsv')
 for eid in ass.endpoint.unique():
  page('7. Complete primary comparison: '+ENDS[eid])
  p('Common-method two-sided donor permutation P and BH q over eligible programs; mean-score difference and stratified percentile 95% CI. O = legacy/inherited; R = repaired/stat-free. Repaired/inherited branch is fully retained in the accompanying TSV. These are not the historical one-sided manuscript P/q.')
  rows=[]
  for mod in LABELS:
   for ver,abbr in [('legacy_inherited','O'),('repaired_statfree','R')]:
    z=ass[(ass.endpoint==eid)&(ass.module_id==mod)&(ass.version==ver)].iloc[0];c=cv[(cv.endpoint==eid)&(cv.module_id==mod)&(cv.version==ver)].iloc[0]
    interval='NE' if z.status!='COMPLETED' else f'{z.mean_difference:.4f} [{z.ci_low:.4f}, {z.ci_high:.4f}]'
    rows.append([LABELS[mod],abbr,f'{c.n_up}/{c.n_down}',interval,fmt(z.p),fmt(z.q_BH_common_method)])
  table(['Program','Ver.','Up/down','Difference [95% CI]','P','BH q'],rows,[105,28,52,190,65,W-440])
 be=read('results/REAL_BENCHMARK_ALL.tsv')
 for eid in be.endpoint.unique():
  page('8. Complete late benchmark: '+ENDS[eid])
  rows=[]
  for z in be[(be.endpoint==eid)&(be.status=='COMPLETED')].itertuples():
   test={'parent_full':'Whole','parent_restricted':'Core+residual union','split_joint_maxT':'Joint component maxT','core_minus_residual':'Core-minus-residual'}[z.test]
   rows.append([z.background,z.family.replace('late_',''),test,z.p,z.holm_p])
  table(['Background','Parent','Test','P','Holm P'],rows,[85,70,190,80,W-425])
  p('Four severity tests for each of the two backgrounds are retained as NE in the complete TSV. No severity test is silently omitted or replaced by a lower coverage threshold. Whole/restricted program tests are association tests, not superiority tests; core-minus-residual and the matched partition reference do not by themselves establish clinical added value.')
 page('9. Secondary dependencies and negative results')
 p('Historical neural gene and sample layers, JCI TPM scores, rat DRG and ocular parent/component projections, PBMC projections, GSE302658 fixed clinical analyses, five protein comparisons, paired nerve scores, reference localization and proxy models have separate complete outputs. Original patient/animal/sample keys and numerical interfaces are reused; no new clinical or RNA-RRBS links are inferred.')
 table(['Dependent output','Impact / boundary'],[
 ['ALC marker-adjustment models','72 finite q among 80 rows; no q<0.10. Marker expression proxies are not cell fractions.'],
 ['TF program-exclusion sensitivity','449 TF group tests, no q<0.10. Original ABC/AB primary screens and prior TF universe were not rebuilt.'],
 ['Cell reference','Three CIAP vs CTRL findings at q=0.07033 at each of two thresholds; no q<0.05. Covariate sensitivity: no q<0.10.'],
 ['Source protein','All-cell residual retains matched-null q=0.000300; source tissue evidence is not independent clinical replication.'],
 ['Clinical pain','Primary baseline severity score gate remains false; longitudinal unadjusted association does not meet the original full confirmation gate.'],
 ['Spatial','Repaired ten programs: 0/10 eligible under the unchanged panel rule. No new large spatial download.']],[155,W-155])
 p('GSE285983 has only two DPN donors among 37 total; CIAP/CIDP inferences are not DPN validation. GSE168243 has five donors after merging hDRG3/hDRG5. Spatial donors overlap the reference and do not add independent DPN validation.')
 page('10. What must not be carried into the manuscript')
 p('Do not carry forward the former severity-core axonal-loss result as a corrected finding. Do not claim all former associations disappeared: weak late-program evidence and the historical parent gate changes are retained. Do not claim decomposition superiority, prospective prediction, disease-specific cell fractions or a mechanism chain.')
 p('All original 405 candidate records are preserved with canonical membership flags. The original combined statistic formula is recalculated only as a conditional historical index. Deposited statistic/shrinkage provenance remains unresolved, so its P/q-like columns do not establish formal gene validation or a new candidate list.')
 p('Source functional annotation uses the archived libraries and the original overlap-filtered testing family (overlap at least three before BH). This limitation is disclosed; the annotation is a membership-impact comparison, not a newly calibrated pathway claim. No new libraries, pathways, TFs or networks were screened.')
 p('The former official-spelling-only mapping sensitivity is not equivalent after canonicalization; it is marked not comparable, not falsely repeated as validation. Full original manuscripts, supplementary display-only PCA/distance/AUC-bootstrap panels and claim registries have not been regenerated. The bundle contains replacement main scores/effects, not authorization to reuse every old display.')
 p('Unaffected whole-gene DE, correlation-aware pathways, WGCNA/preservation, primary ABC/AB TF results, module-only cell localization and individual-gene expression were not rerun. Existing unknown clinical/RNA-RRBS linkage remains unknown. No new patient count or independence assumption was introduced.')
 page('11. Verification, provenance and handoff')
 num=json.loads((B/'reports/REPORT_NUMBERS.json').read_text())
 hc=json.loads((B/'reports/HANDOFF_CHECKS.json').read_text())
 p(f"The final consolidated view includes {num['main_association_rows']} primary association rows, {num['donor_score_rows']} donor-score rows including legacy comparators, {num['benchmark_rows']} benchmark task rows and {num['matched_partitions']} matched random partitions. The current completed test suites contain {hc['all_passes']:,} passing checks. Four earlier resolver failures are retained and resolved in the strict addendum; the first implementation is not falsely certified. Most checks are per-record or hash checks, not independent scientific validations.")
 p('Independent checks cover original score and historical P/q parity, every source eligibility decision, canonical exclusivity and source direction, weighted parent reconstruction, new score recomputation, every stored matched partition, multiple-testing calculations and protected-file hashes. Code/test success does not establish scientific correctness or sample independence.')
 p('Actual failures remain in logs: executable lookup, R/rlang binary compatibility and a post-computation parent-context join-key mismatch. The corrected R 4.6.1 run and successful keyed joins are recorded separately; no successful result was discarded to hide a failed attempt. No global dependency installation occurred.')
 p('Use STATUS.tsv, manifest hashes and actual script logs to distinguish completed analyses, unaffected dependencies and non-evaluable work. This local bundle requires the original authorized source workbooks/reference caches and archived helper code; it is not represented as a standalone clean-environment public reproduction. Legacy date suffixes are interface names, not this run date.')
 p('Original v2 and Batch01-04 remain unchanged. This report is an impact review, not a revised submission package. Nothing was uploaded to GitHub or Zenodo. A future manuscript must be created as another version and explicitly acknowledge the retrospective correction and weakened/non-evaluable claims.')
 def footer(canvas,doc):
  canvas.setFillColor(colors.black);canvas.setFont('Helvetica',8);canvas.drawCentredString(A4[0]/2,24,str(doc.page))
 dst=B/'reports/REPAIR_REVIEW_REPORT.pdf';SimpleDocTemplate(str(dst),pagesize=A4,rightMargin=44,leftMargin=44,topMargin=36,bottomMargin=38,title='DPN canonical gene-definition repair v1',author='').build(story,onFirstPage=footer,onLaterPages=footer)
 print(dst,flush=True)
if __name__=='__main__':run(main)
