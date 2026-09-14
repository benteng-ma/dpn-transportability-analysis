from __future__ import annotations

import os

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from docx import Document
from docx.enum.section import WD_SECTION
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, RGBColor

P = Path(os.environ.get("DPN_PROJECT_ROOT", Path(__file__).resolve().parents[3]))
C=P/"analysis_v4/closeout_20260912T062015Z"
T=P/"analysis_v4/targeted_repair_2026-09-12"
V3=P/"analysis_v3/manuscript_integration_v1/documents"
FIG=C/"07_figures"; DOC=C/"08_manuscript"; FIG.mkdir(exist_ok=True); DOC.mkdir(exist_ok=True)

def black_runs(doc):
    for p in doc.paragraphs:
        for r in p.runs: r.font.color.rgb=RGBColor(0,0,0)
    for table in doc.tables:
        for row in table.rows:
            for cell in row.cells:
                for p in cell.paragraphs:
                    for r in p.runs: r.font.color.rgb=RGBColor(0,0,0)

def three_line(table):
    tblPr=table._tbl.tblPr
    old=tblPr.find(qn("w:tblBorders"))
    if old is not None: tblPr.remove(old)
    borders=OxmlElement("w:tblBorders")
    for edge in ("top","bottom"):
        el=OxmlElement(f"w:{edge}"); el.set(qn("w:val"),"single"); el.set(qn("w:sz"),"8"); el.set(qn("w:color"),"000000"); borders.append(el)
    for edge in ("left","right","insideH","insideV"):
        el=OxmlElement(f"w:{edge}"); el.set(qn("w:val"),"nil"); borders.append(el)
    tblPr.append(borders)
    for cell in table.rows[0].cells:
        tcPr=cell._tc.get_or_add_tcPr(); b=tcPr.find(qn("w:tcBorders"))
        if b is None: b=OxmlElement("w:tcBorders"); tcPr.append(b)
        for edge in ("top","left","right"):
            el=OxmlElement(f"w:{edge}"); el.set(qn("w:val"),"nil"); b.append(el)
        el=OxmlElement("w:bottom"); el.set(qn("w:val"),"single"); el.set(qn("w:sz"),"8"); el.set(qn("w:color"),"000000"); b.append(el)

def add_table(doc, df, max_rows=None, font=7.5):
    d=df if max_rows is None else df.head(max_rows)
    t=doc.add_table(rows=1,cols=len(d.columns)); t.autofit=True
    for j,x in enumerate(d.columns): t.rows[0].cells[j].text=str(x)
    for row in d.itertuples(index=False,name=None):
        cells=t.add_row().cells
        for j,x in enumerate(row):
            if pd.isna(x): s="NE"
            elif isinstance(x,float): s=f"{x:.4g}"
            else: s=str(x)
            cells[j].text=s
    three_line(t)
    for row in t.rows:
        for cell in row.cells:
            for p in cell.paragraphs:
                for r in p.runs: r.font.size=Pt(font); r.font.color.rgb=RGBColor(0,0,0)
    return t

def insert_before(doc,target,text="",style=None,picture=None,width=None):
    p=doc.add_paragraph(style=style)
    if picture:
        p.alignment=WD_ALIGN_PARAGRAPH.CENTER; p.add_run().add_picture(str(picture),width=width)
    else: p.add_run(text)
    target._p.addprevious(p._p)
    return p

def make_figures():
    # Figure 7: donor distributions, arms, real coordinates and leave-one-donor-out ranges.
    don=pd.read_csv(C/"04_M03_spatial/SPATIAL_DONOR_EFFECTS.tsv",sep="\t")
    don=don[don.eligible.astype(str).str.lower().eq("true")]
    arm=pd.read_csv(C/"04_M03_spatial/M03_TWO_ARM_TESTS.tsv",sep="\t")
    oracle=pd.read_csv(C/"04_M03_spatial/helper_spatial_oracle/SPATIAL_AUDIT.tsv",sep="\t")
    coord=pd.read_csv(C/"04_M03_spatial/M03_REAL_COORDINATES_WITH_AUTHOR_ROI.tsv.gz",sep="\t")
    rep=json.loads((C/"04_M03_spatial/M03_CLOSEOUT_SUMMARY.json").read_text())["representative_section"]
    order=[f"P{i}" for i in range(1,11) if i!=9]
    fig,axes=plt.subplots(2,2,figsize=(11,8.3))
    ax=axes[0,0]
    for i,p in enumerate(order):
        x=don[don.program==p].effect.to_numpy(); ax.scatter(np.repeat(i,len(x)),x,s=24,color="#333333",alpha=.8); ax.plot([i-.22,i+.22],[x.mean(),x.mean()],color="#C44E52",lw=2)
    ax.axhline(0,color="black",lw=.7); ax.set_xticks(range(len(order)),order); ax.set_ylabel("Nageotte minus neuronal score"); ax.text(-.10,1.03,"A",transform=ax.transAxes,fontweight="bold",fontsize=12)
    ax=axes[0,1]; pv=arm.pivot(index="program",columns="arm",values="effect").reindex(order); x=np.arange(len(order))
    ax.bar(x-.18,pv.source_up,.36,color="#4C72B0",label="ΔU"); ax.bar(x+.18,-pv.source_down,.36,color="#DD8452",label="−ΔD"); ax.axhline(0,color="black",lw=.7); ax.set_xticks(x,order); ax.set_ylabel("Contribution to ΔS"); ax.legend(frameon=False,fontsize=8); ax.text(-.10,1.03,"B",transform=ax.transAxes,fontweight="bold",fontsize=12)
    ax=axes[1,0]; d=coord[coord.geo_accession==rep]; colors={"Nageotte":"#C44E52","Neuronal":"#4C72B0","mixed_excluded":"#999999"}
    for lab,g in d.groupby("roi"): ax.scatter(g.pixel_col,g.pixel_row,s=18,color=colors.get(lab,"#D9D9D9"),label=lab,alpha=.85,linewidths=0)
    ax.invert_yaxis(); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([]); ax.legend(frameon=False,fontsize=8); ax.text(-.10,1.03,"C",transform=ax.transAxes,fontweight="bold",fontsize=12)
    ax=axes[1,1]; q=oracle[oracle.status=="COMPLETED"].set_index("program").reindex(order); y=np.arange(len(q));
    ax.errorbar(q.effect,y,xerr=[q.effect-q.lodo_min_effect,q.lodo_max_effect-q.effect],fmt="o",color="#333333",ecolor="#7A8A93",capsize=2); ax.axvline(0,color="black",lw=.7); ax.set_yticks(y,order); ax.set_xlabel("Mean effect with leave-one-donor-out range"); ax.text(-.10,1.03,"D",transform=ax.transAxes,fontweight="bold",fontsize=12)
    fig.tight_layout(); fig.savefig(FIG/"Figure_7_spatial_ROI.png",dpi=450); plt.close(fig)

    # Fixed-region evidence.
    d=pd.read_csv(C/"03_M02_fixed_regions/PROPENG_retest_v2/FIXED_REGION_RETEST_ALL.tsv",sep="\t")
    fig,axes=plt.subplots(1,2,figsize=(10,4.2))
    axes[0].scatter(d.discovery_direction*d.measured_fraction,d.effect_beta_adjusted,s=16,color="#555555",alpha=.65)
    axes[0].axhline(0,color="black",lw=.7); axes[0].axvline(0,color="black",lw=.7); axes[0].set_xlabel("Discovery direction × test probe coverage"); axes[0].set_ylabel("PROPENG adjusted beta difference"); axes[0].text(-.10,1.03,"A",transform=axes[0].transAxes,fontweight="bold",fontsize=12)
    axes[1].scatter(d.effect_M,-np.log10(d.P_two_sided),s=18,color="#4C72B0",alpha=.7); axes[1].axhline(-np.log10(.05),color="#999999",lw=.8,ls="--"); axes[1].set_xlabel("PROPENG painful minus painless M difference"); axes[1].set_ylabel("−log10 two-sided P"); axes[1].text(-.10,1.03,"B",transform=axes[1].transAxes,fontweight="bold",fontsize=12)
    fig.tight_layout(); fig.savefig(FIG/"Supplementary_Figure_S10_fixed_regions.png",dpi=450); plt.close(fig)

    # M04 feature flow.
    f=pd.read_csv(C/"01_M04_flow/helper_flow_validation_retry/COUNT_FLOW.tsv",sep="\t")
    vals=dict(zip(f.stage,f["count"])); labs=["Raw features","Unique mapping","Eligible rows","Eligible GeneIDs"]; nums=[vals["all"],vals["unique_mapping"],vals["eligible_preaggregate"],vals["eligible_rank"]]
    fig,ax=plt.subplots(figsize=(7,3.8)); ax.barh(labs[::-1],nums[::-1],color=["#4C72B0","#55A868","#C44E52","#8172B2"])
    for i,n in enumerate(nums[::-1]): ax.text(n+500,i,f"{n:,}",va="center",fontsize=9)
    ax.set_xlabel("Records"); ax.spines[["top","right"]].set_visible(False); fig.tight_layout(); fig.savefig(FIG/"Supplementary_Figure_S11_M04_flow.png",dpi=450); plt.close(fig)

    # Full block-deletion descriptive distribution.
    b=pd.read_csv(C/"05_block_summary/BLOCK_INFLUENCE_STANDARDIZED.tsv",sep="\t"); b=b[b.eligible_after.astype(str).str.lower().eq("true")].copy(); b["change"]=b.deleted_effect-b.original_effect
    fig,ax=plt.subplots(figsize=(8,4.2)); groups=[b.loc[b.analysis_id==x,"change"].to_numpy() for x in ["M01","M04","M05"]]
    ax.boxplot(groups,labels=["M01 neural endpoints","M04 Morton source comparison","M05 muscle comparisons"],showfliers=False); ax.axhline(0,color="black",lw=.7); ax.set_ylabel("Deleted minus original score effect"); fig.tight_layout(); fig.savefig(FIG/"Supplementary_Figure_S12_block_deletions.png",dpi=450); plt.close(fig)
    return d

def build_manuscript():
    src=V3/"Manuscript.docx"; doc=Document(src)
    # Revised abstract: adds closeout evidence without inflating its evidential status.
    doc.paragraphs[6].text=("Sensory ganglion transcriptomic programs may differ from those measured in distal nerve. We examined five human dorsal root ganglion programs and their components across public neural, ocular, blood, methylation and spatial datasets. After auditing gene aliases and conflicting source records, we corrected definitions and reassessed dependent results; this correction occurred after the original results were known. Under common two-sided donor tests with backgrounds independent of target differential statistics, five programs had negative sural neuropathy associations at q<0.05, whereas positive neuronal associations met q<0.10 but not q<0.05. No program met q<0.10 in ganglion, axonal-loss or two acetyl-L-carnitine data-product comparisons. A retrospective test of 167 German discovery methylation regions in 92 British patients produced no BH q<0.10 region. In spatial ganglia from six donors with diabetes history, seven of nine evaluable programs differed between author-annotated Nageotte and nearby neuronal regions at q=0.0402, but explicit DPN diagnosis, tissue images and calibrated boundaries were unavailable. The previous severity-core axonal-loss finding became unevaluable because its corrected downregulated arm contained only seven genes. These results show context-dependent tissue correspondence and define limits on claims of conserved neuropathy biology, prediction or mechanism.")
    doc.paragraphs[7].text="Keywords: diabetic peripheral neuropathy; dorsal root ganglion; sural nerve; transcriptomics; gene identifiers; program scores"
    # Correct the accession typo and expand data statement.
    for p in doc.paragraphs:
        if p.text.startswith("The original public datasets"):
            p.text=p.text.replace("GSE285984","GSE295206")+" The blood methylation series was GSE286347."
    target=next(p for p in doc.paragraphs if p.text.strip()=="Discussion")
    inserts=[
        ("Blood methylation regions did not transfer across known cohorts","Heading 2",None),
        ("The GSE286347 audit retained 139 German and 92 British DPN participants after the pre-existing detection-P quality rules. The available public metadata supplied pain group and sex but not verified age or a technical batch suitable for the planned primary model. Therefore, the age- and batch-adjusted analysis remained blocked; the sex-adjusted branch was retained as exploratory [23].",None,None),
        ("Reconstruction of the installed DMRcate workflow exactly reproduced 167 German and 660 British regions from the archived complete EWAS tables, despite no single CpG meeting BH q<0.05. This establishes computational provenance, not biological validity. We then fixed all 167 German coordinates and their discovery-eligible probe memberships before testing average regional M values in the British participants. All regions were evaluable with at least three independently quality-qualified test probes. Four had nominal P<0.05, but none met BH q<0.10 (minimum q=0.3027); 70 of 167 adjusted beta effects matched the German direction (Supplementary Fig. S10; Supplementary Tables S17–S18). Because both cohort-level results were known before this retest, it is a retrospective fixed-region assessment rather than a new blinded validation.",None,None),
        ("Donor-level spatial contrasts localized programs to pathology-annotated regions","Heading 2",None),
        ("The GSE295206 archive contained 16 sections from seven DRG and six donors with diabetes history [24]. We used author-supplied Nageotte and nearby neuronal ROI labels and retained section-to-DRG-to-donor equal weighting. Nine corrected programs were evaluable; P9 remained NE because its source-down arm contained seven genes. Seven programs had exact two-sided donor sign-flip q=0.04018 for Nageotte-minus-neuronal scores. P5 and P7 were lower in Nageotte ROIs, whereas P2, P4, P6, P8 and P10 were higher (Fig. 7A). These paired spatial effects were stable in direction in leave-one-donor-out summaries, but n=6 leaves only 64 exact sign assignments.",None,None),
        ("The up- and down-arm decomposition showed that similar signed scores could arise through different expression patterns (Fig. 7B; Supplementary Table S19). For example, the positive P4 effect combined a lower source-up arm with a larger decrease in the source-down arm. Gene-level contributions were retained for every measured corrected member rather than selecting top genes by P value. Real Visium coordinates and author ROI labels are shown without a substituted tissue image (Fig. 7C; Supplementary Fig. S13). The archive lacked histology images, scalefactors_json, verified spot diameter and a lesion boundary, so a calibrated distance gradient remained NE. Donors had diabetes history, but explicit DPN diagnosis was unavailable; these findings are pathology-region localization within diabetic ganglia, not independent DPN validation.",None,None),
    ]
    for text,style,pic in inserts: insert_before(doc,target,text,style)
    insert_before(doc,target,picture=FIG/"Figure_7_spatial_ROI.png",width=Inches(6.6))
    insert_before(doc,target,"Figure 7. Donor-level spatial localization in human dorsal root ganglia with diabetes history. (A) Six donor effects for nine evaluable corrected programs; red bars indicate donor means. (B) Mean source-up contribution (ΔU) and the negative source-down contribution (−ΔD), whose sum is ΔS. The complete 18-test arm family is reported separately. (C) Real array coordinates with author ROI labels for the first quality-complete section by accession order; no histology image or calibrated boundary was available. (D) Mean effects and leave-one-donor-out ranges. P9 remained NE. These analyses localize expression relative to annotated pathology regions but do not establish DPN specificity, prospective prediction or causal mechanism.")
    # Methods additions before Evidence traceability.
    mt=next(p for p in doc.paragraphs if p.text.strip()=="Evidence traceability and availability")
    insert_before(doc,mt,"Methylation runtime audit and fixed-region retest","Heading 2")
    insert_before(doc,mt,"The GSE286347 beta matrix was linked by explicit participant codes to public pain-group and sex fields. Samples passed when no more than 1% of assayed probes had detection P>0.01 or non-finite detection P. Within each cohort, autosomal uniquely annotated probes with finite beta values in [0,1] and no more than 5% detection failures were transformed to M values after clipping at 10−6. The exploratory model contained painful status and sex. The planned age- and batch-adjusted model was not run because those covariates were unavailable or unverified. DMRcate 3.8.0 received a CpGannotated object reconstructed from the complete stored t, P, beta-difference and BH-q fields and was called with lambda=1000, C=2, pcutoff=0.05 and min.cpgs=3. Exact coordinate and numerical reproduction was required before follow-up.")
    insert_before(doc,mt,"For the known-results follow-up, all 167 German coordinates and all German QC-qualified EPIC probes within each hg19 interval were fixed before reading British regional results. British probes had to pass the independent British probe rules; no imputation was used. Each region required at least three test probes. Participant-level mean M value was fitted with limma using painful status and sex, robust empirical Bayes moderation and no trend term. Two-sided P values were BH-adjusted over all 167 selected regions, including regions that would have become NE. This average-region test is not the original smoothed DMRcate hypothesis and was explicitly recorded as a retrospective implementation.")
    insert_before(doc,mt,"Spatial ROI analysis","Heading 2")
    insert_before(doc,mt,"Filtered spatial count matrices, real array coordinates and author ROI labels were read from all 16 GSE295206 sections. Mixed spots were excluded. Gene labels were resolved against the same archived NCBI GeneInfo snapshot. Within each section, expression was summed to unique GeneID, all-zero genes were removed and percentile ranks were calculated across the eligible whole-gene background. Corrected signed scores required at least ten measured genes per arm. Nageotte-minus-neuronal differences were averaged first across sections within DRG and then across DRG within donor. Exact two-sided sign-flip tests used the donor mean and all 64 sign assignments; BH covered the nine eligible programs. The 18 arm tests were adjusted as one separate family. Leave-one-donor-out effects and 2,000-draw donor bootstraps were descriptive sensitivity analyses. No spot was treated as an independent participant.")
    # References 23–24 before Acknowledgements.
    ack=next(p for p in doc.paragraphs if p.text.strip()=="Acknowledgements")
    insert_before(doc,ack,"23. National Center for Biotechnology Information. Gene Expression Omnibus series GSE286347: whole-blood DNA methylation in painful and painless diabetic neuropathy. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE286347 (source files accessed 12 September 2026).")
    insert_before(doc,ack,"24. National Center for Biotechnology Information. Gene Expression Omnibus series GSE295206: spatial transcriptomic profiling of human dorsal root ganglia. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE295206 (source files accessed 12 September 2026).")
    for t in doc.tables: three_line(t)
    black_runs(doc)
    doc.save(DOC/"DPN_Closeout_Author_Review_Manuscript.docx")

def build_supplement():
    doc=Document(V3/"Supplementary_Information.docx")
    doc.add_page_break(); doc.add_heading("Closeout evidence supplement",level=1)
    doc.add_paragraph("The following figures and reader tables were generated in the 12 September 2026 targeted closeout. They do not overwrite the historical v3.1 supplement. Complete machine-readable records are supplied in the review bundle.")
    for label,path,caption in [
        ("Supplementary Fig. S10",FIG/"Supplementary_Figure_S10_fixed_regions.png","Retrospective fixed-coordinate assessment of 167 German discovery regions in 92 British participants. No region met BH q<0.10."),
        ("Supplementary Fig. S11",FIG/"Supplementary_Figure_S11_M04_flow.png","M04 raw-feature to canonical-GeneID flow. Counts distinguish raw rows from final rank-universe GeneIDs."),
        ("Supplementary Fig. S12",FIG/"Supplementary_Figure_S12_block_deletions.png","All 632 evaluable fixed-module deletions. Changes are descriptive and do not imply robustness merely because a row was evaluable."),
        ("Supplementary Fig. S13",FIG/"M03_all_sections_real_coordinates.png","Real array coordinates and author ROI labels for all 16 sections. These are not reconstructed histology images."),
    ]:
        doc.add_heading(label,level=2); p=doc.add_paragraph(); p.alignment=WD_ALIGN_PARAGRAPH.CENTER; p.add_run().add_picture(str(path),width=Inches(6.5)); doc.add_paragraph(f"{label}. {caption}")
    doc.add_heading("Supplementary Table S17 M04 feature and canonical-gene flow",level=2); add_table(doc,pd.read_csv(C/"01_M04_flow/helper_flow_validation_retry/COUNT_FLOW.tsv",sep="\t"))
    doc.add_heading("Supplementary Table S18 M02 DMR and fixed-region audit",level=2)
    d=pd.read_csv(C/"03_M02_fixed_regions/PROPENG_retest_v2/FIXED_REGION_RETEST_ALL.tsv",sep="\t"); add_table(doc,d[["region_id","chr","start","end","n_reference_probes","n_test_probes","effect_M","P_two_sided","BH_full_discovery_family","direction_matches_discovery"]].sort_values("P_two_sided"),max_rows=25)
    doc.add_paragraph("The 25 smallest two-sided P values are printed for readability; all 167 prespecified rows are in 03_M02_fixed_regions/PROPENG_retest_v2/FIXED_REGION_RETEST_ALL.tsv and all 167 entered BH correction.")
    doc.add_heading("Supplementary Table S19 M03 complete donor program tests",level=2); add_table(doc,pd.read_csv(C/"04_M03_spatial/helper_spatial_oracle/SPATIAL_AUDIT.tsv",sep="\t")[["program","n_donors","effect","P_two_sided","BH_q_eligible_program_family","lodo_min_effect","lodo_max_effect","status","reason"]])
    doc.add_heading("Supplementary Table S20 M03 complete two-arm family",level=2); add_table(doc,pd.read_csv(C/"04_M03_spatial/M03_TWO_ARM_TESTS.tsv",sep="\t"))
    doc.add_heading("Supplementary Table S21 Fixed module deletion summary",level=2); add_table(doc,pd.read_csv(C/"05_block_summary/helper_block_summary/BLOCK_SUMMARY.tsv",sep="\t"))
    doc.add_heading("Supplementary Table S22 M04 and M05 whole-gene interpretation",level=2); add_table(doc,pd.read_csv(C/"06_remaining_status/M04_M05_DE_INTERPRETATION.tsv",sep="\t"))
    doc.add_heading("Supplementary Table S23 Closeout task status",level=2); add_table(doc,pd.read_csv(C/"00_admin/CLOSEOUT_STATUS.tsv",sep="\t"))
    doc.add_heading("Supplementary Table S24 Evidence boundaries",level=2); add_table(doc,pd.read_csv(C/"00_admin/CLAIM_EVIDENCE_REGISTER.tsv",sep="\t"),font=7)
    for t in doc.tables: three_line(t)
    black_runs(doc)
    doc.save(DOC/"DPN_Closeout_Combined_Supplementary_Information.docx")

def build_cover():
    doc=Document(); sec=doc.sections[0]; sec.top_margin=Inches(.85); sec.bottom_margin=Inches(.85)
    p=doc.add_paragraph(); p.alignment=WD_ALIGN_PARAGRAPH.RIGHT; p.add_run("12 September 2026")
    doc.add_paragraph("Dear Editors,")
    doc.add_paragraph("Please consider our manuscript, “Sensory ganglion transcriptomic programs show divergent associations in diabetic peripheral neuropathy,” for publication as an Article in Scientific Reports.")
    doc.add_paragraph("This study tests whether directionally defined human sensory-ganglion programs retain their structure and orientation across independently processed neural, ocular, blood, methylation and spatial data. A canonical-gene audit changed component membership and made one previously emphasized severity-core result unevaluable. We report the correction explicitly and retain positive, negative and non-evaluable outcomes rather than selecting the most favorable convention.")
    doc.add_paragraph("The closeout analysis adds two bounded lines of evidence. First, an exact reconstruction of the installed DMRcate workflow confirmed the provenance of 167 German discovery regions; a retrospective fixed-coordinate test in 92 British participants found no region at BH q<0.10. Second, donor-weighted spatial analysis of 16 human dorsal-root-ganglion sections localized seven of nine evaluable corrected programs to author-annotated Nageotte versus nearby neuronal regions. This spatial result is presented with its limits: six donors, diabetes history without explicit DPN diagnosis, and no available tissue image or calibrated lesion boundary. It is therefore not described as independent DPN validation.")
    doc.add_paragraph("The principal contribution is a transparent account of how source direction, canonical identifiers, expression background and tissue context affect cross-tissue interpretation. The work does not claim clinical prediction, causal mechanism, superiority of program decomposition or independent cellular validation. Complete code, full results, failed-run logs and versioned evidence registers are included for review. The historical public repository and DOI are identified as historical and have not been updated by this closeout.")
    doc.add_paragraph("The manuscript contains no new recruitment, intervention or access to identifiable participant information. The authors declare no competing interests.")
    doc.add_paragraph("Sincerely,\nBaihua Chen, MD, PhD\nCorresponding author\nDepartment of Ophthalmology, The Second Xiangya Hospital, Central South University")
    black_runs(doc); doc.save(DOC/"DPN_Closeout_Cover_Letter_Draft.docx")

if __name__=="__main__":
    make_figures(); build_manuscript(); build_cover(); build_supplement(); print("manuscript, cover and supplement built")
