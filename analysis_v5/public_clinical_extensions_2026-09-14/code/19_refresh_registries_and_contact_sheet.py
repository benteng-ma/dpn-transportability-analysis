from __future__ import annotations

import hashlib
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from PIL import Image, ImageOps, ImageDraw

RUN = Path(__file__).resolve().parents[1]
INTEGRATION = RUN / "06_integration"

source_map = {
    "M01_F1_structure.png": "01_GSE302658/results/GSE302658_CLINICAL_ANALYSIS_WIDE.tsv",
    "M01_F2_baseline.png": "01_GSE302658/results/BASELINE_SYMPTOM_ASSOCIATIONS.tsv",
    "M01_F3_longitudinal.png": "01_GSE302658/results/LONGITUDINAL_SYMPTOM_ASSOCIATIONS.tsv",
    "M01_F4_interaction.png": "01_GSE302658/results/TREATMENT_INTERACTION_OMNIBUS.tsv",
    "M02_F1_cohorts.png": "02_GSE286347/results/M02_PUBLIC_METADATA_AUDIT.tsv",
    "M02_F2_MODEL_EFFECT_COMPARISON.png": "02_GSE286347/results/M02_EWAS_MODEL_COMPARISON.tsv",
    "M02_F3_fixed_regions.png": "02_GSE286347/results/M02_FIXED_REGION_EFFECT_COMPARISON.tsv",
    "M02_F4_burden.png": "02_GSE286347/results/M02_PROPENG_BURDEN_ASSOCIATIONS.tsv",
    "M03_F1_groups.png": "03_GSE14806x/results/M03_RRBS_SAMPLE_AUDIT.tsv",
    "M03_F2_RNA_PROMOTER_EFFECTS.png": "03_GSE14806x/results/M03_CROSSOMIC_GENE_EFFECTS.tsv.gz",
    "M03_F3_PROGRAM_EVIDENCE_MATRIX.png": "03_GSE14806x/results/M03_PROGRAM_CROSSOMIC_SUMMARY.tsv",
    "M03_F4_GLOBAL_SAMPLE_QC.png": "03_GSE14806x/results/M03_RRBS_GLOBAL_SAMPLE_QC.tsv",
    "M03_DIAGNOSTIC_BETA_VS_M.png": "03_GSE14806x/results/M03_RRBS_DIFFERENTIAL_RESULTS.tsv.gz",
    "M04_F1_PUBLIC_FIELD_AVAILABILITY.png": "04_JCI184075/results/M04_DIRECT_KEY_AUDIT.tsv",
    "M05_F1_real_coordinates.png": "05_GSE295206/results/M05_NAGEOTTE_BURDEN_BY_SECTION.tsv",
    "M05_F2_burden.png": "05_GSE295206/results/M05_NAGEOTTE_BURDEN_BY_DONOR.tsv",
    "M05_F3_program_burden.png": "05_GSE295206/results/M05_NAGEOTTE_BURDEN_PROGRAM_ASSOCIATIONS.tsv",
    "M05_F4_ARM_SPECIFIC_BURDEN.png": "05_GSE295206/results/M05_NAGEOTTE_BURDEN_ARM_ASSOCIATIONS.tsv",
    "M05_F5_MARKER_PROXY_ROI_DIFFERENCES.png": "05_GSE295206/results/M05_SPATIAL_MARKER_PROXY_SUMMARY.tsv",
    "M05_F6_TOPOLOGICAL_RING_TRENDS.png": "05_GSE295206/results/M05_TOPOLOGICAL_RING_ANALYSIS.tsv",
}

figures = []
for module_dir in sorted(RUN.glob("0[1-5]_*/figures")):
    for figure in sorted(module_dir.glob("*.png")):
        rel = figure.relative_to(RUN).as_posix()
        module = figure.name[1:3]
        src = source_map.get(figure.name, "")
        figures.append(
            {
                "figure": figure.name,
                "figure_path": rel,
                "source_tables": src,
                "source_exists": bool(src) and (RUN / src).exists(),
                "display_role": "diagnostic" if "DIAGNOSTIC" in figure.name else "formal",
                "data_status": "actual public-data result",
                "sha256": hashlib.sha256(figure.read_bytes()).hexdigest(),
                "module": module,
            }
        )

df = pd.DataFrame(figures)
df[["figure", "figure_path", "source_tables", "source_exists", "display_role", "data_status", "sha256"]].to_csv(
    INTEGRATION / "FIGURE_SOURCE_REGISTER.tsv", sep="\t", index=False
)
df[["figure", "module", "source_tables", "source_exists", "display_role", "data_status"]].rename(
    columns={"source_tables": "source_record"}
).to_csv(INTEGRATION / "PANEL_SOURCE_REGISTER.tsv", sep="\t", index=False)

formal = df[df["display_role"] == "formal"].copy()
thumb_w, thumb_h = 480, 330
header = 38
cols = 3
rows = (len(formal) + cols - 1) // cols
sheet = Image.new("RGB", (cols * thumb_w, rows * (thumb_h + header)), "white")
draw = ImageDraw.Draw(sheet)
for idx, row in formal.reset_index(drop=True).iterrows():
    im = Image.open(RUN / row["figure_path"]).convert("RGB")
    im.thumbnail((thumb_w - 16, thumb_h - 16))
    canvas = Image.new("RGB", (thumb_w, thumb_h), "white")
    x = (thumb_w - im.width) // 2
    y = (thumb_h - im.height) // 2
    canvas.paste(im, (x, y))
    canvas = ImageOps.expand(canvas, border=1, fill="#b0b0b0")
    col = idx % cols
    rr = idx // cols
    px = col * thumb_w
    py = rr * (thumb_h + header)
    sheet.paste(canvas, (px, py + header))
    draw.text((px + 8, py + 10), row["figure"], fill="black")

contact = RUN / "logs" / "FIGURE_CONTACT_SHEET_FINAL.png"
sheet.save(contact, dpi=(150, 150))

summary = {
    "figure_count": len(df),
    "formal_figure_count": len(formal),
    "diagnostic_figure_count": int((df["display_role"] == "diagnostic").sum()),
    "all_sources_exist": bool(df["source_exists"].all()),
    "contact_sheet": contact.relative_to(RUN).as_posix(),
}
(RUN / "tests" / "FIGURE_REGISTRY_STATUS.txt").write_text(
    "\n".join(f"{k}={v}" for k, v in summary.items()) + "\n", encoding="utf-8"
)
print(summary)
