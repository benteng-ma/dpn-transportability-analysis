import os
from pathlib import Path
import json
import pandas as pd
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import mm
from reportlab.pdfbase.cidfonts import UnicodeCIDFont
from reportlab.pdfbase import pdfmetrics
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,PageBreak,Image,KeepTogether

P = Path(os.environ.get("DPN_PROJECT_ROOT", Path(__file__).resolve().parents[3]))
C=P/"analysis_v4/closeout_20260912T062015Z"; OUT=C/"10_readout/DPN_v4_CLOSEOUT_REVIEW_REPORT.pdf"; OUT.parent.mkdir(exist_ok=True)
pdfmetrics.registerFont(UnicodeCIDFont("STSong-Light"))
styles=getSampleStyleSheet();
body=ParagraphStyle("cn",parent=styles["BodyText"],fontName="STSong-Light",fontSize=9.2,leading=14,spaceAfter=5)
h1=ParagraphStyle("h1cn",parent=styles["Heading1"],fontName="STSong-Light",fontSize=15,leading=20,spaceBefore=10,spaceAfter=8)
h2=ParagraphStyle("h2cn",parent=styles["Heading2"],fontName="STSong-Light",fontSize=11.5,leading=16,spaceBefore=8,spaceAfter=5)
title=ParagraphStyle("titlecn",parent=styles["Title"],fontName="STSong-Light",fontSize=19,leading=25,alignment=TA_CENTER)
small=ParagraphStyle("smallcn",parent=body,fontSize=7.5,leading=10)

def pt(x,style=body): return Paragraph(str(x).replace("&","&amp;"),style)
def tri(df,widths=None,fontsize=6.5):
    data=[[pt(c,small) for c in df.columns]]+[[pt("NE" if pd.isna(v) else (f"{v:.4g}" if isinstance(v,float) else v),small) for v in row] for row in df.itertuples(index=False,name=None)]
    t=Table(data,colWidths=widths,repeatRows=1,hAlign="LEFT")
    t.setStyle(TableStyle([("FONTNAME",(0,0),(-1,-1),"STSong-Light"),("FONTSIZE",(0,0),(-1,-1),fontsize),("VALIGN",(0,0),(-1,-1),"TOP"),("LINEABOVE",(0,0),(-1,0),.8,colors.black),("LINEBELOW",(0,0),(-1,0),.8,colors.black),("LINEBELOW",(0,-1),(-1,-1),.8,colors.black),("LEFTPADDING",(0,0),(-1,-1),2),("RIGHTPADDING",(0,0),(-1,-1),2),("TOPPADDING",(0,0),(-1,-1),2),("BOTTOMPADDING",(0,0),(-1,-1),2)]))
    return t

def footer(canvas,doc):
    canvas.saveState(); canvas.setFont("STSong-Light",8); canvas.drawCentredString(A4[0]/2,12*mm,f"DPN v4 定向收尾审阅报告   第 {doc.page} 页"); canvas.restoreState()

story=[Paragraph("DPN v4 定向收尾审阅报告",title),Paragraph("phase0_6_human_dpn_stage_projection · 2026-09-12",ParagraphStyle("sub",parent=body,alignment=TA_CENTER)),Spacer(1,8)]
story += [Paragraph("结论摘要",h1),pt("本轮没有重启全部 v4，也没有扩大显著性搜索。C01 证实 M04 注释缓存来自固定 NCBI GeneInfo 快照，修正的是报告分母而不是表达矩阵；C02 在当前 R 环境中精确重现 167 个德国和 660 个英国 DMR 坐标；C03 用全部 167 个德国固定区域在 92 名英国患者中作已知结果后的事后复测，没有区域达到 BH q<0.10；C05 完成 6 名供者的真实空间 ROI、两臂、逐基因贡献和留一供者分析，7/9 个可评价程序 q=0.04018，但供者只有糖尿病史而非明确 DPN 诊断，且缺少组织图像、尺度和病灶边界。P9 全 RNA 评分继续不可评价。FinnGen 与 ALC/EDC 仍受权限或正式患者键阻断。"),pt("这些结果支持一条更克制也更清晰的叙事：修复后的程序在不同组织和病理空间中呈现方向依赖、成分依赖和背景依赖的对应关系；它们尚不支持临床预测、因果机制、分解优越性或独立 DPN 细胞验证。")]

status=pd.read_csv(C/"00_admin/CLOSEOUT_STATUS.tsv",sep="\t")
story += [Paragraph("C01–C08 实际状态",h1),tri(status,[12*mm,32*mm,23*mm,107*mm])]
story += [Paragraph("C01：M04 两级流与注释缓存",h1)]
flow=pd.read_csv(C/"01_M04_flow/helper_flow_validation_retry/COUNT_FLOW.tsv",sep="\t")
story += [tri(flow,[35*mm,55*mm,35*mm]),Spacer(1,5),pt("固定 NCBI GeneInfo 快照重建得到的 38,551 个 Ensembl–GeneID 对与历史缓存完全一致。原报告把不同层级计数相减导致口径混乱；正确流为 58,174 个原始 feature、35,749 个唯一映射原始行、27,506 个合格聚合前行和 27,506 个最终排名 GeneID。没有证据要求重算 M04 评分或全基因差异分析。"),Image(str(C/"07_figures/Supplementary_Figure_S11_M04_flow.png"),width=150*mm,height=81*mm)]

story += [PageBreak(),Paragraph("C02–C04：甲基化方法审计与固定区域复测",h1)]
rep=pd.read_csv(C/"02_M02_audit/DMR_exact_reproduction_v2/DMR_REPRODUCTION_AUDIT.tsv",sep="\t")
story += [tri(rep,[28*mm,25*mm,22*mm,24*mm,30*mm,42*mm]),pt("实际调用为 dmrcate(CpGannotated, lambda=1000, C=2, pcutoff=0.05, min.cpgs=3)。输入对象由完整 EWAS 的 t、P、beta 差、BH q 和 is.sig=(q<0.05) 重建。两队列均无单 CpG q<0.05，但这本身不构成 DMR 错误；独立重跑精确复现坐标和数值。"),Image(str(C/"07_figures/Supplementary_Figure_S10_fixed_regions.png"),width=170*mm,height=71*mm)]
fx=pd.read_csv(C/"03_M02_fixed_regions/PROPENG_retest_v2/FIXED_REGION_RETEST_ALL.tsv",sep="\t")
summary=pd.DataFrame([{"固定区域":len(fx),"可评价":int((fx.status=="COMPLETED").sum()),"名义P<0.05":int((fx.P_two_sided<.05).sum()),"BH q<0.10":int((fx.BH_full_discovery_family<.10).sum()),"最小q":fx.BH_full_discovery_family.min(),"方向一致":int(fx.direction_matches_discovery.fillna(False).sum())}])
story += [tri(summary),pt("复测使用全部 167 个德国区域作为校正家族，英国矩阵只保留其自身 QC 合格探针，不插补。该检验是区域平均 M 值的性别调整模型，不等同于 DMRcate 平滑检验；两队列结果事前均已知，因此不能表述为新盲法验证。年龄和可靠技术批次仍缺失，主要模型保持受阻。细胞组成、临床联接以及未实际执行的分臂/启动子敏感性没有被伪装成已完成。")]

story += [PageBreak(),Paragraph("C05：供者级空间 ROI 证据",h1),Image(str(C/"07_figures/Figure_7_spatial_ROI.png"),width=175*mm,height=132*mm)]
sp=pd.read_csv(C/"04_M03_spatial/helper_spatial_oracle/SPATIAL_AUDIT.tsv",sep="\t")
story += [tri(sp[["program","n_donors","effect","P_two_sided","BH_q_eligible_program_family","lodo_min_effect","lodo_max_effect","status"]],[16*mm,16*mm,22*mm,20*mm,25*mm,25*mm,25*mm,28*mm]),pt("所有 16 张切片先形成切片内配对差，再在 DRG 内和供者内等权。P9 因 source-down 仅 7 个成员保持 NE。空间分数、两臂和 129,537 条逐切片成员贡献均来自真实 H5、固定成员和作者 ROI 标签。80 个 tar 成员中未发现组织学图像、scalefactors_json、验证过的 spot 直径或病灶边界，因此距离梯度保持 NE。")]

story += [PageBreak(),Paragraph("C06–C08：模块删除、全基因结果与访问边界",h1),Image(str(C/"07_figures/Supplementary_Figure_S12_block_deletions.png"),width=170*mm,height=89*mm)]
bs=pd.read_csv(C/"05_block_summary/helper_block_summary/BLOCK_SUMMARY.tsv",sep="\t")
story += [pt(f"固定模块删除共有 632 条可评价记录，覆盖 {bs.shape[0]} 个来源–终点–程序组合。完整表保留删除前后效应、剩余上下臂成员、方向保持和 NE；不把‘可评价’自动写成‘稳健’，也不按删除后较小 P 值挑块。"),pt("M04 的 primary 与 sex-adjusted DE 使用不同设计矩阵，filterByExpr 可产生不同基因宇宙，因此两者作为独立敏感性家族报告。M05 的 AH 对比有 1,288 个 q<0.05 基因，MG 和诊断×肌肉交互为 0；AH 显著不能推出交互成立。缺少真实疼痛或病理键的子任务继续 NE。FinnGen 和 ALC/EDC 没有授权输入或正式患者映射；未发送邮件、未代填协议、未猜键。")]

story += [Paragraph("代码包与验证边界",h1),pt("用户下载目录中只找到 10 个散文件，没有完整 DPN_v4_Closeout_CodeKit 文件夹。缺失 seal_bundle.py、tests/、smoke_fixed_region.R、docx_three_line.py 和 REPORTED_BASELINE.md，因此不能声称原 CodeKit 封印校验或 24 项本地单元测试通过。两个 R 脚本均完成解析；运行时和 CpGannotated 对象已审计；helper 的 flow、blocks、spatial 和 tar-inventory 已在真实适配输入上执行。缺失的 baseline 只生成了当日重建说明，未倒填。"),Paragraph("稿件整合",h1),pt("作者审阅稿在 v3.1 严格背景稿基础上加入甲基化和空间收尾证据，标题不含连字符；Figure 7 展示供者分布、两臂、真实坐标和留一供者范围。合并补充材料加入 S10–S13 及 S17–S24。正文和补表均为黑色文字及三线表。历史 GitHub/Zenodo 仅按旧版本引用，本轮未发布、未更新 DOI、未投稿。")]

doc=SimpleDocTemplate(str(OUT),pagesize=A4,rightMargin=18*mm,leftMargin=18*mm,topMargin=17*mm,bottomMargin=20*mm,title="DPN v4 Closeout Review Report")
doc.build(story,onFirstPage=footer,onLaterPages=footer)
print(OUT)
