from pathlib import Path
import json, pandas as pd
from docx import Document
from docx.shared import Inches,Pt,RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.section import WD_SECTION
from docx.enum.table import WD_TABLE_ALIGNMENT,WD_CELL_VERTICAL_ALIGNMENT
from docx.oxml import OxmlElement
from docx.oxml.ns import qn

R = Path(__file__).resolve().parents[1]
OUT=R/"reports/DPN_v4_Targeted_Repair_and_Continuation_Report.docx";OUT.parent.mkdir(parents=True,exist_ok=True)
def read(rel):return pd.read_csv(R/rel,sep="\t")
def fmt(x):
 if pd.isna(x):return ""
 if isinstance(x,float):
  if x!=0 and abs(x)<0.001:return f"{x:.3e}"
  return f"{x:.4f}".rstrip('0').rstrip('.')
 return str(x)
def borders(table):
 tblPr=table._tbl.tblPr; old=tblPr.find(qn('w:tblBorders'))
 if old is not None:tblPr.remove(old)
 tb=OxmlElement('w:tblBorders')
 for edge,val,sz in [('top','single','10'),('left','nil','0'),('bottom','single','10'),('right','nil','0'),('insideH','nil','0'),('insideV','nil','0')]:
  e=OxmlElement('w:'+edge);e.set(qn('w:val'),val);e.set(qn('w:sz'),sz);e.set(qn('w:color'),'000000');tb.append(e)
 tblPr.append(tb)
 for cell in table.rows[0].cells:
  tcPr=cell._tc.get_or_add_tcPr();cb=OxmlElement('w:tcBorders');bot=OxmlElement('w:bottom');bot.set(qn('w:val'),'single');bot.set(qn('w:sz'),'8');bot.set(qn('w:color'),'000000');cb.append(bot);tcPr.append(cb)
def table(doc,headers,rows,widths=None):
 t=doc.add_table(rows=1,cols=len(headers));t.alignment=WD_TABLE_ALIGNMENT.CENTER;t.autofit=True
 for j,h in enumerate(headers):
  c=t.rows[0].cells[j];c.text=str(h);c.vertical_alignment=WD_CELL_VERTICAL_ALIGNMENT.CENTER
  for run in c.paragraphs[0].runs:run.bold=True
 for row in rows:
  cells=t.add_row().cells
  for j,x in enumerate(row):cells[j].text=fmt(x);cells[j].vertical_alignment=WD_CELL_VERTICAL_ALIGNMENT.CENTER
 for row in t.rows:
  for c in row.cells:
   for p in c.paragraphs:
    p.paragraph_format.space_after=Pt(0);p.paragraph_format.line_spacing=1.0
    for run in p.runs:run.font.size=Pt(8);run.font.color.rgb=RGBColor(0,0,0);run.font.name='Microsoft YaHei';run._element.rPr.rFonts.set(qn('w:eastAsia'),'Microsoft YaHei')
 borders(t);doc.add_paragraph().paragraph_format.space_after=Pt(1);return t
def p(doc,text,bold=False,italic=False):
 x=doc.add_paragraph();x.paragraph_format.first_line_indent=Pt(21);x.paragraph_format.space_after=Pt(5);x.paragraph_format.line_spacing=1.25
 r=x.add_run(text);r.bold=bold;r.italic=italic;return x
def head(doc,text,level=1):
 x=doc.add_heading(text,level=level);x.paragraph_format.keep_with_next=True;return x
def cap(doc,text):
 x=doc.add_paragraph();x.alignment=WD_ALIGN_PARAGRAPH.CENTER;x.paragraph_format.space_after=Pt(8);r=x.add_run(text);r.bold=True;r.font.size=Pt(9);return x
def page_no(paragraph):
 paragraph.alignment=WD_ALIGN_PARAGRAPH.CENTER
 run=paragraph.add_run();begin=OxmlElement('w:fldChar');begin.set(qn('w:fldCharType'),'begin');instr=OxmlElement('w:instrText');instr.set(qn('xml:space'),'preserve');instr.text=' PAGE ';end=OxmlElement('w:fldChar');end.set(qn('w:fldCharType'),'end');run._r.extend([begin,instr,end])

doc=Document();sec=doc.sections[0];sec.page_width=Inches(8.27);sec.page_height=Inches(11.69);sec.top_margin=Inches(.7);sec.bottom_margin=Inches(.65);sec.left_margin=Inches(.72);sec.right_margin=Inches(.72)
styles=doc.styles
for s in ['Normal','Title','Subtitle','Heading 1','Heading 2','Heading 3']:
 st=styles[s];st.font.name='Microsoft YaHei';st._element.rPr.rFonts.set(qn('w:eastAsia'),'Microsoft YaHei');st.font.color.rgb=RGBColor(0,0,0)
styles['Normal'].font.size=Pt(10);styles['Title'].font.size=Pt(22);styles['Heading 1'].font.size=Pt(15);styles['Heading 2'].font.size=Pt(12)
page_no(sec.footer.paragraphs[0])

x=doc.add_paragraph();x.alignment=WD_ALIGN_PARAGRAPH.CENTER;x.paragraph_format.space_before=Pt(80);r=x.add_run('DPN v4 Targeted Repair and Continuation Report');r.bold=True;r.font.size=Pt(22)
x=doc.add_paragraph();x.alignment=WD_ALIGN_PARAGRAPH.CENTER;x.add_run('定向排障 修复与已授权未完成任务续作').font.size=Pt(14)
x=doc.add_paragraph();x.alignment=WD_ALIGN_PARAGRAPH.CENTER;x.paragraph_format.space_before=Pt(25);x.add_run('版本日期 2026年9月12日\n基线 v3.1 strict background addendum 与 v4 历史执行包\n性质 已知结果后的纠错与探索性续作').font.size=Pt(10)
p(doc,'本报告不是新一轮扩展选题，也不是整篇稿件重写。历史 v4、v3.1 和更早文件保持只读；所有新输出均位于独立的 targeted_repair_2026-09-12 目录。',italic=True)
doc.add_page_break()

head(doc,'执行结论',1)
p(doc,'本轮完成了两个 P0 实现修复、一个计数显示修复、R 项目环境恢复、M04/M05 全基因差异分析，以及 M02 公开甲基化矩阵和 M03 公开空间转录组容器的实际获取与可行续作。结果不能概括为“全部阴性”：M04 固定程序与全基因比较存在明显来源差异，M05 的 AH 全基因对比存在信号，M03 的 Nageotte 结节与邻近神经元 ROI 存在供者一致的程序定位差异；但这些结果分别受到解剖和手术来源混杂、小样本组织特异性、以及缺乏明确 DPN 诊断与供者独立性证明的限制。')
p(doc,'M02 简化 beta 分支在两个 DPN 队列均未得到探针级 BH 支持；DMRcate 分别得到 167 和 660 个区域，但跨队列只有 3 对坐标重叠且方向一致为 0。PROPENG 中两个固定程序的区域集合富集未在预定发现队列 PROPGER 出现，因此不能称为复制。原需年龄和已验证技术批次的主要模型仍为受阻。')
doc.add_picture(str(R/'08_integration/TARGETED_REPAIR_SUMMARY.png'),width=Inches(6.7));cap(doc,'图1  定向修复与续作的核心数量结果。红色仅编码相应分析家族内 q<0.05。')

head(doc,'任务状态与边界',1)
st=read('STATUS.tsv');table(doc,['任务','状态','说明'],st[['task','status','detail']].itertuples(index=False,name=None))
p(doc,'REVIEW_CN.md、AUDIT_METRICS.json、evidence 目录和原始 SOURCE_VERIFICATION_NOTES.md 的独立字节未在附件或项目中找到。用户消息中的详细审核内容被作为执行规范使用，但没有被冒充为已恢复的原始文件。唯一可逐字节恢复的附件为 7,566 字节 EXECUTION_CONFIG_TEMPLATE.yaml，收件时间和哈希已单独登记。')

head(doc,'P0 修复结果',1)
head(doc,'M05 全表达排名背景',2)
b=read('02_M05_background_repair/results/M05_BACKGROUND_AUDIT.tsv')
table(doc,['审计指标','结果'],[
 ('原始特征行','56,609'),('实际 RNA 文库','15'),('可规范映射原始行','36,685'),
 ('全零原始行','13,420'),('聚合前纳入行','29,408'),('合格唯一 GeneID','29,373'),
 ('其中程序成员','3,993'),('其中非程序基因','25,380'),('历史程序限定背景','3,728'),
 ('重复特征聚合','不看结局的中位数'),('排名定义','全合格表达背景内平均百分位秩')])
p(doc,'历史实现先用程序成员表过滤 56,609 个原始特征，再进行百分位排名，实际背景只有 3,728 个程序成员。修复分支独立规范映射所有原始行，排除歧义、非有限和全零特征，以 GeneID 聚合后得到 29,373 个合格背景基因，其中 25,380 个为非程序基因。实际 RNA 文库为 15，length 未作为样本列。')
m5=read('02_M05_background_repair/results/M05_CORRECTED_EFFECTS_BY_MUSCLE.tsv');m5=m5[m5.component.eq('S')][['comparison','module_id','n_positive','n_negative','effect_positive_minus_negative','p','q_family']]
table(doc,['比较','程序','n正','n负','效应','P','家族q'],m5.itertuples(index=False,name=None))
p(doc,'修复后两个肌肉层面的 P5 完整评分比较均为名义 P=0.0286，但同一完整评分家族 BH q=0.257。其余完整评分同样无 q<0.05；供者配对交互也未获得校正支持。该结论保留所有正负方向，不因新背景而挑选较小 P 值。')
doc.add_picture(str(R/'02_M05_background_repair/figures/M05_corrected_program_effects.png'),width=Inches(6.7));cap(doc,'图2  M05 修复背景后的完整程序与两臂效应。')

head(doc,'固定模块映射与块删除',2)
mm=read('03_module_block_repair/results/MODULE_MAPPING_SUMMARY.tsv');table(doc,['模块','排除原因','GeneID数'],mm.itertuples(index=False,name=None))
p(doc,'原实现的第一列和最后一列静默回退已取消。MODULE_MEMBERS.tsv 被显式校验为 gene 与 module 两列，符号经统一 NCBI 规则映射；grey 作为未分配集合排除。九个固定非灰模块共映射 4,368 个 GeneID。M01、M04、M05 共 632 个块删除结果均可评价，原 72 行全 unassigned 且两臂为空的输出保留为历史错误路径，不再作为有效敏感性结果。')

head(doc,'M04 映射背景核查',2)
a=read('01_inputs/M04_MAPPING_AND_BACKGROUND_AUDIT.tsv')
table(doc,['审计指标','结果'],[
 ('原始特征行','58,174'),('规范映射后唯一 GeneID','35,749'),('全零 GeneID','15,430'),
 ('修复后纳入 GeneID','27,506'),('合格排名背景 GeneID','27,506'),('映射缓存记录','38,551'),
 ('不在 GSE302658 中的映射记录','1,198'),('缓存性质结论','一般注释宇宙，不是血液筛选子集'),
 ('实际修复触发','历史路径保留全零行')])
m4=read('01_inputs/M04_CORRECTED_EFFECTS.tsv');m4=m4[m4.component.eq('S')][['module_id','effect_positive_minus_negative','ci_low','ci_high','p','q_family']]
table(doc,['程序','效应','CI低','CI高','P','家族q'],m4.itertuples(index=False,name=None))
p(doc,'GSE302658_ensembl_to_ncbi 映射缓存包含 1,198 条不在 GSE302658 中的映射记录，证实其是一般注释表而非血液筛选子集。需要修复的是历史路径保留全零行，而不是缓存文件名本身。修复后 P6 与 P5 的家族 q 分别为 0.0009 和 0.0288；跨模块 133 项敏感性中仅 P6 保持 q<0.05。病例为足底趾神经，控制来自其他上下肢神经，诊断、解剖与手术来源完全混杂。')
doc.add_picture(str(R/'02_M05_background_repair/figures/M04_corrected_program_effects.png'),width=Inches(6.7));cap(doc,'图3  M04 一般背景审计后的固定程序与两臂效应。')

head(doc,'P1 计数与语义验证',1)
c=read('04_M03_count_repair/results/M03_COUNT_RECONCILIATION.tsv');table(doc,list(c.columns),c.itertuples(index=False,name=None))
p(doc,'供者 sub UTD 181030DHB 有两个 DRG，每个 DRG 各 2 张切片。修正表按 Sample ID 逐 DRG 计数，总数为 16 张切片、7 个 DRG、6 名供者；历史展示把供者总数 4 重复到两个 DRG 行，导致表内相加为 20。')

head(doc,'全基因差异分析续作',1)
d4=read('06_whole_gene_DE/results/M04_DE_SUMMARY.tsv');d5=read('06_whole_gene_DE/results/M05_DE_SUMMARY.tsv')
table(doc,['模块','模型或对比','检验基因','q<0.05','q<0.10','供者相关'],[(r.module,r.model,r.tested_genes,r.q_lt_005,r.q_lt_010,'') for r in d4.itertuples(index=False)]+[(r.module,r.contrast,r.tested_genes,r.q_lt_005,r.q_lt_010,r.duplicate_correlation) for r in d5.itertuples(index=False)])
p(doc,'M04 主要组别模型检验 16,702 个基因，其中 3,959 个 q<0.05；性别敏感性模型检验 17,321 个基因，其中 2,456 个 q<0.05。如此广泛的差异与来源完全混杂相容，不能称为 Morton 或共同神经损伤的疾病特异转录组。')
p(doc,'M05 供者阻断模型的 duplicate correlation 为 0.3918。AH 对比有 1,288 个基因 q<0.05，MG 与诊断乘肌肉交互均为 0。AH 单组织差异未被 MG 或交互支持，且诊断组供者数很小，因此只能作为组织特异的探索性来源差异。所有有限基因结果均保存在完整表中，没有只导出显著行。')
doc.add_picture(str(R/'06_whole_gene_DE/figures/M04_M05_DE_diagnostic.png'),width=Inches(6.7));cap(doc,'图4  M04 来源比较与 M05 AH 对比的全基因效应和原始 P 值。')

head(doc,'M02 甲基化输入恢复与探索分析',1)
p(doc,'GSE286347 beta 矩阵通过 GEO 官方下载入口完整取得。文件为 2,699,597,157 字节，SHA256 为 6c3b6a7f4481463d5cb9c7f320734582c2c284917e10e5a54ec8e7eb7dd5f0f4，gzip 全流校验通过。矩阵包含 315 个样本的 Beta 与 Detection P 成对列；231 个 DPN 样本均通过每样本 detection P 失败比例不超过 1% 的阈值。EPIC v1 hg19 注释版本 0.6.0 共 865,859 行。')
e=read('07_resource_recovery/M02_results/M02_EWAS_SUMMARY.tsv');table(doc,list(e.columns),e.itertuples(index=False,name=None))
ds=read('07_resource_recovery/M02_results/M02_DMR_STATUS.tsv');table(doc,['队列','状态','说明','DMR数'],ds.itertuples(index=False,name=None))
x=json.loads((R/'07_resource_recovery/M02_results/M02_DMR_CROSS_COHORT_SUMMARY.json').read_text());table(doc,['指标','值'],[(k,v) for k,v in x.items()])
p(doc,'两个队列均无 EWAS BH q<0.10 探针。DMRcate 在 PROPGER 和 PROPENG 分别产生 167 与 660 个区域，但只有 3 对区域坐标重叠，且方向一致为 0。固定 CpG 集合分析因预定 q<0.05 探针种子为零而不可评价，阈值未放宽。区域集合分析在 PROPENG 的 P2 与 P8 获得校正支持，在发现队列 PROPGER 均未支持，因此是 held out only 信号，不是方向复制。')
p(doc,'最重要的限制没有被输入恢复消除：GEO 不提供实际年龄和可复用、非混杂的技术批次，故原预设主要模型仍为 BLOCKED_RESOURCE。本轮 sex adjusted beta branch 在取得数据并完成完整性检查后、查看关联结果前单独锁定，只能作为观察性探索。')

head(doc,'M03 病理空间定位续作',1)
p(doc,'GSE295206_RAW.tar 通过官方入口取得，精确大小 9,698,846,720 字节，SHA256 为 1216d4d0fdfa0b26d8522ea8292793a9bc5fe5c5a97f50272773bbeb52c108c5。归档包含 80 个文件；按历史候选清单提取的 48 个 filtered H5、坐标和作者 ROI 文件全部与 tar 成员大小一致并逐个哈希。')
m3=read('07_resource_recovery/M03_results/M03_PROGRAM_TESTS.tsv');table(doc,['程序','供者n','效应','CI低','CI高','精确P','BH q'],m3[['program','n_donors','effect','CI_low','CI_high','P','q']].itertuples(index=False,name=None))
p(doc,'分析单位为供者：spot 先汇总到切片，再到 DRG，最后对供者等权；双侧符号翻转在 6 名供者上有 64 种精确排列。P2、P4、P5、P6、P7、P8 和 P10 的 BH q 为 0.0402。P5 和 P7 为负方向，其余为正方向。所有方向均保留，不能把它们合并为单一“神经元转录本减少”机制。')
p(doc,'作者 ROI 标签中的 Nageotte also touching neuron 被预先排除为混合区域；single 与 multiple neuron near Nageotte 合并为邻近神经元 ROI。首次实现因文字匹配未纳入复数 neurons 标签，质量复核后按锁定意图修复并全量重跑；首次输出保存在 superseded 目录。由于文件没有可验证的 spot diameter 校准，预定距离梯度未用近似单位替代，保持 NOT_EVALUABLE。')
p(doc,'该结果支持病理空间背景，而不是独立 DPN 疾病验证：公开资料仅说明糖尿病史，没有明确 DPN 诊断；六名供者是否与程序发现源重叠无法从公开键确定。因此不能与既有疾病关联拼接成因果机制链。')
doc.add_picture(str(R/'07_resource_recovery/M03_program_effects.png'),width=Inches(6.5));cap(doc,'图5  六名供者等权的 Nageotte 结节减邻近神经元 ROI 程序评分。')

head(doc,'验证 追溯与剩余问题',1)
t=read('tests/ALL_TARGETED_REPAIR_TESTS.tsv')
table(doc,['检查汇总','数量'],[('定向检查总数',34),('通过',33),('未通过',1)])
table(doc,['关键语义检查','观察值','状态'],[
 ('M05 RNA 文库计数','15，不含 length','通过'),('M05 非程序背景基因','25,380','通过'),
 ('固定模块块删除','632 条可评价','通过'),('M03 推断层级','16切片／7 DRG／6供者','通过'),
 ('M02 下载大小与哈希','完整匹配','通过'),('M03 下载大小与哈希','完整匹配','通过'),
 ('M02 预设主要模型','资源不足，保持受阻','通过'),('M03 距离梯度','无校准，保持不可评价','通过'),
 ('历史原始清单','161/162','未通过'),('排除打包后写日志的历史科学条目','161/161','通过')])
p(doc,'34 项定向检查中 33 项通过。唯一失败项是历史 v4 自身清单中的 logs/13_build_review_bundle.log：清单记录空文件哈希，但当前历史目录中的打包日志非空。排除这一打包后写日志后，161/161 个历史科学和管理条目均与清单一致。该差异被保留而非改写成全绿。')
p(doc,'M07 仍因 FinnGen 访问与端点定义缺失而 BLOCKED_ACCESS；M08 正式 ALC 与 EDC 患者键仍为 BLOCKED_LINKAGE，申请草稿没有发送。GitHub、Zenodo、DOI 与正式投稿未在本轮更新。')

head(doc,'对稿件整合的直接含义',1)
p(doc,'可纳入后续修正版稿件的新增证据包括：M01 两臂与真实固定模块删除的解释性敏感性；M04 修复后的来源比较，但必须在标题或正文中明确非 DPN、解剖与手术混杂；M05 修复背景后的无校正评分支持和 AH only 全基因差异；M02 双队列探针级无支持、跨队列 DMR 方向不一致以及 held out only 区域集合富集；M03 供者级病理 ROI 定位及其诊断和独立性边界；M06 原试验的有限无支持结果保持不变。')
p(doc,'不应恢复原 severity core 轴突丢失阳性，不应宣称程序分解优越、临床预测、因果机制、共同神经损伤已经独立证实或获得独立 DPN 空间验证。新增结果提高了分析完整性和方法审计可信度，但没有弥补独立临床验证不足。')

head(doc,'关键文件索引',1)
files=[
('科学锁定与修订','00_admin/'),('状态与科学摘要','STATUS.tsv；08_integration/TARGETED_REPAIR_SCIENTIFIC_SUMMARY.tsv'),
('Panel 与验证登记','08_integration/PANEL_SOURCE_REGISTER.tsv；tests/'),
('M05 旧新对照','02_M05_background_repair/results/M05_OLD_NEW_EFFECT_COMPARISON.tsv'),
('固定模块块删除','03_module_block_repair/results/'),('全基因 DE','06_whole_gene_DE/results/'),
('M02／M03 与日志','07_resource_recovery/M02_results；M03_results；logs/')]
table(doc,['内容','相对路径'],files)

# Word requires a terminal paragraph after a table. Compress the spacer so it
# remains on the preceding page instead of producing a trailing blank page.
if doc.paragraphs and not doc.paragraphs[-1].text:
 tail=doc.paragraphs[-1];tail.paragraph_format.space_before=Pt(0);tail.paragraph_format.space_after=Pt(0);tail.paragraph_format.line_spacing=Pt(1)
 rr=tail.add_run(' ');rr.font.size=Pt(1)

doc.core_properties.title='DPN v4 Targeted Repair and Continuation Report';doc.core_properties.subject='Targeted repair and continuation audit';doc.core_properties.author='DPN analysis team'
doc.save(OUT);print(OUT)
