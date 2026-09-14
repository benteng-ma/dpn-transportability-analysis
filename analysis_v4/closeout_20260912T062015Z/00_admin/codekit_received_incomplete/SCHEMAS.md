# 本地输入适配与代码接口

这些脚本不猜作者项目里的列名。Codex需读取实际源表和已执行代码，写一个小型adapter将真实字段转成以下接口。adapter必须保留原路径、hash、行键、转换原因；不得用报告里的四舍五入数字代替真实值。数据适配不是重新定义科学对象。

## 1. 安全调用方式
在代码包目录先执行：
```text
python -B code/seal_bundle.py verify .
python -B -m pytest -q -p no:cacheprovider tests
```
已有DPN环境满足依赖时不要重新安装。各命令输出目录必须不存在；历史目录不得作为输出。缺列/重复键/非法布尔会直接失败，不会用空值假装完成。TSV.gz亦可作为Python输入。

## 2. M04两级过滤审计
### RAW_FEATURE_DISPOSITION.tsv
一行一原始表达feature，必需列：
- feature_id：唯一原始行键，不是允许重复的gene symbol。
- gene_id：仅唯一映射行填写；歧义或无映射为空。
- mapping_status：unique / ambiguous / unmapped。
- finite_all：真实RNA列是否全部有限，true/false。
- all_zero：此原始行是否在全部真实RNA列为零，true/false。
- eligible_preaggregate：符合已声明唯一映射、finite、非全零规则，true/false。
- exclusion_reason：每个排除行唯一主原因；纳入行留空。

### CANONICAL_GENE_DISPOSITION.tsv
一行一个所有唯一映射得到的GeneID，**包括最后未纳入基因**：
- gene_id
- finite_after_aggregate / all_zero_after_aggregate / eligible_rank：true/false。
- n_source_rows：所有唯一映射到该基因的原始行数。
- n_eligible_source_rows：进入聚合的合格原始行数。

本函数验证适合此次报告的严格行级过滤→聚合规则。若已记录的项目算法在零值/有限值上顺序不同，应先查明事实、记录差异再适配；不能只为通过测试制造两个去向表。表达矩阵聚合本身要另用真实数据独立检查。

```text
python -B code/dpn_closeout.py flow --features RAW_FEATURE_DISPOSITION.tsv --genes CANONICAL_GENE_DISPOSITION.tsv --out new_M04_flow_audit
```

## 3. 固定模块删除汇总
标准化BLOCK_INFLUENCE.tsv的一行是一项真实已规划删除。必需列：
`analysis_id source endpoint program block original_effect deleted_effect eligible_before eligible_after rank_universe_same remaining_up remaining_down reason`

前三个来源/终点键与program/block共同唯一。program为P1..P10；block只能为固定九个非灰模块。NE的deleted_effect空或NA，reason说明。remaining_up/down是实际剩余数，不能填原成员数。rank_universe_same必须由实际缓存/代码核对，不能批量硬填true。

```text
python -B code/dpn_closeout.py blocks --input BLOCK_INFLUENCE.tsv --zero-tolerance 1e-10 --out new_block_summary
```
1e-10只作数值接近零和比值显示处理，不是生物效应阈值，不决定P值或样本删除。

## 4. M03供者级独立数值核查
标准化SPATIAL_DONOR_EFFECTS.tsv必需列：
`donor_id program effect up_effect down_effect eligible reason`

一行一供者×程序；必须从真实section配对→DRG→donor等权结果得到。`effect = up_effect - down_effect`。本helper不负责从原H5重新分区/过滤，它检查已适配供者输出。完整H5/ROI链条由Codex按任务C05重建。

```text
python -B code/dpn_closeout.py spatial --input SPATIAL_DONOR_EFFECTS.tsv --statistic mean --bootstrap-draws 2000 --seed 20260912 --minimum-donors 3 --out new_spatial_oracle
```
上面的mean仅演示参数位置。必须先确认原实际统计量；若原为studentized，传studentized。bootstrap是独立核查实现，随机流不同于旧代码可能导致CI轻微不同，不能默默替换旧CI。helper返回全部十程序；P9和缺失输入保持NE；当前有限家族按实际合格数BH。

## 5. 归档资产检查
```text
python -B code/dpn_closeout.py tar-inventory --tar GSE295206_RAW.tar --out new_tar_inventory
```
只列目录、大小、尺度/图像/ROI候选类型，不解压，不生成标定。未知位置的候选输入盘点：
```text
python -B code/dpn_closeout.py inventory --root "ACTUAL_DPN_ROOT" --out new_candidate_inventory
```
输出候选不代表自动选中版本。

## 6. M02安装实现审计
```text
Rscript code/audit_methylation_runtime.R new_R_runtime_audit
Rscript code/audit_methylation_runtime.R new_R_object_audit ACTUAL_CpGannotated.rds
```
第二种调用仅在原对象真实存在时用。脚本记录安装版本的formals/body和对象字段；不会虚构历史调用参数、不会执行DMR。

## 7. 固定区域复测 JOB.rds 的完整接口
`fixed_region_retest.R`的输入是list：

### beta
数字matrix；行名为通过固定QC的唯一CpG，列名为唯一PROPENG患者；只需提取全部德国固定区域相关CpG，不必把2.7GB全文再装入内存。必须从官方beta列取得，不得混入Detection P列。要求有限0..1，不自行插补。

### metadata
一行一实际患者，必须含`sample_id,cohort,pain_group`，和config声明的所有协变量。cohort全部PROPENG；pain_group仅painful/painless；严格匹配beta列名，禁止依顺序补键。sex值照真实编码，factor_levels显式列出；年龄缺失不伪造。

### regions
一行一个全部有效PROPGER发现区域，含：
`region_id,chr,start,end,discovery_direction,genome_build`
坐标1-based closed。discovery_direction为-1/0/1/NA，使用与英国显示同定义的德国beta效应；混合/无法验证不随意取最显著CpG方向。

### membership
一行一个固定region与探针的关系：
`region_id,probe_id,probe_chr,probe_position,genome_build`
保留全部德国原区域探针，不仅英国可测者。允许同一CpG落在不同固定区域，但region/probe必须唯一。位置须来自已核实EPIC注释，所有坐标单位一致。

### config
必须有以下全部字段，**示例不是已锁定运行配置**：
```r
list(
 locked_utc = "ACTUAL_UTC_TIMESTAMP_FROM_LOCAL_LOCK",
 role = "RETROSPECTIVE_FIXED_REGION_RETEST",
 discovery_cohort = "PROPGER", test_cohort = "PROPENG",
 genome_build = "hg19",
 coordinates_verified = TRUE, discovery_dmr_audit_passed = TRUE,
 qc_verified = TRUE, phenotype_columns_verified = TRUE,
 annotation_verified = TRUE,
 both_cohort_results_already_known = TRUE,
 covariates = c("sex"),
 covariate_types = c(sex="factor"),
 factor_levels = list(sex=c("ACTUAL_REFERENCE_LEVEL", "ACTUAL_OTHER_LEVEL")),
 branch = "SEX_ADJUSTED_EXPLORATORY",
 epsilon = 1e-6, min_cpgs = 3,
 min_probe_fraction = 0,  # 不额外新增覆盖百分比门槛；仍逐区报告覆盖
 min_residual_df = 1,
 ebayes_robust = TRUE, ebayes_trend = FALSE,
 input_provenance = list(beta="ACTUAL_FILE_AND_HASH", regions="ACTUAL_FILE_AND_HASH")
)
```
所有verified字段须有真实证据后填写，不能为了绕开gate直接设TRUE。若原锁定有明确region model/covariates/参数，优先重用其实际方法；上述参数仅为原来未具体规定区域统计量时的收尾实现建议，必须注明事后新增实现。M值均值的检验不是DMRcate原平滑检验。

```text
Rscript code/fixed_region_retest.R VERIFIED_PROPENG_FIXED_REGION_JOB.rds new_fixed_region_retest
```
脚本对整个德国区域名单作BH，包括NE造成的未检测区域计入family长度；P/CI为有效区域真实结果，NE仍NA。测试目标是患者层的平均区域M差异；调整后的beta效应是显示量，不做另一套可供择优的beta显著性筛选。

## 8. Word三线表辅助与封包
```text
python -B code/docx_three_line.py Input.docx Output_ThreeLine.docx
```
仅供一行表头的简单表块。复杂合并表头及跨页续表由Codex按实际布局处理。代码只检查文本保持和显式边框；必须再渲染检查，不得用“脚本运行”替代视觉验收。

所有日志关闭、图表与正文定稿后：
```text
python -B code/seal_bundle.py seal NEW_REVIEW_DIRECTORY --zip OUTSIDE_REVIEW_DIRECTORY.zip
```
ZIP不能在被打包目录内；manifest无自引用；会检查zip字节和原文件、源目录是否在打包中被改动，另写.sha256。封包日志只能在源目录外写。

## 9. 本包未做到的事
没有接入作者本机最新H5/beta/DMR对象；没有执行新的DPN科学分析；没有执行本包R脚本；没有使FinnGen或ALC自动获权。提供的Python单元测试是合成fixture，具体测试范围见verification，不是对真实供者结果的确认。
