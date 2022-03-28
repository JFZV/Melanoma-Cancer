library(TCGAbiolinks)
library(SummarizedExperiment)
# 下载的基因表达、DNA甲基化、基因突变数据都是关于hg38
# 下载基因表达数据
query_exp <- GDCquery(project = "TCGA-SKCM", 
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "HTSeq - FPKM") # 下载FPKM数据
GDCdownload(query_exp) 
data_exp <- GDCprepare(query_exp) 
exp_assay <- as.data.frame(assay(data_exp))

saveRDS(exp_assay, "./data/TCGA-SKCM/EXP/HTSeq_FPKM.RData") # 基因表达矩阵
exp_colData <- as.data.frame(colData(data_exp))
saveRDS(exp_colData, "./data/TCGA-SKCM/EXP/colData.RData") # 病人注释（472个病人）
exp_rowRanges <- as.data.frame(rowRanges(data_exp))
saveRDS(exp_rowRanges, "./data/TCGA-SKCM/EXP/rowRanges.RData") # 基因注释

# 下载甲基化数据
query_mty <- GDCquery(project = "TCGA-SKCM",
                      data.category = "DNA Methylation",
                      platform = "Illumina Human Methylation 450")
GDCdownload(query_mty)
data_mty <- GDCprepare(query_mty)
mty_assay <- as.data.frame(assay(data_mty))
mty_colData <- as.data.frame(colData(data_mty))
saveRDS(mty_colData, "./data/TCGA-SKCM/MTY/colData.RData") # 病人注释（475个病人）
mty_rowRanges <- as.data.frame(rowRanges(data_mty))
saveRDS(mty_rowRanges, "./data/TCGA-SKCM/MTY/rowRanges.RData") # cg探针注释（探针对应着多个基因）

# 下载基因突变数据
maf<-GDCquery_Maf("SKCM", pipelines = "muse") # 下载的是mutation annotation文件
saveRDS(maf, "./data/TCGA-SKCM/SNV/maf.RData")


# 下载临床数据
clinical <- GDCquery_clinic(project = "TCGA-SKCM", type = "clinical")
saveRDS(clinical, "./data/TCGA-SKCM/clinical.RData")

query <- GDCquery(project = "TCGA-SKCM", 
                  data.category = "Clinical",
                  data.type = "Clinical Supplement",
                  data.format = "BCR Biotab"
)
GDCdownload(query)
sckm.tab.all <- GDCprepare(query)

therapy <- sckm.tab.all$clinical_drug_skcm
saveRDS(therapy, "./data/TCGA-SKCM/therapy.RData")
therapy$pharmaceutical_therapy_type # therapy types

radiation <- sckm.tab.all$clinical_radiation_skcm
saveRDS(radiation, "./data/TCGA-SKCM/radiation.RData")