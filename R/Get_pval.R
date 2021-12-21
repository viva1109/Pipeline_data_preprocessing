args<-commandArgs(trailingOnly = TRUE)
print(args)
wd<-args[1]
setwd(wd)
db<-args[2]

nperm<-as.numeric(args[3])
methods_bf<-args[4]
ncol<-as.numeric(args[5])
ncore<-as.numeric(args[6])
conti<-try(as.logical(as.numeric(args[7])))
ancom_FDR<-as.numeric(args[8])

# wd<-"D:/study/school/20160615 paper/Analysis/20171215NewSimul/realAnalysis/"
# # wd<-"~/analysis/paper_20160615/Analysis/20171215NewSimul/realAnalysis/paperReal/PRJEB6070NEW"

# methods_bf<-"OMiAT_wilcoxon_TMAT15"

# methods_bf<-c("TMAT15_TMAT16_ANCOM_edgeR_oMirkat_aMiSPU_OMiAT_wilcoxon")
# db<-"silva"
# nperm<-10
# ncol<-3
# ncore<-1
# conti<-F


source(paste0("/data/sharedcode/kjkim/RFunctions/FunctionsTMAT/Stat_Pipelines.R"),encoding = "UTF-8")
get_pval(1)

# library(stringr)
# wd<-"~/analysis/paper_20160615/Analysis/20171215NewSimul/realAnalysis/paperReal/PRJEB6070NEW_0204/"
# dir.create(wd,recursive=T)
# server_home<-"/home2/kjkim/analysis/paper_20160615/"
# data_d<-paste0(server_home,"data/Data_1111/PRJEB6070/")
# methods_bf<-"TMAT15"
# methods<-str_split(methods_bf,"_")[[1]]
# file_tree_ez<-paste0(server_home,"./data/Data_1111/ez/ez_sina_fastmp.tre")
# file_tree_silva<-paste0(server_home,"./data/Data_1111/silva128/97_otus.tre")
# setwd(wd)
# phenofile<-"pheno_togo.csv"
# tt_db<-"silva"
# ncol<-3
# ncore<-24
# conti<-F
# db<-"silva"
# nperm<-10

###Windows
# wd<-"D:/study/school/20160615 paper/Analysis/20181222TMATFunc/"
# dir.create(wd,recursive=T)
# server_home<-"D:/study/school/20160615 paper/"
# data_d<-paste0(server_home,"data/Data_1111/PRJEB6070/")
# methods_bf<-"OMiAT_wilcoxon_TMAT15"
# methods_bf<-"TMAT15"
# methods<-str_split(mths,"_")[[1]]
# file_tree_ez<-paste0(server_home,"./data/Data_1111/ez/ez_sina_fastmp.tre")
# file_tree_silva<-paste0(server_home,"./data/Data_1111/silva128/97_otus.tre")
# setwd(wd)
# phenofile<-"pheno_togo.csv"
# tt_db<-"ez"
# ncol<-3
# ncore<-1
# conti<-F
# db<-"ez"
# nperm<-10
