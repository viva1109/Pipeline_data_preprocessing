

args<-commandArgs(trailingOnly = TRUE)
print(args)
wd<-args[1]

setwd(wd)
data_d<-args[2]
mths<-args[3]
print(data_d)
dir.create(data_d,recursive=T)
file_tree_ez<-args[4]
file_tree_silva<-args[5]
phenofile<-args[6]
tt_db<-args[7]
input_rank<-args[8]
print(tt_db)
# source("http://viva1109.iptime.org/RFunctions/Tree_based_methods.R")
# source(paste0("http://viva1109.iptime.org/RFunctions/FunctionsTMAT/Stat_Pipelines.R"),encoding = "UTF-8")
# source(paste0("http://viva1109.iptime.org/RFunctions/FunctionsTMAT/simulation_funcs.R"),encoding = "UTF-8")
source("/data/sharedcode/kjkim/RFunctions/Tree_based_methods.R")
source(paste0("/data/sharedcode/kjkim/RFunctions/FunctionsTMAT/Stat_Pipelines.R"),encoding = "UTF-8")
source(paste0("/data/sharedcode/kjkim/RFunctions/FunctionsTMAT/simulation_funcs.R"),encoding = "UTF-8")



get_set(1)
# ##########
# library(stringr)
# wd<-"D:/study/school/20160615 paper/Analysis/20171215NewSimul/realAnalysis/paperReal/RealTimeCal/"
# data_d<-"D:/study/school/GK/20180719Matching/data/aaaaa/se"
# mths<-"OMiAT_ANCOM_wilcoxon"
# methods<-str_split(mths,"_")[[1]]
# file_tree_ez<-"D:/study/school/20160615 paper/data/Data_1111/ez/ez_sina_fastmp.tre"
# file_tree_silva<-"D:/study/school/20160615 paper/data/Data_1111/silva128/97_otus.tre"
# dir.create(wd)
# setwd(wd)
# tt_db<-"ez"
# phenofile<-"pheno_togo_se.csv"
# ##########

# ##########
# library(stringr)
# wd<-"D:/study/school/20160615 paper/Analysis/20171215NewSimul/realAnalysis/paperReal/RealTimeCal/"
# data_d<-"D:/study/school/GK/20180719Matching/data/aaaaa/se"
# mths<-"OMiAT_ANCOM_wilcoxon"
# methods<-str_split(mths,"_")[[1]]
# file_tree_ez<-"D:/study/school/20160615 paper/data/Data_1111/ez/ez_sina_fastmp.tre"
# file_tree_silva<-"D:/study/school/20160615 paper/data/Data_1111/silva128/97_otus.tre"
# dir.create(wd)
# setwd(wd)
# tt_db<-"ez"
# phenofile<-"pheno_togo_se.csv"
# ##########


#########server
#wd<-"~/analysis/paper_20160615/Analysis/20190306_Real/PRJEB13092/"
#dir.create(wd,recursive=T)
#server_home<-"/home2/kjkim/analysis/paper_20160615/"
#data_d<-paste0(server_home,"data/Data_1111/PRJEB13092/")
#mths<-"TMAT15_TMAT16_ANCOM_edgeR_oMirkat_aMiSPU_OMiAT_wilcoxon"
#methods<-str_split(mths,"_")[[1]]
#file_tree_ez<-paste0(server_home,"./data/Data_1111/ez/ez_sina_fastmp.tre")
#file_tree_silva<-paste0(server_home,"./data/Data_1111/silva128/97_otus.tre")
#setwd(wd)
#phenofile<-"pheno_togo.csv"
#tt_db<-"silva"
#########


##########
# library(stringr)
# wd<-"D:/study/school/20160615 paper/Analysis/20181222TMATFunc/"
# data_d<-"D:/study/school/20160615 paper/data/Data_1111/PRJEB6070"
# mths<-"TMAT15"
# methods<-str_split(mths,"_")[[1]]
# file_tree_ez<-"D:/study/school/20160615 paper/data/Data_1111/ez/ez_sina_fastmp.tre"
# file_tree_silva<-"D:/study/school/20160615 paper/data/Data_1111/silva128/97_otus.tre"
# dir.create(wd)
# setwd(wd)
# tt_db<-"ez"
# phenofile<-"pheno_togo.csv"
##########

