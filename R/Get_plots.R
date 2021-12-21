
args<-commandArgs(trailingOnly = TRUE)

folder<-args[1]

nm_y_bf<-args[2]

ncol<-args[3]

db_bf<-args[4]
methods_bf<-args[5]

# methods_bf<-"wilcoxon_TMAT15_OMiAT_ANCOM"
# methods_bf<-"TMAT15_OMiAT_wilcoxon"
# methods_bf<-"ANCOM_TMAT15"
# folder<- "~/analysis/paper_20160615/Analysis/20171215NewSimul/realAnalysis/paperReal/PRJEB6070NEW/"
# folder<- "~/analysis/paper_20160615/Analysis/20190306_Real/PRJEB6070/"
# nm_y_bf<-"Control_Case"
# ncol<-"3"
# db_bf<- "silva"

source(paste0("/data/sharedcode/kjkim/RFunctions/FunctionsTMAT/Stat_Pipelines.R"),encoding = "UTF-8")
get_plots(1)



# methods_bf<-"wilcoxon_OMiAT_TMAT15"
# folder<- "~/analysis/paper_20160615/Analysis/20171215NewSimul/realAnalysis/paperReal/PRJEB6070NEW/"
# # folder<- "/home2/kjkim/analysis/20180115CHOI/data/output/SP_Cont3"
# nm_y_bf<-"Control_Case"
# ncol<-"3"
# db_bf<- "ez"
     