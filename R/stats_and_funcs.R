wilcoxon_output<-function(y,X_P_rel3){
	wilcoxon<-sapply(data.frame(X_P_rel3),function(data,ty){
	rree<-suppressWarnings( wilcox.test(data ~ ty) )
	return(rree$p.value)
	 },y)
	wilcoxon<-wilcoxon[!is.nan(wilcoxon)]
	wilcoxon_minp<-exact_pval(wilcoxon)
	wilcoxon_fisher<-fisher_pval(wilcoxon)
	return(c(wilcoxon_minp,wilcoxon_fisher))
}

fisher_pval<-function(pval_1){
  len_tt<-length(pval_1)
  tt_1<- -2*sum(log(pval_1))
  1-pchisq(tt_1,df=2*len_tt)
}

getTtest<-function(t_y,t_Z){
  tre<-t.test(t_Z~t_y)
  tre$p.value
}

logcpm<-function(data){
  nomi<-data+0.5
  denomi<-apply(data,1,sum)+1
  output<-log2(nomi/denomi*10^10)
  return(output)
}

get_naive_t<-function(data,y,A,yAy){
  # browser()
  new_Z<-matrix(data,ncol=1)
  new_t_re<-getT(y,A,new_Z,yAy)
  1-pchisq(new_t_re,df=1)
}



real_0921<-function(t_pheno,subtree,datainput,nperm=2000,type,candi_Anal,index_whatTo,include_Toprank=F,total.reads,cov=NULL,t_table=NULL,simul=NULL){
 #browser()
  whatToAnal<-candi_Anal[index_whatTo]
  pruned.tree_perphy<-subtree[[1]]
  lili<-call_whattodo_v2(pruned.tree_perphy)
  ind_gogo<-match(subtree[[2]]$code,colnames(datainput))
  X_P_bf<-datainput
  X_P_bf3<-datainput[,ind_gogo]
  n_indv<-dim(X_P_bf3)[1]
  p<-dim(cov)[1]
 
  if(!is.null(simul)){
    beta_per<-t_table$beta_per_list
    y <- t_pheno[simul$pheno_index]
    if(t_table$plus_minus) beta_per<-sample(c(1,-1),1)*beta_per
    X_P_bf3<-generateH1(y,X_P_bf3,simul$ind_causal,beta_per)

    #X_P_bf[,ind_gogo]<-X_P_bf3
  }else{
	  y <- t_pheno
  }

  if("wilcoxon_minp"%in%whatToAnal){
    X_P_rel<-X_P_bf/apply(X_P_bf,1,sum)
    X_P_rel3<-X_P_rel[,ind_gogo]
  }

  if("oMirkat"%in%whatToAnal){
    otu.tab.rff <- Rarefy(X_P_bf)$otu.tab.rff
    otu.tab.rff_genus <- otu.tab.rff[,ind_gogo]
  }
  if("aMiSPU"%in%whatToAnal){
    X = as.matrix(X_P_bf3)
    ind_rowsumnot0<-rowSums(X)!=0
    X_fixed<-X[ind_rowsumnot0,]
  }

  if("Min_p_chi"%in%whatToAnal){
	if(type%in%c(1,3)){
	  X_P_fbf<-X_P_bf
	}else if (type%in%c(2,4,5,9:11)){
	  X_P_fbf<-t(cpm(t(X_P_bf),log=TRUE))
	}else{
	  X_P_fbf<-logcpm(X_P_bf)
	}
	
	X_P<-X_P_fbf[,ind_gogo]
  }
  if(include_Toprank){
      X_P_comp<-apply(X_P_fbf[,-ind_gogo],1,sum)
  }
  out_omiat_p<-NA
  out_mirkat_p<-NA
  out_mispu_p<-NA

  if("wilcoxon_minp"%in%whatToAnal){
    toutput<-wilcoxon_output(y,X_P_rel3)
    wilcoxon_minp<-toutput[1]
    wilcoxon_fisher<-toutput[2]
  }

  if("oMirkat"%in%whatToAnal){
    out_mirkat <- out_MirKAT_result(otu.tab.rff_genus,pruned.tree_perphy,y,nperm=nperm)
    out_mirkat_p <- out_mirkat$omnibus_p
  }
  if("aMiSPU"%in%whatToAnal){
    phenos_fixed<-y[ind_rowsumnot0]
    out_mispu <-MiSPU(phenos_fixed,X_fixed, pruned.tree_perphy,model = "binomial", pow = c(2:8, Inf), n.perm = nperm,cov=cov)
    out_mispu_p<-out_mispu$aMiSPU$pvalue
  }
  if("OMiAT"%in%whatToAnal){
    ind_out<-apply(X_P_bf3,2,sum)==0
    if(sum(!ind_out)==1){
	out_omiat_p<-NA
    }else{
        pruned.tree_perphy_omiat<-drop.tip(pruned.tree_perphy,colnames(X_P_bf3)[ind_out])
	 out_omiat<-OMiAT(Y=y, otu.tab=X_P_bf3[,!ind_out], tree=pruned.tree_perphy_omiat,model="binomial", n.perm=nperm,total.reads=total.reads, cov = cov)
	 out_omiat_p<-out_omiat$OMiAT.pvalue
    }
  }
  if("Min_p_chi"%in%	whatToAnal){
    tp_pv<- get_PMT_out_new0918(lili,y,X_P,NULL,t_type=type,cov=cov,total.reads=total.reads)
    pval_1<-tp_pv[1,]
    pval_2<-tp_pv[2,]
    if(include_Toprank){
      C1<-apply(X_P,1,sum)
      C2<-X_P_comp
      dm<- type_to_dm(C1=C1,C2=C2,type=type,total.reads=total.reads)
      genus_t<-getT2(y,dm[,1],cov=cov)
      genus_pval<-1-pchisq(genus_t,df=1)
      pval_1<-c(pval_1,genus_pval)
      pval_2<-c(pval_2,1-pf(genus_pval,df1 = 1,df2 = n_indv-p))
    }
    pval_minp_ch<-exact_pval(pval_1)
    pval_minp_f<-fisher_pval(pval_2)
    pval_fisher_ch<-fisher_pval(pval_1)
    pval_fisher_f<-fisher_pval(pval_2)
  }else{
    pval_minp_ch<-NA
    pval_minp_f<-NA
    pval_fisher_ch<-NA
    pval_fisher_f<-NA
  }

  tmp_pval_1<-sapply(data.frame(X_P),function(data,ty){
      t_re_1<-getT2(y,data,cov=cov)
      tpval_1<-pchisq(t_re_1,df=1)
      return(tpval_1)
    },y)

  score_pval<-exact_pval(c(tmp_pval_1,genus_pval))
  X_P_rel<-X_P_bf/apply(X_P_bf,1,sum)
  X_P_rel3<-X_P_rel[,ind_gogo]
  data_wil<-apply(X_P_rel3,1,sum)
  genus_wilcoxon<-suppressWarnings( wilcox.test(data_wil ~ y)$p.value     )

if("treeBased_0" %in% whatToAnal){
  tmp_1<-get_PMT_treebased_0914(lili,t_pheno,X_P,NULL,noCrossTerm = F,tree=pruned.tree_perphy,useWeight = T)
  treeBased_0<-tmp_1[1]
  treeBased_1<-tmp_1[2]
  
  tmp_1<-get_PMT_treebased_0914(lili,t_pheno,X_P,NULL,noCrossTerm = T,tree=pruned.tree_perphy,useWeight = T)
  treeBased_0_NC<-tmp_1[1]
  treeBased_1_NC<-tmp_1[2] 
  
  tmp_1<-get_PMT_treebased_0914(lili,t_pheno,X_P,NULL,noCrossTerm = F,tree=pruned.tree_perphy,useWeight = F)
  treeBased_0_nW<-tmp_1[1]
  treeBased_1_nW<-tmp_1[2] 
  
  tmp_1<-get_PMT_treebased_0914(lili,t_pheno,X_P,NULL,noCrossTerm = T,tree=pruned.tree_perphy,useWeight = F)
  treeBased_0_NC_nW<-tmp_1[1]
  treeBased_1_NC_nW<-tmp_1[2] 
}else{
    treeBased_0<-NA
    treeBased_1<-NA
    treeBased_0_NC<-NA
    treeBased_1_NC<-NA
    treeBased_0_nW<-NA
    treeBased_1_nW<-NA
    treeBased_0_NC_nW<-NA
    treeBased_1_NC_nW<-NA 
}
  # result_tmp<-c( NA,exact_pval(pval_1),exact_pval(pval_2),exact_pval(pval_3),exact_pval(wilcoxon),min(p.adjust(pval_1,method = "fdr")),min(p.adjust(pval_2,method = "fdr")),min(p.adjust(pval_3,method = "fdr")),min(p.adjust(wilcoxon,method = "fdr")),min(p.adjust(pval_1,method = "bonferroni")),min(p.adjust(pval_2,method = "bonferroni")),min(p.adjust(pval_3,method = "bonferroni")),min(p.adjust(wilcoxon,method = "bonferroni")),genus_pval,out_mirkat$omnibus_p,out_mispu$aMiSPU$pvalue,out_omiat$OMiAT.pvalue)
  
  result_tmp<-c( NA,pval_minp_ch,pval_minp_f,pval_fisher_ch,pval_fisher_f,out_mirkat_p,out_mispu_p,out_omiat_p,score_pval,wilcoxon_minp,wilcoxon_fisher,genus_pval,genus_wilcoxon,treeBased_0,treeBased_1,treeBased_0_NC,treeBased_1_NC,treeBased_0_nW,treeBased_1_nW,treeBased_0_NC_nW,treeBased_1_NC_nW)[c(1,(index_whatTo+1))]
  return(result_tmp)
}

real_1129_findWhere<-function(t_pheno,subtree,datainput,nperm=2000,type,candi_Anal,index_whatTo,include_Toprank=F,total.reads,cov=NULL){
  #browser()
  whatToAnal<-candi_Anal[index_whatTo]
  pruned.tree_perphy<-subtree[[1]]
  lili<-call_whattodo_v2(pruned.tree_perphy)
  ind_gogo<-match(subtree[[2]]$code,colnames(datainput))
  X_P_bf<-datainput
  # X_P_bf2<-datainput[,ind_gogo]
  X_P_bf3<-datainput[,ind_gogo]
  n_indv<-dim(X_P_bf3)[1]
  v_one_n<-rep(1,n_indv)
  A <- diag(n_indv)-v_one_n%*%solve(t(v_one_n)%*%v_one_n)%*%t(v_one_n)
  y <- t_pheno
  yAy <- t(y)%*%A%*%y

  
  if("wilcoxon_minp"%in%whatToAnal){
    X_P_rel<-X_P_bf/apply(X_P_bf,1,sum)
    X_P_rel3<-X_P_rel[,ind_gogo]
    
    wilcoxon<-sapply(data.frame(X_P_rel3),function(data,ty){
      rree<-suppressWarnings( wilcox.test(data ~ ty) )
      return(rree$p.value)
    },y)
    wilcoxon_minp<-exact_pval(wilcoxon)
    wilcoxon_fisher<-fisher_pval(wilcoxon)
  }else{
    wilcoxon_minp<-NA
    wilcoxon_fisher<-NA
  }
  
  
  
  out_omiat_p<-NA
  out_mirkat_p<-NA
  out_mispu_p<-NA
  
  if("oMirkat"%in%whatToAnal){
    out_mirkat<-out_MirKAT_result(X_P_bf,pruned.tree_perphy,y,nperm=nperm)
    out_mirkat_p<-out_mirkat$omnibus_p
  }
  if("aMiSPU"%in%whatToAnal | "OMiAT"%in%whatToAnal){
    X = as.matrix(X_P_bf3)
    ind_rowsumnot0<-rowSums(X)!=0
    X_fixed<-X[ind_rowsumnot0,]
    phenos_fixed<-y[ind_rowsumnot0]
  }
  if("aMiSPU"%in%whatToAnal){
    out_mispu <-MiSPU(phenos_fixed,X_fixed, pruned.tree_perphy,model = "binomial", pow = c(2:8, Inf), n.perm = nperm)
    out_mispu_p<-out_mispu$aMiSPU$pvalue
  }
  if("OMiAT"%in%whatToAnal){
    ind_out<-apply(X_P_bf3,2,sum)==0
    if(sum(!ind_out)==1){
	out_omiat_p<-NA
    }else{
        pruned.tree_perphy_omiat<-drop.tip(pruned.tree_perphy,colnames(X_P_bf3)[ind_out])
	 out_omiat<-OMiAT(Y=y, otu.tab=X_P_bf3[,!ind_out], tree=pruned.tree_perphy_omiat,model="binomial", n.perm=nperm,total.reads=total.reads, cov = cov)
	 out_omiat_p<-out_omiat$OMiAT.pvalue
    }
  }
  
  #get the selected one
  if(type%in%c(1,3)){
    X_P_fbf<-X_P_bf
  }else if (type%in%c(2,4,5,9:11)){
    X_P_fbf<-t(cpm(t(X_P_bf),log=TRUE))
  }else{
    X_P_fbf<-logcpm(X_P_bf)
  }
  
  X_P<-X_P_fbf[,ind_gogo]
  
  if("Min_p_chi"%in%whatToAnal){
    tp_pv<- get_PMT_out_new0918(lili,t_pheno,X_P,NULL,t_type=type,cov=cov)
    pval_1<-tp_pv[1,]
    pval_2<-tp_pv[2,]
    pval_3<-tp_pv[3,]
    
    if(include_Toprank){
      X_P_comp<-apply(X_P_fbf[,-ind_gogo],1,sum)
      C1<-apply(X_P,1,sum)
      C2<-X_P_comp
      dm<-matrix(nrow=length(C1),ncol=2)
      if(type%in%c(1:2,7)){
        Cout<-cbind(C1,C2)
        dm<-t(cpm(t(Cout),log=TRUE))
      }else if (type%in%c(3,5:6)){
        Cout<-cbind(C1,C2)
        dm<-logcpm(Cout)
      }else if (type%in%c(4,8)){
        dm[,1]<-log(C1/(C1+C2))
        dm[,2]<-log(C2/(C1+C2))
      }else if (type%in%c(9)){
        dm[,1]<-(C1/(C1+C2))
        dm[,2]<-(C2/(C1+C2))
      }else if (type%in%c(10)){
        dm[,1]<-C1/C2
        dm[,2]<-C2/C1
      }else if (type%in%c(11)){
        dm[,1]<-log(C1/C2)
        dm[,2]<-log(C2/C1)
        # dm[,1]<-(C1/(C1+C2))
        # dm[,2]<-(C2/(C1+C2))
      }
      t_re_1<-getT2(y,dm[,1])
      pval_1<-c(pval_1,1-pchisq(t_re_1,df=1))
      pval_2<-c(pval_2,1-pf(t_re_1,df1 = 1,df2 = n_indv-1))
      meanre<-mean_group<-tapply(dm[,1],y,mean)
      pval_3<-c(pval_3,meanre[1]<meanre[2])
    }
    pval_minp_ch<-exact_pval(pval_1)
    pval_minp_f<-fisher_pval(pval_2)
    pval_fisher_ch<-fisher_pval(pval_1)
    pval_fisher_f<-fisher_pval(pval_2)
  }else{
    pval_minp_ch<-NA
    pval_minp_f<-NA
    pval_fisher_ch<-NA
    pval_fisher_f<-NA
  }
  
  
  
  # genus_pval<-tp_pv[1,1]
  			  if("treeBased_0" %in% whatToAnal){
  tmp_1<-get_PMT_treebased_0914(lili,t_pheno,X_P,NULL,noCrossTerm = F,tree=pruned.tree_perphy,useWeight = T)
  treeBased_0<-tmp_1[1]
  treeBased_1<-tmp_1[2]
  
  tmp_1<-get_PMT_treebased_0914(lili,t_pheno,X_P,NULL,noCrossTerm = T,tree=pruned.tree_perphy,useWeight = T)
  treeBased_0_NC<-tmp_1[1]
  treeBased_1_NC<-tmp_1[2] 
  
  tmp_1<-get_PMT_treebased_0914(lili,t_pheno,X_P,NULL,noCrossTerm = F,tree=pruned.tree_perphy,useWeight = F)
  treeBased_0_nW<-tmp_1[1]
  treeBased_1_nW<-tmp_1[2] 
  
  tmp_1<-get_PMT_treebased_0914(lili,t_pheno,X_P,NULL,noCrossTerm = T,tree=pruned.tree_perphy,useWeight = F)
  treeBased_0_NC_nW<-tmp_1[1]
  treeBased_1_NC_nW<-tmp_1[2] 
}else{
    treeBased_0<-NA
    treeBased_1<-NA
    treeBased_0_NC<-NA
    treeBased_1_NC<-NA
    treeBased_0_nW<-NA
    treeBased_1_nW<-NA
    treeBased_0_NC_nW<-NA
    treeBased_1_NC_nW<-NA 
}
  # result_tmp<-c( NA,exact_pval(pval_1),exact_pval(pval_2),exact_pval(pval_3),exact_pval(wilcoxon),min(p.adjust(pval_1,method = "fdr")),min(p.adjust(pval_2,method = "fdr")),min(p.adjust(pval_3,method = "fdr")),min(p.adjust(wilcoxon,method = "fdr")),min(p.adjust(pval_1,method = "bonferroni")),min(p.adjust(pval_2,method = "bonferroni")),min(p.adjust(pval_3,method = "bonferroni")),min(p.adjust(wilcoxon,method = "bonferroni")),genus_pval,out_mirkat$omnibus_p,out_mispu$aMiSPU$pvalue,out_omiat$OMiAT.pvalue)
  
  result_tmp<-c( NA,pval_minp_ch,pval_minp_f,pval_fisher_ch,pval_fisher_f,out_mirkat_p,out_mispu_p,out_omiat_p,wilcoxon_minp,wilcoxon_fisher,treeBased_0,treeBased_1,treeBased_0_NC,treeBased_1_NC,treeBased_0_nW,treeBased_1_nW,treeBased_0_NC_nW,treeBased_1_NC_nW)[c(1,(index_whatTo+1))]
  return(list(pval=result_tmp,min_index=pval_1,testnodes=lili,dir=pval_3))
}



getOutput<-function(ind_what,lout_list,t_list_upper,prefix){
  # browser()
  What<-t_list_upper[ind_what]
  lout<-lout_list[[ind_what]]
  l2<-lapply(1:length(lout),function(ind){
    ttt<-apply(lout[[ind]][,-1],2,as.numeric)
    dim(ttt)[2]
    df<-matrix(ncol=dim(ttt)[2],nrow=dim(ttt)[2])
    for(i in 1:(dim(ttt)[2]-1)){
      for ( k in (i+1):dim(ttt)[2]){
        df[i,k]<-cor(ttt[,i],ttt[,k])
      }
    }
    colnames(df)<-paste0("T",1:dim(ttt)[2])
    # write.csv(df,paste0("Genus_pval_corr",ind,".csv"),row.names = F)
    return(df)
  })
  
  lll1<-lapply(l2,function(data){
    dim(data)<-NULL
    return(data[!is.na(data)])
  })
  
  func_lll1_2<-function(ind){
    # browser()
    write(paste0(ind,"\n"),paste0(filedir,"/midoutput/ongoing/debug_func_lll1_2.txt"),append = T)
    data<-l2[[ind]]
    dim(data)<-NULL
    out<-data[!is.na(data)]
    oott<-cbind(paste0("Genus",ind),data.frame(out))
    names(oott)<-c("Name","correlation")
    return(oott)
  }
  
  lll1_2<-lapply(1:length(l2),func_lll1_2)
  
  lll1_out<-do.call("rbind",lll1_2)
  head(lll1_out)
  
  png(paste0(filedir,"/midoutput/",prefix,"_",What,"_box_byTaxon.png"),width = 500, height = 800)  
  boxplot(correlation~Name,data=lll1_out)
  abline(h=0,col="blue")
  dev.off()
  
  lll2<-do.call("c",lll1)
  png(paste0(filedir,"/midoutput/",prefix,"_",What,"_box_allTaxon.png"),width = 500, height = 800)  
  boxplot(lll2)
  abline(h=0,col="blue")
  dev.off()
  
  
  
  outsep_pvales<-lapply(lout,function(data){
    unlist(data[,-1])
  }
  )
  
  outsep_pvales2<-do.call("c",outsep_pvales)
  
  png(paste0(filedir,"/midoutput/",prefix,"_",What,"_rawpval.png"),width = 1000, height = 700)  
  hist(as.numeric(outsep_pvales2),breaks=50)
  dev.off()
  # hist(as.numeric(outsep_pvales2),breaks=100)
  
  out_combpval<-function(data){
    dd<-data
    dd1<-dd
    dd1[,-1]<-apply(dd[,-1],2,as.numeric)
    apply(dd1[,-1],1,exact_pval)
    outdf<-cbind(apply(dd1[,-1],1,exact_pval),apply(dd1[,-1],1,fisher_pval))
  }
  
  outcomb1<-lapply(lout,out_combpval)
  
  
  jjoutcomb1<-do.call("rbind",outcomb1)
  png(paste0(filedir,"/midoutput/",prefix,"_",What,"_Minpval.png"),width = 1000, height = 700)  
  hist(jjoutcomb1[,1],main="Minimum pval",breaks = 50)
  dev.off()
  
  
  png(paste0(filedir,"/midoutput/",prefix,"_",What,"_Fishval.png"),width = 1000, height = 700)  
  hist(jjoutcomb1[,2],main="Fisher's method",breaks = 50)
  dev.off()
  
  
  outcomb2<-lapply(outcomb1,function(data){
    apply(data,2,function(data2){
      c(sum(data2<0.1),sum(data2<0.05),sum(data2<0.01),sum(data2<0.005))/length(data2)
    })
  })
  
  outcomb3<-do.call("rbind",outcomb2)
  Genus<-c()
  alpha<-c()
  for ( i in 1:length(phylums_candidates)){
    Genus<-c(Genus,rep(i,4))
    alpha<-c(alpha,c(0.1,0.05,0.01,0.005))
  }
  
  outcomb4<-cbind(Genus,alpha,outcomb3)
  colnames(outcomb4)[3:4]<-c("exact","fisher")
  
  
  outcomb5_bf<-apply(jjoutcomb1,2,function(data2){
    c(sum(data2<0.1),sum(data2<0.05),sum(data2<0.01),sum(data2<0.005))/length(data2)
  })
  
  
  outcomb5<-cbind("AllGenus",c(0.1,0.05,0.01,0.005),outcomb5_bf)
  colnames(outcomb5)<-c("WhatGenus","alpha","exact","fisher")
  
  
  write.csv(outcomb4,paste0(filedir,"/midoutput/csv/",prefix,"_",What,"_outcomb4.csv"))
  
  write.csv(outcomb5,paste0(filedir,"/midoutput/csv/",prefix,"_",What,"outcomb5.csv"))
  
}



# ww<-tmp_ind$ww_list
# n_sim=n_sim
# beta_per=tmp_ind$beta_per_list
# propor_factor=tmp_ind$propor_factor_list
# permu=go_phylum_setting[[ind_data]]
# nperm=npermSet

gotest_real_2<-function(ww,n_sim,beta_per,propor_factor,permu,nperm,type,index_whatTo,candi_Anal){
  write(paste0(ww," "),paste0(filedir,"/midoutput/ongoing/ong_",beta_per,"_",propor_factor,"_",q_input,".txt"),append = T)
  return(t(sapply(permu$permu_indcausal_pheno,gogo_power_0925,permu$subtree_tmp,sorted_df_table_formed,beta_per,nperm,type=type,candi_Anal=candi_Anal,index_whatTo=index_whatTo)))
}



powerTest<-function(ww,n_sim,beta_per,propor_factor,permu,nperm,type,index_whatTo,candi_Anal,include_Toprank=F){
  write(paste0(ww," "),paste0(filedir,"/midoutput/ongoing/ong_",beta_per,"_",propor_factor,"_",q_input,".txt"),append = T)
  return(t(sapply(permu$permu_indcausal_pheno,get_statsitics_0925,permu$subtree_tmp,sorted_df_table_formed,beta_per,nperm,type=type,candi_Anal=candi_Anal,index_whatTo=index_whatTo,include_Toprank=include_Toprank,simul=T)))
}


powerTest_H0<-function(ww,n_sim,permu_list,nperm,type,index_whatTo,candi_Anal){
  write(paste0(ww," "),paste0(filedir,"/midoutput/ongoing/H0_ong_","_",q_input,".txt"),append = T)
  subtree_tmp<-gettheSubtree(ww,phylum_list_perphy,phylums_candidates,pruned.tree2)
  return(t(sapply(permu_list,get_statsitics_0925,subtree_tmp,sorted_df_table_formed,0,nperm,type=type,candi_Anal=candi_Anal,index_whatTo=index_whatTo,include_Toprank=F,simul=F,H0=T)))
}


real_analysis<-function(index_whatTo,candi_Anal,tt,n_core,nperm,phylum_list_perphy,phylums_candidates,pruned.tree2,taxo_sorted_df_table,sorted_df_table_formed,phenos_original,total.reads){
  #browser()
  whatToAnal<-candi_Anal[index_whatTo]
  real_output<-data.frame(matrix(nrow=length(phylums_candidates),ncol=length(whatToAnal)+2))
  names(real_output)<-c("Genus","num_species",whatToAnal)
  tt<-11
  if(str_detect(sessionInfo()[[1]]$os,"linux")){
    library(parallel)
    result_all_list<-mclapply(1:length(phylums_candidates),function(ww){
      subtree_tmp<-gettheSubtree(ww,phylum_list_perphy,phylums_candidates,pruned.tree2,taxo_sorted_df_table=taxo_sorted_df_table)
      if(!is.rooted(subtree_tmp[[1]])){
        subtree_tmp[[1]]<-root(subtree_tmp[[1]], 1, r = TRUE)
      }
      real_output[ww,0:length(whatToAnal)+2]<-real_0921(phenos_original,subtree_tmp,sorted_df_table_formed,nperm=nperm,type=tt,candi_Anal = candi_Anal,index_whatTo = index_whatTo,include_Toprank = T,total.reads = total.reads)
    },mc.cores=n_core)
  }else{
    library(foreach)
    require(doSNOW)
    getDoParWorkers()
    getDoParName()
    cl<-makeCluster(4, type = "SOCK")
    registerDoSNOW(cl)
    result_all_list<-foreach(ww = 1:length(phylums_candidates)) %dopar% {
      library(GUniFrac)
      library(ape);
      library(stringr);
      library(edgeR);
      library(MiSPU)
      library(MiRKAT)
      library(OMiAT)
      library(ecodist)
      source("./code/FunctionsTMAT/functions_tree_down.R",encoding = "UTF-8")
      source("./code/FunctionsTMAT/stats_and_funcs.R",encoding = "UTF-8")

      subtree_tmp<-gettheSubtree(ww,phylum_list_perphy,phylums_candidates,pruned.tree2,taxo_sorted_df_table=taxo_sorted_df_table)
      if(!is.rooted(subtree_tmp[[1]])){
        subtree_tmp[[1]]<-root(subtree_tmp[[1]], 1, r = TRUE)
      }
      real_output[ww,0:length(whatToAnal)+2]<-real_0921(phenos_original,subtree_tmp,sorted_df_table_formed,nperm=nperm,type=tt,candi_Anal = candi_Anal,index_whatTo = index_whatTo,include_Toprank = T,total.reads = total.reads)
    }
    stopCluster(cl)
  }
  real_output[,0:length(whatToAnal)+2]<-do.call("rbind",result_all_list)
  real_output[,1]<-phylums_candidates
  real_output[,2]<-howmany_list
  adj_TMAT<-p.adjust(real_output$Min_p_chi,method="fdr")
  adj_wil<-p.adjust(real_output$wilcoxon_minp,method="fdr")
  adj_omiat<-p.adjust(real_output$OMiAT,method="fdr")
  real_output<-cbind(real_output,adj_TMAT,adj_omiat,adj_wil)

  temptax_rank<-"D_1"
  patten_find<-paste0("(?<=",temptax_rank,"__).+?(?=;)")
  phylum_taxonomy<-str_extract(taxo_sorted_df_table[,2][match(phylums_candidates,phylum_list)],patten_find)
  output_real<-cbind(phylum_taxonomy,real_output)
  write.csv(output_real,paste0(filedir,"output1218_",filter1,"_",filter2,"_",nperm,"_",tax_rank,".csv"))
  return(output_real)
}


draw_pic<-function(ww){
  #browser()
  subtree_tmp<-gettheSubtree(ww,phylum_list_perphy,phylums_candidates,pruned.tree2,taxo_sorted_df_table=taxo_sorted_df_table)
  if(!is.rooted(subtree_tmp[[1]])){
    subtree_tmp[[1]]<-root(subtree_tmp[[1]], 1, r = TRUE)
  }
  # real_0921(phenos_original,subtree_tmp,sorted_df_table_formed,nperm=nperm,type=tt,candi_Anal = candi_Anal,index_whatTo = index_whatTo,include_Toprank = T,total.reads = total.reads)
  output<-real_1129_findWhere(phenos_original,subtree_tmp,sorted_df_table_formed,nperm=2,type=11,candi_Anal = candi_Anal,index_whatTo = index_whatTo,include_Toprank = T,total.reads = total.reads)
  # output$testnodes[[1]]
  tree.Groups_toplot<-subtree_tmp[[1]]
  tree.Groups_toplot$tip.label<-str_extract(subtree_tmp[[2]]$full,"(?<=;)D_6__[\\w .\\d]+$")
  colpick3<-rainbow(3)
  nlabel<- vector("character",Nnode(tree.Groups_toplot))
  npch<- vector("character",Nnode(tree.Groups_toplot))  

  #colpick3<-brewer.pal(3, "Spectral")
  if(which.min(output$min_index)==length(output$min_index)){
    #no<-1
    #npch[no]<-"★"
    l_pch<-""
    txt<-paste0("Test of Whole Genus Group with minimum p-value: ",round(output$min_index[which.min(output$min_index)]*10000)/10000)
    
  }else{
    target<-output$testnodes[[which.min(output$min_index)]]
    no<-target$root-Ntip(tree.Groups_toplot)
    #no<-1
    npch[no]<-"★"
    l_pch<-"★"
    txt<-paste0("TestNode with minimum p-value: ",round(output$min_index[which.min(output$min_index)]*10000)/10000)        
  }
  
  layout(matrix(c(1,2), ncol=1, byrow=TRUE), heights=c(1, 9))
  
  #par(mai=rep(0.4, 4))
  par(mai=c(0.1,0.1,0.4,0.1))
  plot.new()
  title(main = paste0("Genus: ",phylums_candidates[ww]),cex.main = 1.5)
  legend("left",pch=l_pch,col=colpick3[2],txt, cex=1.2,bty="n" )
  if(which.min(output$min_index)==length(output$min_index)){
    plot(tree.Groups_toplot)
    mycol<-rep(colpick3[c(1,3)][output$dir[which.min(output$min_index)]+1],length(target$group1))
  }else{
    
    plot(tree.Groups_toplot)
    mycol<-rep("black",length(target$group1))
    mycol[target$group1]<-c(colpick3[c(1,3)])[(output$dir[which.min(output$min_index)]+1)]
    mycol[target$group2]<-setdiff(colpick3[c(1,3)],colpick3[c(1,3)][(output$dir[which.min(output$min_index)]+1)])
  }

  tiplabels(paste0(tree.Groups_toplot$tip.label," (",subtree_tmp[[2]]$code,")"),pch="", col="white", adj=0, bg=mycol, cex=1.3)
  #nlabel[no]<-txt
  nodelabels(nlabel,pch=npch, col=colpick3[2], bg="white", cex=1.5,frame="none")
}

