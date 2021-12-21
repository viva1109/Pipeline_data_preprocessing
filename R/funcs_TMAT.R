inv_norm_trans<-function(x){
  qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x))) 
}
fisher_pval<-function(pval_1){
  len_tt<-length(pval_1)
  tt_1<- -2*sum(log(pval_1))
  1-pchisq(tt_1,df=2*len_tt)
}
make_data_TMAT<-function(target,type){
  cnt_table<<-plyr::count(target$tax_upper)
  genus<<-as.character(cnt_table[,1])
  lapply(1:length(genus),function(ind_genus){
    ind_otu_go<-which(target$tax_upper==genus[ind_genus])
    dSet<-matrix(target$dataset[,ind_otu_go],ncol=length(ind_otu_go))
    colnames(dSet)<-colnames(target$dataset)[ind_otu_go]
    simulSet_by_methods_fixed(y=target$groups,dataSet=dSet,ttr=target$totalcounts,type=type,method="TMAT")
  })
}
getCores<-function(dm,y,cov,conti,contrast,DetailResult=FALSE){
  # browser()
  betahat<-get_betahat(y,dm,cov=cov,conti=conti)
  n_indv<-length(y)
  if(conti){
    q<-1
  }else{
    q<-length(unique(y))-1
  }
  if(!is.null(cov)){
    t_re_1<-getT2(y,dm,cov=cov,conti=conti,DetailResult =DetailResult)
    p<-dim(cov)[2]
  }else{
    t_re_1<-getT2(y,dm,conti=conti,DetailResult =DetailResult)
    p<-1
  }  
  meanre<-tapply(dm,y,mean)
  # t_re_2<-getT2(y,dm[,2],conti=conti)
  if(DetailResult){
    result_bf<-c(t_re_1$Summary$Statistic[1],1-pchisq(t_re_1$Summary$Statistic[1],df=q),1-pf(t_re_1$Summary$Statistic[1],df1 = q,df2 = n_indv-p),meanre[contrast[1]]<meanre[contrast[2]],betahat)
    names(result_bf)<-c("statistic","p.chi","p.F","direction","betahat")
    result<-list(Node_summary=result_bf,Summary_Detail=t_re_1)
  }else{
    result<-c(t_re_1,1-pchisq(t_re_1,df=q),1-pf(t_re_1,df1 = q,df2 = n_indv-p),meanre[contrast[1]]<meanre[contrast[2]],betahat)
    names(result)<-c("statistic","p.chi","p.F","direction","betahat")  
  }
  return(result)
}

TMAT_func_dm<-function(dataSet,indi_tAnal,nperm=10000,cov=NULL,total.reads=NULL,type=NULL,include_Toprank=TRUE,permu=FALSE,conti=TRUE,multi=FALSE,contrast=c(1,2),getBeta=FALSE,getNodeStat=FALSE){
  # browser()
  n_indv<-dim(dataSet$simData$X_P)[1]
  if(is.null(cov)){
    p<-1
  }else{
    cov<-cbind(1,cov)
    p<-dim(cov)[2]
  }
  
  if(length(unique(dataSet$y[!is.na(dataSet$y)]))>2){
    if(conti){
      qq<-1
    }else{
      qq<-length(unique(dataSet$y[!is.na(dataSet$y)]))-1
    }
  }else{
    qq<-1
  }
  betahat<-NULL
  if(permu==T){
    permu_y<-lapply(1:nperm,function(x){
      sample(dataSet$y)
    })
    if(dim(dataSet$simData$X_P)[2]==1){
      lili<-NULL
      C1<-apply(dataSet$simData$X_P,1,sum)
      C2<-dataSet$simData$X_P_comp
      dm<- type_to_dm(C1,C2,type,total.reads)
      # genus_t<-getT2(dataSet$y,dm[,1],cov=cov)
      t_H0<-sapply(permu_y,function(py){
        genus_t<-getT2(py,dm[,1],cov=cov,conti=conti)
      })
      genus_t_real<-getT2(dataSet$y,dm[,1],cov=cov,conti=conti)
      pval_1<-get_permu_p(genus_t_real[1,1],t_H0,type="upper")
      meanre<-mean_group<-tapply(dm[,1],dataSet$y,mean)
      direction<-meanre[contrast[1]]<meanre[contrast[2]]
      out<-c(pval_1,pval_1,pval_1,pval_1)
    }else{
      t_H0<-t(sapply(permu_y,function(py){
        lili<-call_whattodo_v2(indi_tAnal$subtree_tmp[[1]])
        tp_pv<- get_PMT_out_new0918(lili,py,dataSet$simData$X_P,NULL,t_type=type,cov=cov,total.reads=total.reads,conti=conti)
        pval_1<-tp_pv[1,]
        pval_2<-tp_pv[2,]
        if(include_Toprank){
          C1<-apply(dataSet$simData$X_P,1,sum)
          C2<-dataSet$simData$X_P_comp
          dm<- type_to_dm(C1,C2,type,total.reads)
          genus_t<-getT2(py,dm[,1],cov=cov,conti=conti)
          genus_pval<-1-pchisq(genus_t,df=qq)
          pval_1<-c(pval_1,genus_pval)
          pval_2<-c(pval_2,1-pf(genus_t,df1 = qq,df2 = n_indv-p))
        }
        pval_minp_ch<-fisher_pval(pval_1)
        pval_minp_f<-fisher_pval(pval_2)
        out<-c(pval_minp_ch,pval_minp_f)
      }))
      
      lili<-call_whattodo_v2(indi_tAnal$subtree_tmp[[1]])
      tp_pv<- get_PMT_out_new0918(lili,dataSet$y,dataSet$simData$X_P,NULL,t_type=type,cov=cov,total.reads=total.reads)
      pval_1<-tp_pv[1,]
      pval_2<-tp_pv[2,]
      direction<-tp_pv[3,]
      
      if(include_Toprank){
        C1<-apply(dataSet$simData$X_P,1,sum)
        C2<-dataSet$simData$X_P_comp
        dm<- type_to_dm(C1,C2,type,total.reads)
        genus_t<-getT2(dataSet$y,dm[,1],cov=cov,conti=conti)
        genus_pval<-1-pchisq(genus_t,df=qq)
        pval_1<-c(pval_1,genus_pval)
        pval_2<-c(pval_2,1-pf(genus_t,df1 = qq,df2 = n_indv-p))
        genus_meanre<-mean_group<-tapply(dm[,1],dataSet$y,mean)
        direction<-c(direction,genus_meanre[contrast[1]]<genus_meanre[contrast[2]])
      }
      pval_minp_ch<-fisher_pval(pval_1)
      pval_minp_f<-fisher_pval(pval_2)
      out_bf<-c(pval_minp_ch,pval_minp_f)
      out<-c(get_permu_p(out_bf[1],t_H0[,1]),get_permu_p(out_bf[2],t_H0[,2]),get_permu_p(out_bf[3],t_H0[,3]),get_permu_p(out_bf[4],t_H0[,4]))
    }
  }else{
    
    if(dim(dataSet$simData$X_P)[2]==1){
      lili<-NULL
      C1<-apply(dataSet$simData$X_P,1,sum)
      C2<-dataSet$simData$X_P_comp
      # ddm<- type_to_dm(C1,C2,type,total.reads)
      tmp_dm<-type_to_dm(C1,C2,type,total.reads)  
      if(getNodeStat){
        if(!is.null(cov)){
          t_re_1<-getT2(dataSet$y,tmp_dm[,1],multi = multi,cov=cov,conti=conti)
          p<-dim(cov)[2]
        }else{
          t_re_1<-getT2(dataSet$y,tmp_dm[,1],multi = multi,conti=conti)
          p<-1
        }
        output<-c(t_re_1,p)
      }else{
        output<-tmp_dm
      }
      
    }else{
      lili<-call_whattodo_v2(indi_tAnal$subtree_tmp[[1]])
      if(getNodeStat){
        output<- get_PMT_out_new0918(lili,dataSet$y,dataSet$simData$X_P,NULL,t_type=type,cov=cov,total.reads=total.reads,contrast=contrast,getBeta=getBeta,getNodeStat = T)        
      }else{
        output<- get_PMT_out_new0918(lili,dataSet$y,dataSet$simData$X_P,NULL,t_type=type,cov=cov,total.reads=total.reads,contrast=contrast,getBeta=getBeta,getdm=T)
      }
      
      
      if(include_Toprank){
        C1<-apply(dataSet$simData$X_P,1,sum)
        C2<-dataSet$simData$X_P_comp
        tmp_dm2<-type_to_dm(C1,C2,type,total.reads)  
        if(getNodeStat){
          if(!is.null(cov)){
            t_re_2<-getT2(dataSet$y,tmp_dm2[,1],multi = multi,cov=cov,conti=conti)
            p2<-dim(cov)[2]
          }else{
            t_re_2<-getT2(dataSet$y,tmp_dm2[,1],multi = multi,conti=conti)
            p2<-1
          }
          output<-cbind(output,c(t_re_2,p2))
        }else{
          output<-cbind(output,dm)
        }
        
        
      }

    }
  }
  return(output)
}



TMAT_func_dm2<-function(dataSet,indi_tAnal,nperm=10000,cov=NULL,total.reads=NULL,type=NULL,include_Toprank=TRUE,permu=FALSE,conti=TRUE,contrast=c(1,2),getBeta=FALSE,DetailResult=FALSE){
    #browser()
  n_indv<-dim(dataSet$simData$X_P)[1]
  if(is.null(cov)){
    p<-1
  }else{
    cov<-cbind(1,cov)
    p<-dim(cov)[2]
  }
  if(length(unique(dataSet$y))>2){
    if(conti){
      qq<-1
    }else{
      qq<-length(unique(dataSet$y))-1
    }
  }else{
    qq<-1
  }
  betahat<-NULL
  
  if(dim(dataSet$simData$X_P)[2]==1){
    lili<-NULL
    C1<-apply(dataSet$simData$X_P,1,sum)
    C2<-dataSet$simData$X_P_comp
    dm<- type_to_dm(C1,C2,type,total.reads)
    # print("dm ,y")
    # print(dm[,1])
    # print(dataSet$y)
    
    if(DetailResult){
      list_cores_detail<-list(getCores(dm[,1],dataSet$y,cov=cov,conti=conti,DetailResult=DetailResult,contrast=contrast))
      cores<-list_cores_detail[[1]]$Node_summary
    }else{
      cores<-getCores(dm[,1],dataSet$y,cov=cov,conti=conti,DetailResult=DetailResult,contrast=contrast)
    }
    # genus_t<-cores["statistic"]
    if(getBeta){
      betahat<-cores["betahat"]
    }
    pval_1<-cores["p.chi"]
    pval_2<-cores["p.F"]
    direction<-cores["direction"]
    out<-c(cores[c("p.F")],cores[c("p.chi")])
    out_chi<-c(cores[c("p.F")],cores[c("p.chi")])
    dms<-dm[,1]
    if(type==16){
    dms<-inv_norm_trans(C1)
    }else{
    dms<-C1#change for biascor
    }
    
    


  }else{
    lili<-call_whattodo_v2(indi_tAnal$subtree_tmp[[1]])
    
    # tp_pv<- get_PMT_out_new0918(lili,dataSet$y,dataSet$simData$X_P,NULL,t_type=type,cov=cov,total.reads=total.reads,contrast=contrast,getBeta=getBeta)  
    
    dms<-lapply(lili,get_DM_tree,dataSet$y,dataSet$simData$X_P,t_X_P_comp=NULL,type=type,cov=cov,total.reads=total.reads,conti=conti,contrast=contrast,getBeta=getBeta,getdm=getdm,getNodeStat=getNodeStat)
    if(include_Toprank){
      C1<-apply(dataSet$simData$X_P,1,sum)
      #change 20210429 to correct bias
      C2<-dataSet$simData$X_P_comp
      dm<- type_to_dm(C1,C2,type,total.reads)
      dms<-c(dms,list(dm[,1]))
      
      if(type==16){
		dms<-c(dms,list(inv_norm_trans(C1)))
      }else{
		dms<-c(dms,list(C1))#change for biascor
      }

    }
    
    
    if(DetailResult){
      list_cores_detail<-lapply(dms,getCores,dataSet$y,cov,conti,DetailResult,contrast=contrast)
      list_cores<-sapply(list_cores_detail,function(data){
        data$Node_summary
      })
    }else{
      list_cores<-sapply(dms,getCores,dataSet$y,cov,conti,DetailResult,contrast=contrast)
    }  

    pval_1<-list_cores["p.chi",]
    pval_2<-list_cores["p.F",]
    direction<-list_cores["direction",]
    if(getBeta){
      betahat<-list_cores["betahat",]
    }
    
    p1<-exact_pval(pval_2)    #f-test,min_p
    p2<-fisher_pval(pval_2) #f-test,fisher
    p3<-exact_pval(pval_1)    #chi-test,min_p
    p4<-fisher_pval(pval_1) #chi-test,fisher
    out<-c(p1,p3)
    out_chi<-c(p2,p4)
  }
   
  
  if(DetailResult){
    return(list(pval=out,chi_pval=out_chi,min_index=pval_1,testnodes=lili,dir=direction,betahat=betahat,Details=list_cores_detail))
  }else{
    return(list(pval=out,chi_pval=out_chi,min_index=pval_1,testnodes=lili,dir=direction,betahat=betahat,dms=dms))
    
  }  
  # }else{
  #   if(dim(dataSet$simData$X_P)[2]==1){
  #     lili<-NULL
  #     C1<-apply(dataSet$simData$X_P,1,sum)
  #     C2<-dataSet$simData$X_P_comp
  #     dm<- type_to_dm(C1,C2,type,total.reads)
  #     genus_t<-getT2(dataSet$y,dm[,1],cov=cov,conti=conti)  
  #     if(getBeta){
  #       betahat<-get_betahat(dataSet$y,dm[,1],cov=cov,conti=conti)
  #     }
  #     genus_pval<-1-pchisq(genus_t,df=qq)
  #     pval_1<-genus_pval
  #     pval_2<-1-pf(genus_t,df1 = qq,df2 = n_indv-p)
  #     meanre<-mean_group<-tapply(dm[,1],dataSet$y,mean)
  #     direction<-meanre[contrast[1]]<meanre[contrast[2]]
  #     out<-c(pval_1,pval_2)
  #   }else{
  #     lili<-call_whattodo_v2(indi_tAnal$subtree_tmp[[1]])
  #     tp_pv<- get_PMT_out_new0918(lili,dataSet$y,dataSet$simData$X_P,NULL,t_type=type,cov=cov,total.reads=total.reads,contrast=contrast,getBeta=getBeta)  
  #     pval_1<-tp_pv[1,]
  #     pval_2<-tp_pv[2,]
  #     direction<-tp_pv[3,]
  #     if(getBeta){
  #       betahat<-tp_pv[4,]
  #     }
  #     if(include_Toprank){
  #       C1<-apply(dataSet$simData$X_P,1,sum)
  #       C2<-dataSet$simData$X_P_comp
  #       dm<- type_to_dm(C1,C2,type,total.reads)
  #       genus_t<-getT2(dataSet$y,dm[,1],cov=cov,conti=conti)
  #       if(getBeta){
  #         genus_betahat<-get_betahat(dataSet$y,dm[,1],cov=cov,conti=conti)
  #         betahat<-c(betahat,genus_betahat)
  #       }
  #       genus_pval<-1-pchisq(genus_t,df=qq)
  #       pval_1<-c(pval_1,genus_pval)
  #       pval_2<-c(pval_2,1-pf(genus_t,df1 = qq,df2 = n_indv-p))
  #       meanre<-mean_group<-tapply(dm[,1],dataSet$y,mean)
  #       direction<-c(direction,meanre[contrast[1]]<meanre[contrast[2]])
  #       
  #     }
  #     pval_minp_ch<-exact_pval(pval_1)
  #     pval_minp_f<-exact_pval(pval_2)
  #     out<-c(pval_minp_ch,pval_minp_f)
  #   }
  #   return(list(pval=out,min_index=pval_2,testnodes=lili,dir=direction,betahat=betahat))
  # }
 
}


TMAT_func_dm_mTMAT<-function(dataSet,indi_tAnal,nperm=10000,cov=NULL,total.reads=NULL,type=NULL,include_Toprank=TRUE,permu=FALSE,conti=TRUE,contrast=c(1,2),getBeta=FALSE,DetailResult=FALSE){
    #browser()
  n_indv<-dim(dataSet$simData$X_P)[1]
  if(is.null(cov)){
    p<-1
  }else{
    cov<-cbind(1,cov)
    p<-dim(cov)[2]
  }
  if(length(unique(dataSet$y))>2){
    if(conti){
      qq<-1
    }else{
      qq<-length(unique(dataSet$y))-1
    }
  }else{
    qq<-1
  }
  betahat<-NULL
  
  if(dim(dataSet$simData$X_P)[2]==1){
    lili<-NULL
    C1<-apply(dataSet$simData$X_P,1,sum)
    C2<-dataSet$simData$X_P_comp
    dm<- type_to_dm(C1,C2,type,total.reads)
    # print("dm ,y")
    # print(dm[,1])
    # print(dataSet$y)
    
    if(DetailResult){
      list_cores_detail<-list(getCores(dm[,1],dataSet$y,cov=cov,conti=conti,DetailResult=DetailResult,contrast=contrast))
      cores<-list_cores_detail[[1]]$Node_summary
    }else{
      cores<-getCores(dm[,1],dataSet$y,cov=cov,conti=conti,DetailResult=DetailResult,contrast=contrast)
    }
    # genus_t<-cores["statistic"]
    if(getBeta){
      betahat<-cores["betahat"]
    }
    pval_1<-cores["p.chi"]
    pval_2<-cores["p.F"]
    direction<-cores["direction"]
    out<-c(cores[c("p.F")],cores[c("p.chi")])
    out_chi<-c(cores[c("p.F")],cores[c("p.chi")])
    dms<-dm[,1]
    #if(type==16){
    #dms<-inv_norm_trans(C1)
    #}else{
    #dms<-C1#change for biascor
    #}
    
    


  }else{
    lili<-call_whattodo_v2(indi_tAnal$subtree_tmp[[1]])
    
    # tp_pv<- get_PMT_out_new0918(lili,dataSet$y,dataSet$simData$X_P,NULL,t_type=type,cov=cov,total.reads=total.reads,contrast=contrast,getBeta=getBeta)  
    
    dms<-lapply(lili,get_DM_tree,dataSet$y,dataSet$simData$X_P,t_X_P_comp=NULL,type=type,cov=cov,total.reads=total.reads,conti=conti,contrast=contrast,getBeta=getBeta,getdm=getdm,getNodeStat=getNodeStat)
    if(include_Toprank){
      C1<-apply(dataSet$simData$X_P,1,sum)
      #change 20210429 to correct bias
      C2<-dataSet$simData$X_P_comp
      dm<- type_to_dm(C1,C2,type,total.reads)
      dms<-c(dms,list(dm[,1]))
      
      #if(type==16){
	#	dms<-c(dms,list(inv_norm_trans(C1)))
      #}else{
	#	dms<-c(dms,list(C1))#change for biascor
      #}

    }
    
    
    if(DetailResult){
      list_cores_detail<-lapply(dms,getCores,dataSet$y,cov,conti,DetailResult,contrast=contrast)
      list_cores<-sapply(list_cores_detail,function(data){
        data$Node_summary
      })
    }else{
      list_cores<-sapply(dms,getCores,dataSet$y,cov,conti,DetailResult,contrast=contrast)
    }  

    pval_1<-list_cores["p.chi",]
    pval_2<-list_cores["p.F",]
    direction<-list_cores["direction",]
    if(getBeta){
      betahat<-list_cores["betahat",]
    }
    
    p1<-exact_pval(pval_2)    #f-test,min_p
    p2<-fisher_pval(pval_2) #f-test,fisher
    p3<-exact_pval(pval_1)    #chi-test,min_p
    p4<-fisher_pval(pval_1) #chi-test,fisher
    out<-c(p1,p3)
    out_chi<-c(p2,p4)
  }
   
  
  if(DetailResult){
    return(list(pval=out,chi_pval=out_chi,min_index=pval_1,testnodes=lili,dir=direction,betahat=betahat,Details=list_cores_detail))
  }else{
    return(list(pval=out,chi_pval=out_chi,min_index=pval_1,testnodes=lili,dir=direction,betahat=betahat,dms=dms))
    
  }  
  # }else{
  #   if(dim(dataSet$simData$X_P)[2]==1){
  #     lili<-NULL
  #     C1<-apply(dataSet$simData$X_P,1,sum)
  #     C2<-dataSet$simData$X_P_comp
  #     dm<- type_to_dm(C1,C2,type,total.reads)
  #     genus_t<-getT2(dataSet$y,dm[,1],cov=cov,conti=conti)  
  #     if(getBeta){
  #       betahat<-get_betahat(dataSet$y,dm[,1],cov=cov,conti=conti)
  #     }
  #     genus_pval<-1-pchisq(genus_t,df=qq)
  #     pval_1<-genus_pval
  #     pval_2<-1-pf(genus_t,df1 = qq,df2 = n_indv-p)
  #     meanre<-mean_group<-tapply(dm[,1],dataSet$y,mean)
  #     direction<-meanre[contrast[1]]<meanre[contrast[2]]
  #     out<-c(pval_1,pval_2)
  #   }else{
  #     lili<-call_whattodo_v2(indi_tAnal$subtree_tmp[[1]])
  #     tp_pv<- get_PMT_out_new0918(lili,dataSet$y,dataSet$simData$X_P,NULL,t_type=type,cov=cov,total.reads=total.reads,contrast=contrast,getBeta=getBeta)  
  #     pval_1<-tp_pv[1,]
  #     pval_2<-tp_pv[2,]
  #     direction<-tp_pv[3,]
  #     if(getBeta){
  #       betahat<-tp_pv[4,]
  #     }
  #     if(include_Toprank){
  #       C1<-apply(dataSet$simData$X_P,1,sum)
  #       C2<-dataSet$simData$X_P_comp
  #       dm<- type_to_dm(C1,C2,type,total.reads)
  #       genus_t<-getT2(dataSet$y,dm[,1],cov=cov,conti=conti)
  #       if(getBeta){
  #         genus_betahat<-get_betahat(dataSet$y,dm[,1],cov=cov,conti=conti)
  #         betahat<-c(betahat,genus_betahat)
  #       }
  #       genus_pval<-1-pchisq(genus_t,df=qq)
  #       pval_1<-c(pval_1,genus_pval)
  #       pval_2<-c(pval_2,1-pf(genus_t,df1 = qq,df2 = n_indv-p))
  #       meanre<-mean_group<-tapply(dm[,1],dataSet$y,mean)
  #       direction<-c(direction,meanre[contrast[1]]<meanre[contrast[2]])
  #       
  #     }
  #     pval_minp_ch<-exact_pval(pval_1)
  #     pval_minp_f<-exact_pval(pval_2)
  #     out<-c(pval_minp_ch,pval_minp_f)
  #   }
  #   return(list(pval=out,min_index=pval_2,testnodes=lili,dir=direction,betahat=betahat))
  # }
 
}



TMAT_func<-function(dataSet,indi_tAnal,nperm=10000,cov=NULL,total.reads=NULL,type=NULL,include_Toprank=TRUE,permu=FALSE,conti=TRUE,contrast=c(1,2),getBeta=TRUE,DetailResult=FALSE){
    #browser()
  n_indv<-dim(dataSet$simData$X_P)[1]
  if(is.null(cov)){
    p<-1
  }else{
    cov<-cbind(1,cov)
    p<-dim(cov)[2]
  }
  if(length(unique(dataSet$y))>2){
    if(conti){
      qq<-1
    }else{
      qq<-length(unique(dataSet$y))-1
    }
  }else{
    qq<-1
  }
  betahat<-NULL
  
  if(dim(dataSet$simData$X_P)[2]==1){
    lili<-NULL
    C1<-apply(dataSet$simData$X_P,1,sum)
    C2<-dataSet$simData$X_P_comp
    dm<- type_to_dm(C1,C2,type,total.reads)
    # print("dm ,y")
    # print(dm[,1])
    # print(dataSet$y)
    
    if(DetailResult){
      list_cores_detail<-list(getCores(dm[,1],dataSet$y,cov=cov,conti=conti,DetailResult=DetailResult,contrast=contrast))
      cores<-list_cores_detail[[1]]$Node_summary
    }else{
      cores<-getCores(dm[,1],dataSet$y,cov=cov,conti=conti,DetailResult=DetailResult,contrast=contrast)
    }
    # genus_t<-cores["statistic"]
    if(getBeta){
      betahat<-cores["betahat"]
    }
    pval_1<-cores["p.chi"]
    pval_2<-cores["p.F"]
    direction<-cores["direction"]
    out<-c(cores[c("p.F")],cores[c("p.chi")])
    out_chi<-c(cores[c("p.F")],cores[c("p.chi")])
  }else{
    lili<-call_whattodo_v2(indi_tAnal$subtree_tmp[[1]])
    
    # tp_pv<- get_PMT_out_new0918(lili,dataSet$y,dataSet$simData$X_P,NULL,t_type=type,cov=cov,total.reads=total.reads,contrast=contrast,getBeta=getBeta)  
    
    dms<-lapply(lili,get_DM_tree,dataSet$y,dataSet$simData$X_P,t_X_P_comp=NULL,type=type,cov=cov,total.reads=total.reads,conti=conti,contrast=contrast,getBeta=getBeta,getdm=getdm,getNodeStat=getNodeStat)
    if(include_Toprank){
      C1<-apply(dataSet$simData$X_P,1,sum)
      C2<-dataSet$simData$X_P_comp
      dm<- type_to_dm(C1,C2,type,total.reads)
      dms<-c(dms,list(dm[,1]))
    }
    
    
    if(DetailResult){
      list_cores_detail<-lapply(dms,getCores,dataSet$y,cov,conti,DetailResult,contrast=contrast)
      list_cores<-sapply(list_cores_detail,function(data){
        data$Node_summary
      })
    }else{
      list_cores<-sapply(dms,getCores,dataSet$y,cov,conti,DetailResult,contrast=contrast)
    }  

    pval_1<-list_cores["p.chi",]
    pval_2<-list_cores["p.F",]
    direction<-list_cores["direction",]
    if(getBeta){
      betahat<-list_cores["betahat",]
    }
    
    p1<-exact_pval(pval_2)    #f-test,min_p
    p2<-fisher_pval(pval_2) #f-test,fisher
    p3<-exact_pval(pval_1)    #chi-test,min_p
    p4<-fisher_pval(pval_1) #chi-test,fisher
    out<-c(p1,p3) #min_p
    out_chi<-c(p2,p4) #fisher
  }
   
  
  if(DetailResult){
    return(list(pval=out,chi_pval=out_chi,min_index=pval_1,testnodes=lili,dir=direction,betahat=betahat,Details=list_cores_detail))
  }else{
    return(list(pval=out,chi_pval=out_chi,min_index=pval_1,testnodes=lili,dir=direction,betahat=betahat))
    
  }  
  # }else{
  #   if(dim(dataSet$simData$X_P)[2]==1){
  #     lili<-NULL
  #     C1<-apply(dataSet$simData$X_P,1,sum)
  #     C2<-dataSet$simData$X_P_comp
  #     dm<- type_to_dm(C1,C2,type,total.reads)
  #     genus_t<-getT2(dataSet$y,dm[,1],cov=cov,conti=conti)  
  #     if(getBeta){
  #       betahat<-get_betahat(dataSet$y,dm[,1],cov=cov,conti=conti)
  #     }
  #     genus_pval<-1-pchisq(genus_t,df=qq)
  #     pval_1<-genus_pval
  #     pval_2<-1-pf(genus_t,df1 = qq,df2 = n_indv-p)
  #     meanre<-mean_group<-tapply(dm[,1],dataSet$y,mean)
  #     direction<-meanre[contrast[1]]<meanre[contrast[2]]
  #     out<-c(pval_1,pval_2)
  #   }else{
  #     lili<-call_whattodo_v2(indi_tAnal$subtree_tmp[[1]])
  #     tp_pv<- get_PMT_out_new0918(lili,dataSet$y,dataSet$simData$X_P,NULL,t_type=type,cov=cov,total.reads=total.reads,contrast=contrast,getBeta=getBeta)  
  #     pval_1<-tp_pv[1,]
  #     pval_2<-tp_pv[2,]
  #     direction<-tp_pv[3,]
  #     if(getBeta){
  #       betahat<-tp_pv[4,]
  #     }
  #     if(include_Toprank){
  #       C1<-apply(dataSet$simData$X_P,1,sum)
  #       C2<-dataSet$simData$X_P_comp
  #       dm<- type_to_dm(C1,C2,type,total.reads)
  #       genus_t<-getT2(dataSet$y,dm[,1],cov=cov,conti=conti)
  #       if(getBeta){
  #         genus_betahat<-get_betahat(dataSet$y,dm[,1],cov=cov,conti=conti)
  #         betahat<-c(betahat,genus_betahat)
  #       }
  #       genus_pval<-1-pchisq(genus_t,df=qq)
  #       pval_1<-c(pval_1,genus_pval)
  #       pval_2<-c(pval_2,1-pf(genus_t,df1 = qq,df2 = n_indv-p))
  #       meanre<-mean_group<-tapply(dm[,1],dataSet$y,mean)
  #       direction<-c(direction,meanre[contrast[1]]<meanre[contrast[2]])
  #       
  #     }
  #     pval_minp_ch<-exact_pval(pval_1)
  #     pval_minp_f<-exact_pval(pval_2)
  #     out<-c(pval_minp_ch,pval_minp_f)
  #   }
  #   return(list(pval=out,min_index=pval_2,testnodes=lili,dir=direction,betahat=betahat))
  # }
 
}

Out_setTMAT<-function(dataSet,ttr,type){
#print(paste0("type",type))
#browser()
  # print("ttrcheck")
  # print(ttr)
      if(type%in%c(1,3,13)){
        X_P<-dataSet
        X_P_comp<-ttr-apply(dataSet,1,sum)
      }else if (type%in%c(2,4,5,9:11,14)){
        X_P<-t(cpm(t(dataSet),log=TRUE,lib.size = ttr))
        X_P_comp_bf<-ttr-apply(dataSet,1,sum)
        X_P_comp<-t(cpm(t(X_P_comp_bf),log=TRUE,lib.size = ttr))
      }else if (type%in%c(6:8,12)){
        X_P<-logcpm(dataSet)#be careful
      }else if (type%in%c(15:18)){
        X_P_comp_bf<-ttr-apply(dataSet,1,sum)
        dataSet_after<-cbind(dataSet,X_P_comp_bf)
        X_P_fbf0<-t(cpm(t(dataSet_after),log=TRUE,lib.size = ttr))
        X_P_bf<-	log(exp(X_P_fbf0)+1)
        X_P<-matrix(X_P_bf[,-dim(X_P_bf)[2]],ncol=dim(X_P_bf)[2]-1)
        colnames(X_P) <-colnames(dataSet)
        X_P_comp<-matrix(X_P_bf[,dim(X_P_bf)[2]],ncol=1)
      }else if (type%in%c(19)){
        X_P_comp_bf<-ttr-apply(dataSet,1,sum)
        dataSet_after<-cbind(dataSet,X_P_comp_bf)
        X_P_fbf0<-t(cpm(t(dataSet_after),log=TRUE,lib.size = ttr))
        X_P_bf<-X_P_fbf0+1-min(X_P_fbf0)
        X_P<-matrix(X_P_bf[,-dim(X_P_bf)[2]],ncol=dim(X_P_bf)[2]-1)
        colnames(X_P) <-colnames(dataSet)
        X_P_comp<-matrix(X_P_bf[,dim(X_P_bf)[2]],ncol=1)
      }
      outdata<-list(X_P=X_P,X_P_comp=X_P_comp)
      return(outdata)
}

simulSet_by_methods_fixed<-function(y,dataSet,otu.tab.rff=NULL,ttr,type,method){
    # browser()
    outdata<-NULL
    #wilco
    if(method=="wilcoxon"){
      outdata<-dataSet/ttr
    }
    #MirKAT
     if(method=="oMirkat" & !is.null(otu.tab.rff) ){
       outdata <- otu.tab.rff
     }
    
    #MiSPU
    if(method=="aMiSPU"){
      X = as.matrix(dataSet)
      ind_rowsumnot0<-rowSums(X)!=0
      X_fixed<-X[ind_rowsumnot0,]
      outdata<-list(X_fixed=X_fixed,ind_rowsumnot0=ind_rowsumnot0)
    }
    if(method=="OMiAT"){
      outdata<-dataSet
    }
    #TMAT
    if(str_detect(method,"TMAT")){
    	outdata<-Out_setTMAT(dataSet=dataSet,ttr=ttr,type=type)
    }
    
    list(y=y,simData=outdata,ttr=ttr)
}

get_PMT_out_new0918<-function(lili,t_phenos,t_X_P,t_X_P_comp,t_type,cov=NULL,total.reads,conti=TRUE,contrast=c(1,2),getBeta=FALSE,getdm=FALSE,getNodeStat=FALSE){
  #browser()
  mm<-matrix(ncol=length(lili),nrow=2)
  tp_pv<-sapply(lili,get_PMT_new_F0918,t_phenos,t_X_P,t_X_P_comp=t_X_P_comp,type=t_type,cov=cov,total.reads=total.reads,conti=conti,contrast=contrast,getBeta=getBeta,getdm=getdm,getNodeStat=getNodeStat)
  mm<-tp_pv
  return(mm)
}




get_PMT_new_F0918<-function(data_lili,t_phenos,t_X_P,type,statistic=FALSE,multi=FALSE,t_X_P_comp=NULL,cov=NULL,total.reads,conti=TRUE,contrast=c(1,2),getBeta=FALSE,getdm=FALSE,getNodeStat=FALSE){
  # browser()
  if(type==1){
    ind_remain<- (t_X_P%*%as.matrix((data_lili$group1+data_lili$group2),ncol=1))!=0
  }else{
    ind_remain<-1:dim(t_X_P)[1]
  }
  X_P_remain<-t_X_P[ind_remain,]

  y<-t_phenos[ind_remain]
  n_indv<-length(y)
  #yAy <- t(y)%*%A%*%y
  C1<-X_P_remain%*%data_lili$group1
  C2<-X_P_remain%*%data_lili$group2
  dm<- type_to_dm(C1=C1,C2=C2,type=type,total.reads=total.reads)    
  if(getBeta){
	  betahat<-get_betahat(y,dm[,1],cov=cov,conti=conti)
  }else{
	  betahat<-NULL
  }
  if(!is.null(cov)){
     t_re_1<-getT2(y,dm[,1],multi = multi,cov=cov,conti=conti)
     p<-dim(cov)[2]
  }else{
     t_re_1<-getT2(y,dm[,1],multi = multi,conti=conti)
     p<-1
  }  
  meanre<-tapply(dm[,1],y,mean)
  # t_re_2<-getT2(y,dm[,2],conti=conti)
   if(getdm){
	   return(dm[,1])
   }else if(getNodeStat){
	   return(c(t_re_1,p))
   }else{
	  if(statistic){
	    return(c(t_re_1,t_re_1))
	  }else{
	    return(c(1-pchisq(t_re_1,df=1),1-pf(t_re_1,df1 = 1,df2 = n_indv-p),meanre[contrast[1]]<meanre[contrast[2]],betahat))
	  }
   }

  
  
}


get_DM_tree<-function(data_lili,t_phenos,t_X_P,type,statistic=FALSE,multi=FALSE,t_X_P_comp=NULL,cov=NULL,total.reads,conti=TRUE,contrast=c(1,2),getBeta=FALSE,getdm=FALSE,getNodeStat=FALSE){
  # browser()
  if(type==1){
    ind_remain<- (t_X_P%*%as.matrix((data_lili$group1+data_lili$group2),ncol=1))!=0
  }else{
    ind_remain<-1:dim(t_X_P)[1]
  }
  X_P_remain<-t_X_P[ind_remain,]
  
  y<-t_phenos[ind_remain]
  n_indv<-length(y)
  #yAy <- t(y)%*%A%*%y
  C1<-X_P_remain%*%data_lili$group1
  C2<-X_P_remain%*%data_lili$group2
  dm<- type_to_dm(C1=C1,C2=C2,type=type,total.reads=total.reads)  
  dm[,1]
}

type_to_dm<-function(C1,C2,type,total.reads,inv_norm_tr=FALSE){
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
      }else if (type%in%c(12,14)){
        dm[,1]<-C1-C2
        dm[,2]<-C2-C1
        # dm[,1]<-(C1/(C1+C2))
        # dm[,2]<-(C2/(C1+C2))
      }else if (type%in%c(13)){
        Cout<-cbind(C1,C2)
        dm<-t(cpm(t(Cout),lib.size=total.reads,log=TRUE))
      }else if (type%in%c(15,17)){
        dm[,1]<-log((C1)/(C2))
        dm[,2]<-log((C2)/(C1))
        # dm[,1]<-(C1/(C1+C2))
        # dm[,2]<-(C2/(C1+C2))
      }else if (type%in%c(16,18)){
        dm[,1]<-log((C1)/(C2))
        dm[,2]<-log((C2)/(C1))
        dm<-apply(dm,2,inv_norm_trans)
        # dm[,1]<-(C1/(C1+C2))
        # dm[,2]<-(C2/(C1+C2))
      }
      # if(inv_norm_tr){
      # print("invnorm")
      #dm<-apply(dm,2,inv_norm_trans)
      # }
      
      return(dm)
}


getT2<-function(t_y_bf,t_X,multi=FALSE,cov=NULL,conti=TRUE,DetailResult=FALSE){
  # browser()
  N<- length(t_y_bf)
  
  v_one_n<-rep(1,N)
  if(!is.null(cov)){
    t_A <- diag(N) - (cov %*% solve(t(cov)%*%cov) %*% t(cov))
    p<-dim(cov)[2]
  }else{
    t_A <- diag(N)-v_one_n%*%solve(t(v_one_n)%*%v_one_n)%*%t(v_one_n)
    p<-1
  }
  if(!is.null(cov)){
    ready_tt<-as.matrix(t(t_X)%*%t_A%*%t_X,ncol=1)
    t_SIGMA <- (ready_tt/(N-dim(cov)[2]))[1,1]
    # t_SIGMA <- cov(as.matrix(t_X,ncol=1))[1,1]
  }else{
    t_SIGMA <- cov(as.matrix(t_X,ncol=1))[1,1]
  }
  # browser()
  if(conti){
    t_y<-t_y_bf
    qq<-1
  }else{
    t_y<-c()
    if(length(unique(t_y_bf))>2){
      for( i in 1:length(unique(t_y_bf))){
        if(i!=1){
          t_y<-cbind(t_y,as.numeric(t_y_bf==unique(t_y_bf)[i]))
        }
      }
      qq<-dim(t_y)[2]
    }else{
      t_y<-t_y_bf
      qq<-1
    }
  }
  
  y_sq<-(t(t_y)%*%t_y)
  t_S <- t(t_y)%*%t_A%*%t_X
  t_yAy <- t(t_y)%*%t_A%*%t_y

  if(multi){
    if(!is.null(cov)){
      t_varS <- y_sq
    }else{
      t_varS <- t_yAy
    }
    t_varS <- t_SIGMA*t_yAy
    t_T<-t(t_S)%*%solve(t_varS)%*%t_S
    return(t_T)
  }else{
    t_varS <- t_SIGMA*t_yAy
    # print("t_varS")
    # print(t_varS)
    t_T<-t(t_S)%*%solve(t_varS)%*%t_S
    betahat<-solve(y_sq)%*%t_S

    if(DetailResult){
      if(length(unique(t_y_bf))>2){
        y_tmp<-unique(t_y_bf)
        Contrast<-paste0("[Y = ",y_tmp[-1]," <-> ", y_tmp[1],"]")
      }else{
        Contrast<-paste0("[Y = 1 <-> 0]")
      }
      Contrast<-c("F-test","Chi-squared Test",Contrast)
      
      row1<-c(NA,t_T,qq,N-p,1-pf(t_T,df1 = qq,df2 = N-p))
      row2<-c(NA,t_T,qq,NA,1-pchisq(t_T,df=qq))
      outSummary<-rbind(row1,row2,cbind(betahat,NA,NA,NA,NA))
      rownames(outSummary)<-Contrast
      
      colnames(outSummary)<-c("Estimate","Statistic","df1","df2","p.value")
      outSummary<-data.frame(outSummary)
      
      t_T_betahat_list<-lapply(1:dim(cov)[2],function(j){
        # browser()
        G<-matrix(cov[,j],ncol=length(j))
        Zg<-matrix(cov[,-j],ncol=dim(cov)[2]-length(j))
        H<-cbind(Zg,t_y)
        t_A <- diag(N) - (H %*% solve(t(H)%*%H) %*% t(H))
        ready_tt<-as.matrix(t(t_X)%*%t_A%*%t_X,ncol=1)
        t_SIGMA <- (ready_tt/(N-dim(cov)[2]))[1,1]
        # t_SIGMA <- cov(as.matrix(t_X,ncol=1))[1,1]
        t_S <- t(G)%*%t_A%*%t_X
        t_yAy <- t(G)%*%t_A%*%G
        y_sq<-(t(G)%*%G)
        t_varS <-  t_SIGMA*t_yAy
        tt_T<-t(t_S)%*%solve(t_varS)%*%t_S
        betahat<-solve(y_sq)%*%t_S
        if(length(unique(t_y_bf))>2){
          y_tmp<-unique(t_y_bf)
          Contrast<-paste0("[Y = ",y_tmp[-1]," <-> ", y_tmp[1],"]")
        }else{
          Contrast<-paste0("[Y = 1 <-> 0]")
        }
        Contrast<-c("F-test","Chi-squared Test",Contrast)
        
        row1<-c(NA,tt_T,qq,N-p,1-pf(tt_T,df1 = qq,df2 = N-p))
        row2<-c(NA,tt_T,qq,NA,1-pchisq(tt_T,df=qq))
        
        outSummary<-rbind(row1,row2,cbind(betahat,NA,NA,NA,NA))
        rownames(outSummary)<-Contrast
        
        colnames(outSummary)<-c("Estimate","Statistic","df1","df2","p.value")
        outSummary<-data.frame(outSummary)
        return(outSummary)
      })
      return(list(Summary=outSummary,Covariates_Summary=t_T_betahat_list))
      
    }else{

      return(t_T)
    }
  }
  
  

}


getT3<-function(t_y_bf,t_X,multi=FALSE,cov=NULL,conti=TRUE,DetailResult=FALSE){
  
  #browser()
  N<- length(t_y_bf)
  
  v_one_n<-rep(1,N)
  if(!is.null(cov)){
    t_A <- diag(N) - (cov %*% solve(t(cov)%*%cov) %*% t(cov))
    p<-dim(cov)[2]
  }else{
    t_A <- diag(N)-v_one_n%*%solve(t(v_one_n)%*%v_one_n)%*%t(v_one_n)
    p<-1
  }
  if(!is.null(cov)){
    ready_tt<-as.matrix(t(t_X)%*%t_A%*%t_X,ncol=1)
    t_SIGMA <- (ready_tt/(N-dim(cov)[2]))[1,1]
    # t_SIGMA <- cov(as.matrix(t_X,ncol=1))[1,1]
  }else{
    t_SIGMA <- cov(as.matrix(t_X,ncol=1))[1,1]
  }
  # browser()
  if(conti){
    t_y<-t_y_bf
    qq<-1
  }else{
    t_y<-c()
    if(length(unique(t_y_bf))>2){
      for( i in 1:length(unique(t_y_bf))){
        if(i!=1){
          t_y<-cbind(t_y,as.numeric(t_y_bf==unique(t_y_bf)[i]))
        }
      }
      qq<-dim(t_y)[2]
    }else{
      t_y<-t_y_bf
      qq<-1
    }
  }
  
  y_sq<-(t(t_y)%*%t_y)
  t_S <- t(t_y)%*%t_A%*%t_X
  t_yAy <- t(t_y)%*%t_A%*%t_y
  if(multi){
    if(!is.null(cov)){
      t_varS <- y_sq
    }else{
      t_varS <- t_yAy
    }
    t_varS <- t_SIGMA*t_yAy
    t_T<-t(t_S)%*%solve(t_varS)%*%t_S
    return(t_T)
  }else{
    t_varS <- t_SIGMA*t_yAy
    t_T<-t(t_S)%*%solve(t_varS)%*%t_S
    betahat<-solve(y_sq)%*%t_S

    if(DetailResult){
      if(length(unique(t_y_bf))>2){
        y_tmp<-unique(t_y_bf)
        Contrast<-paste0("[Y = ",y_tmp[-1]," <-> ", y_tmp[1],"]")
      }else{
        Contrast<-paste0("[Y = 1 <-> 0]")
      }
      Contrast<-c("F-test","Chi-squared Test",Contrast)
      
      row1<-c(NA,t_T,qq,N-p,1-pf(t_T,df1 = qq,df2 = N-p))
      row2<-c(NA,t_T,qq,NA,1-pchisq(t_T,df=qq))
      outSummary<-rbind(row1,row2,cbind(betahat,NA,NA,NA,NA))
      rownames(outSummary)<-Contrast
      
      colnames(outSummary)<-c("Estimate","Statistic","df1","df2","p.value")
      outSummary<-data.frame(outSummary)
      
      t_T_betahat_list<-lapply(1:dim(cov)[2],function(j){
        # browser()
        G<-matrix(cov[,j],ncol=length(j))
        Zg<-matrix(cov[,-j],ncol=dim(cov)[2]-length(j))
        H<-cbind(Zg,t_y)
        t_A <- diag(N) - (H %*% solve(t(H)%*%H) %*% t(H))
        ready_tt<-as.matrix(t(t_X)%*%t_A%*%t_X,ncol=1)
        t_SIGMA <- (ready_tt/(N-dim(cov)[2]))[1,1]
        # t_SIGMA <- cov(as.matrix(t_X,ncol=1))[1,1]
        t_S <- t(G)%*%t_A%*%t_X
        t_yAy <- t(G)%*%t_A%*%G
        t_varS <-  t_SIGMA*t_yAy
        tt_T<-t(t_S)%*%solve(t_varS)%*%t_S
        betahat<-solve(y_sq)%*%t_S
        if(length(unique(t_y_bf))>2){
          y_tmp<-unique(t_y_bf)
          Contrast<-paste0("[Y = ",y_tmp[-1]," <-> ", y_tmp[1],"]")
        }else{
          Contrast<-paste0("[Y = 1 <-> 0]")
        }
        Contrast<-c("F-test","Chi-squared Test",Contrast)
        
        row1<-c(NA,tt_T,qq,N-p,1-pf(tt_T,df1 = qq,df2 = N-p))
        row2<-c(NA,tt_T,qq,NA,1-pchisq(tt_T,df=qq))
        
        outSummary<-rbind(row1,row2,cbind(betahat,NA,NA,NA,NA))
        rownames(outSummary)<-Contrast
        
        colnames(outSummary)<-c("Estimate","Statistic","df1","df2","p.value")
        outSummary<-data.frame(outSummary)
        return(outSummary)
      })
      return(list(Summary=outSummary,Covariates_Summary=t_T_betahat_list))
      
    }else{

      return(t_T)
    }
  }

}

exact_pval<-function(tmp_pval){
  1-(1-min(tmp_pval))^length(tmp_pval)
}

get_betahat<-function(t_y_bf,t_Z,cov=NULL,conti=TRUE){
  # browser()
  N<- length(t_y_bf)
  v_one_n<-rep(1,N)
  if(!is.null(cov)){
    t_A <- diag(N) - (cov %*% solve(t(cov)%*%cov) %*% t(cov))
  }else{
    t_A <- diag(N)-v_one_n%*%solve(t(v_one_n)%*%v_one_n)%*%t(v_one_n)
  }
    if(conti){
    t_y<-t_y_bf
  }else{
    t_y<-c()
    if(length(unique(t_y_bf))>2){
      for( i in 1:length(unique(t_y_bf))){
        if(i!=1){
          t_y<-cbind(t_y,as.numeric(t_y_bf==unique(t_y_bf)[i]))
        }
      }
    }else{
      t_y<-t_y_bf
    }
  }
  y_sq<-(t(t_y)%*%t_y)
  t_S <- t(t_y)%*%t_A%*%t_Z
  betahat<-solve(y_sq)%*%t_S
  return(betahat)
}

OMiAT<-function (Y, otu.tab, cov = NULL, tree, total.reads = NULL, model = c("gaussian", 
                                                                             "binomial"), pow = c(1:4, Inf), g.unif.alpha = c(0.5), n.perm = 5000) 
{
  if (length(Y) != nrow(otu.tab)) {
    otu.tab <- t(otu.tab)
  }
  if (is.null(total.reads)) {
    prop <- apply(t(apply(otu.tab, 1, function(x) x/sum(x))), 
                  2, scale)
  }
  else {
    prop <- otu.tab/total.reads
    prop <- apply(prop, 2, scale)
  }
  if (is.null(cov)) {
    r <- Y - mean(Y)
  }
  else {
    fit <- glm(Y ~ cov, family = model)
    res <- Y - fitted.values(fit)
    r <- res - mean(res)
  }
  r <- jitter(r)
  r.s<-lapply(1:n.perm,function(ind){sample(r)})
  U <- as.vector(t(prop) %*% r)
  Ts = rep(NA, length(pow))
  Ts[which(pow == Inf)] <- max(U)
  Ts[which(pow != Inf)] <- unlist(lapply(as.list(pow[which(pow != 
                                                             Inf)]), function(x) return(sum(U^x))))
  U0 <- lapply(r.s, function(x) return(as.vector(t(prop) %*% 
                                                   x)))
  T0s <- list()
  pvs <- rep(NA, length(pow))
  for (j in 1:length(pow)) {
    T0s[[j]] <- T0s.e <- sapply(U0, function(x) if (pow[j] < 
                                                    Inf) 
      return(sum(x^pow[j]))
      else return(max(x)))
    pvs[j] <- length(which(abs(T0s.e) > abs(Ts[j])))/n.perm
  }
  T.aspu <- min(pvs)
  T0.aspu<-sapply(1:n.perm,function(l){
    T0s.n<-lapply(T0s,function(data,tl){data[-tl]},l)
    a.Ts<-lapply(T0s,function(data,tl){data[tl]},l)
    a.pvs <- unlist(mapply(function(x, y) length(which(abs(x) >= 
                                                         abs(y)))/(n.perm - 1), T0s.n, a.Ts))
    return(min(a.pvs))
  })
  p.aspu <- length(which(T0.aspu < T.aspu))/n.perm
  Ts <- c(Ts, T.aspu)
  names(Ts) <- c(paste("SPU(", pow, ")", sep = ""), "aSPU")
  spu.pvs <- c(pvs, p.aspu)
  names(spu.pvs) <- c(paste("SPU(", pow, ")", sep = ""), "aSPU")
  if (is.null(total.reads)) {
    unifs <- GUniFrac(otu.tab, tree, alpha = c(g.unif.alpha, 
                                               1))$unifracs
    bray.curtis <- as.matrix(bcdist(otu.tab))
    u.unif <- unifs[, , "d_UW"]
    w.unif <- unifs[, , "d_1"]
    g.unif <- list()
  }
  else {
    unifs <- GUniFrac2(otu.tab, tree, alpha = c(g.unif.alpha, 
                                                1), total.reads = total.reads)$unifracs
    bray.curtis <- as.matrix(bcdist(otu.tab))
    u.unif <- unifs[, , "d_UW"]
    w.unif <- unifs[, , "d_1"]
    g.unif <- list()
  }
  for (k in 1:length(g.unif.alpha)) {
    g.unif[[k]] <- unifs[, , paste("d_", g.unif.alpha[k], 
                                   sep = "")]
  }
  if (sum(length(which(is.na(u.unif)))) > 0) {
    bray.curtis.kern <- D2K(bray.curtis)
    if (model == "gaussian") 
      Q <- as.numeric(t(r) %*% bray.curtis.kern %*% r)
    if (model == "binomial") 
      Q <- as.numeric(t(r) %*% bray.curtis.kern %*% r)
    Q0 <- rep(NA, n.perm)
    for (j in 1:n.perm) {
      if (model == "gaussian") 
        Q0[j] <- t(r.s[[j]]) %*% bray.curtis.kern %*% 
          r.s[[j]]
      if (model == "binomial") 
        Q0[j] <- t(r.s[[j]]) %*% bray.curtis.kern %*% 
          r.s[[j]]
    }
    Q.omni <- p.bc <- length(which(abs(Q0) > abs(Q)))/n.perm
    Qs <- c(Q, Q.omni)
    mirkat.pvs <- c(p.bc, p.bc)
    names(Qs) <- c("Bray-Curtis", "Opt.MiRKAT")
    names(mirkat.pvs) <- c("Bray-Curtis", "Opt.MiRKAT")
    M.omiat0 <- apply(cbind(T0.aspu, Q0), 1, min)
    M.omiat <- min(T.aspu, Q.omni)
    p.omiat <- length(which(M.omiat0 < M.omiat))/n.perm
    names(M.omiat) <- "OMiAT"
    names(p.omiat) <- "OMiAT"
  }
  else {
    bray.curtis.kern <- D2K(bray.curtis)
    u.unif.kern <- D2K(u.unif)
    w.unif.kern <- D2K(w.unif)
    g.unif.kern <- list()
    for (k in 1:length(g.unif.alpha)) {
      g.unif.kern[[k]] <- D2K(g.unif[[k]])
    }
    list.kernels <- c(list(bray.curtis.kern = bray.curtis.kern, 
                           u.unif.kern = u.unif.kern, w.unif.kern = w.unif.kern), 
                      g.unif.kern)
    Qs <- rep(NA, length(list.kernels))
    for (j in 1:length(list.kernels)) {
      if (model == "gaussian") 
        Qs[j] <- t(r) %*% list.kernels[[j]] %*% r
      if (model == "binomial") 
        Qs[j] <- t(r) %*% list.kernels[[j]] %*% r
    }
    Q0s <- list()
    for (j in 1:length(list.kernels)) {
      Q0s.inv <- rep(NA, n.perm)
      for (k in 1:n.perm) {
        if (model == "gaussian") 
          Q0s.inv[k] <- t(r.s[[k]]) %*% list.kernels[[j]] %*% 
            r.s[[k]]
        if (model == "binomial") 
          Q0s.inv[k] <- t(r.s[[k]]) %*% list.kernels[[j]] %*% 
            r.s[[k]]
      }
      Q0s[[j]] <- Q0s.inv
    }
    mirkat.pvs <- rep(NA, length(list.kernels))
    for (j in 1:length(list.kernels)) {
      mirkat.pvs[j] <- length(which(abs(Q0s[[j]]) > abs(Qs[[j]])))/n.perm
    }
    Q.omni <- min(mirkat.pvs)
    Q0.omni <- rep(NA, n.perm)
    for (l in 1:n.perm) {
      Q0s.n <- list()
      for (m in 1:length(list.kernels)) {
        Q0s.n[[m]] <- Q0s[[m]][-l]
      }
      if (model == "gaussian") 
        a.Qs <- unlist(lapply(list.kernels, function(x) return(t(r.s[[l]]) %*% 
                                                                 x %*% r.s[[l]])))
      if (model == "binomial") 
        a.Qs <- unlist(lapply(list.kernels, function(x) return(t(r.s[[l]]) %*% 
                                                                 x %*% r.s[[l]])))
      a.pvs <- unlist(mapply(function(x, y) length(which(abs(x) > 
                                                           abs(y)))/(n.perm - 1), Q0s.n, a.Qs))
      Q0.omni[l] <- min(a.pvs)
    }
    p.omni <- length(which(Q0.omni < Q.omni))/n.perm
    Qs <- c(Qs, Q.omni)
    names(Qs) <- c("Bray-Curtis", "U.UniFrac", "W.UniFrac", 
                   paste("G.UniFrac(", g.unif.alpha, ")", sep = ""), 
                   "Opt.MiRKAT")
    mirkat.pvs <- c(mirkat.pvs, p.omni)
    names(mirkat.pvs) <- c("Bray-Curtis", "U.UniFrac", "W.UniFrac", 
                           paste("G.UniFrac(", g.unif.alpha, ")", sep = ""), 
                           "Opt.MiRKAT")
    M.omiat0 <- apply(cbind(T0.aspu, Q0.omni), 1, min)
    M.omiat <- min(T.aspu, Q.omni)
    p.omiat <- length(which(M.omiat0 < M.omiat))/n.perm
    names(M.omiat) <- "OMiAT"
    names(p.omiat) <- "OMiAT"
  }
  return(list(SPU.pvs = spu.pvs, MiRKAT.pvs = mirkat.pvs, OMiAT.pvalue = p.omiat))
}