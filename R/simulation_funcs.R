permu_p_MiRKATmethod<-function(H1_p,H0_ps,lower=T){
	nperm_t<-length(H0_ps)
	if(lower){
		p_permu<-(sum(H0_ps<H1_p)+sum(H0_ps==H1_p)/2+1/2)/(nperm_t+1)
		#p_permu = base::rank(p_added)[1]/(nperm_t + 1)
	}else{
		#p_permu = 1-(base::rank(p_added)[1]-1)/(nperm_t + 1)
		p_permu<-(sum(H0_ps>H1_p)+sum(H0_ps==H1_p)/2+1/2)/(nperm_t+1)
	}
}

lapply_to_mclapply<-function(input,func,chunk1=1200,n_core=1,method=NULL,indind=NULL,...){
		#browser()
		  len<-length(input)
		  chunk2<-ceiling(len/n_core)
		  chunk<-min(chunk1,chunk2)
		  J<-ceiling(len/chunk)
		  list_input<-vector("list",J)
		print("J")
		print(J)
		  for (j in 1:J){
			  list_input[[j]]<-input[(1+chunk*(j-1)):min(chunk*j,len)]
		  }
		  outputSet2<-mclapply(list_input,function(input_indiv){
		    #outputSet<-lapply(data[[t_ind]][(1+chunk*(j-1)):min(chunk*j,len)],func,indi_tAnal=indi_tAnal[[t_ind]],nperm=nperm,cov=cov,total.reads=total.reads,type=type,indind=t_ind)
		    outputSet<-  lapply(input_indiv,func,...)
		  },mc.cores=n_core)
		#print(outputSet2)
		print(indind)
		if(is.null(method)){
			print("method is null")
			saveRDS(outputSet2, paste0(paste0(filedir,head_label),"/rawP_NONAME_",indind,".rds"))
		}else{
		       saveRDS(outputSet2, paste0(paste0(filedir,head_label),"/rawP_",method,"_",indind,".rds"))
		}
		
		outputSet3<-do.call("c",outputSet2)
		rm(outputSet2)
		return(outputSet3)
}


split_save_mclapply<-function(input,label="nolabel",sublabel=1,func,chunk1=100,n_core=1,...){
		if(n_core!=1){
			chunk1<-n_core
		}
		  J<-ceiling(len/chunk1)
		  list_input<-vector("list",J)
		  system(paste0("mkdir ",paste0(filedir,head_label),"/",label,"_",sublabel2))
		  for (j in 1:J){
			  list_input[[j]]<-input[(1+chunk1*(j-1)):min(chunk1*j,len)]
			  result<-mclapply(list_input[[j]],func,mc.cores=n_core,...)
			  tmpChar<-formatC(j, width = nchar(J),flag = 0)
			  sublabel2<-formatC(sublabel, width = nchar(t_table2[1]),flag = 0)
			  saveRDS(result ,paste0(paste0(filedir,head_label),"/",label,"_",sublabel2,"/","Sim_",tmpChar,".rds"))
		  }
		return(NULL)
}

generateH1<-function(y,sel_data,sel_data_rf,simul=NULL,useMC=T){
  beta_per<-simul$t_table$beta_per_list
  if(as.logical(simul$t_table$plus_minus)){
	if(useMC){
		lapply_to_mclapply(simul$permu_indcausal_pheno,function(list,y,beta_per,sel_data,sel_data_rf){
		    y <- y[list$pheno_index]
		    ind_causal<-list$ind_causal
		    TF_plus<-sapply(ind_causal,function(ind){sample(c(T,F),1)})
		    ind_plus<- ind_causal[TF_plus]
		    ind_minus<-ind_causal[!TF_plus]		
		    if(length(ind_plus)==1){
		      plus_ap_sd<-sd(sel_data[,ind_plus])
		      plus_ap_sd_rf<-sd(sel_data_rf[,ind_plus])
		      sel_data[as.logical(y),ind_plus]<-round(t(t(sel_data[as.logical(y),ind_plus])+beta_per*plus_ap_sd))
		      sel_data_rf[as.logical(y),ind_plus]<-round(t(t(sel_data_rf[as.logical(y),ind_plus])+beta_per*plus_ap_sd_rf))
		    }else if(length(ind_plus)==0){

		    }else{
		      plus_ap_sd<-apply(sel_data[,ind_plus],2,sd)
		      plus_ap_sd_rf<-apply(sel_data_rf[,ind_plus],2,sd)
		      sel_data[as.logical(y),ind_plus]<-round(t(t(sel_data[as.logical(y),ind_plus])+beta_per*plus_ap_sd))
		      sel_data_rf[as.logical(y),ind_plus]<-round(t(t(sel_data_rf[as.logical(y),ind_plus])+beta_per*plus_ap_sd_rf))
		    }
		    if(length(ind_minus)==1){
		      minus_ap_sd<-sd(sel_data[,ind_minus])
		      minus_ap_sd_rf<-sd(sel_data_rf[,ind_minus])
		      sel_data[!as.logical(y),ind_minus]<-round(t(t(sel_data[!as.logical(y),ind_minus])+beta_per*minus_ap_sd))
		      sel_data_rf[!as.logical(y),ind_minus]<-round(t(t(sel_data_rf[!as.logical(y),ind_minus])+beta_per*minus_ap_sd_rf))
		    }else if(length(ind_minus)==0){

		    }else{
		      minus_ap_sd<-apply(sel_data[,ind_minus],2,sd)
		      minus_ap_sd_rf<-apply(sel_data_rf[,ind_minus],2,sd)
		      sel_data[!as.logical(y),ind_minus]<-round(t(t(sel_data[!as.logical(y),ind_minus])+beta_per*minus_ap_sd))
		      sel_data_rf[!as.logical(y),ind_minus]<-round(t(t(sel_data_rf[!as.logical(y),ind_minus])+beta_per*minus_ap_sd_rf))
		    }
		    return(list(y=y,sel_data=sel_data,sel_data_rf=sel_data_rf))
		  },n_core=n_core,y=y,beta_per=beta_per,sel_data=sel_data,sel_data_rf=sel_data_rf)
	}else{
	       	  lapply(simul$permu_indcausal_pheno,function(list,y,beta_per,sel_data,sel_data_rf){
		    y <- y[list$pheno_index]
		    ind_causal<-list$ind_causal
		    TF_plus<-sapply(ind_causal,function(ind){sample(c(T,F),1)})
		    ind_plus<- ind_causal[TF_plus]
		    ind_minus<-ind_causal[!TF_plus]		
		    if(length(ind_plus)==1){
		      plus_ap_sd<-sd(sel_data[,ind_plus])
		      plus_ap_sd_rf<-sd(sel_data_rf[,ind_plus])
		      sel_data[as.logical(y),ind_plus]<-round(t(t(sel_data[as.logical(y),ind_plus])+beta_per*plus_ap_sd))
		      sel_data_rf[as.logical(y),ind_plus]<-round(t(t(sel_data_rf[as.logical(y),ind_plus])+beta_per*plus_ap_sd_rf))
		    }else if(length(ind_plus)==0){

		    }else{
		      plus_ap_sd<-apply(sel_data[,ind_plus],2,sd)
		      plus_ap_sd_rf<-apply(sel_data_rf[,ind_plus],2,sd)
		      sel_data[as.logical(y),ind_plus]<-round(t(t(sel_data[as.logical(y),ind_plus])+beta_per*plus_ap_sd))
		      sel_data_rf[as.logical(y),ind_plus]<-round(t(t(sel_data_rf[as.logical(y),ind_plus])+beta_per*plus_ap_sd_rf))
		    }
		    if(length(ind_minus)==1){
		      minus_ap_sd<-sd(sel_data[,ind_minus])
		      minus_ap_sd_rf<-sd(sel_data_rf[,ind_minus])
		      sel_data[!as.logical(y),ind_minus]<-round(t(t(sel_data[!as.logical(y),ind_minus])+beta_per*minus_ap_sd))
		      sel_data_rf[!as.logical(y),ind_minus]<-round(t(t(sel_data_rf[!as.logical(y),ind_minus])+beta_per*minus_ap_sd_rf))
		    }else if(length(ind_minus)==0){

		    }else{
		      minus_ap_sd<-apply(sel_data[,ind_minus],2,sd)
		      minus_ap_sd_rf<-apply(sel_data_rf[,ind_minus],2,sd)
		      sel_data[!as.logical(y),ind_minus]<-round(t(t(sel_data[!as.logical(y),ind_minus])+beta_per*minus_ap_sd))
		      sel_data_rf[!as.logical(y),ind_minus]<-round(t(t(sel_data_rf[!as.logical(y),ind_minus])+beta_per*minus_ap_sd_rf))
		    }

		    return(list(y=y,sel_data=sel_data,sel_data_rf=sel_data_rf))
		  },y=y,beta_per=beta_per,sel_data=sel_data,sel_data_rf=sel_data_rf)
	}
	  

  }else{
  	if(useMC){
		lapply_to_mclapply(simul$permu_indcausal_pheno,function(list,y,beta_per,sel_data,sel_data_rf){
		    y <- y[list$pheno_index]
		    ind_causal<-list$ind_causal
		    if(length(ind_causal)==1){
		      ap_sd<-sd(sel_data[,ind_causal])
		      ap_sd_rf<-sd(sel_data_rf[,ind_causal])
		    }else{
		      ap_sd<-apply(sel_data[,ind_causal],2,sd)
		      ap_sd_rf<-apply(sel_data_rf[,ind_causal],2,sd)
		    }
		    sel_data[as.logical(y),ind_causal]<-round(t(t(sel_data[as.logical(y),ind_causal])+beta_per*ap_sd))
		    sel_data_rf[as.logical(y),ind_causal]<-round(t(t(sel_data_rf[as.logical(y),ind_causal])+beta_per*ap_sd_rf))
		    return(list(y=y,sel_data=sel_data,sel_data_rf=sel_data_rf))
		  },n_core=n_core,y=y,beta_per=beta_per,sel_data=sel_data,sel_data_rf=sel_data_rf)
	}else{
		  lapply(simul$permu_indcausal_pheno,function(list,y,beta_per,sel_data,sel_data_rf){
		    y <- y[list$pheno_index]
		    ind_causal<-list$ind_causal
		    if(length(ind_causal)==1){
		      ap_sd<-sd(sel_data[,ind_causal])
		      ap_sd_rf<-sd(sel_data_rf[,ind_causal])
		    }else{
		      ap_sd<-apply(sel_data[,ind_causal],2,sd)
		      ap_sd_rf<-apply(sel_data_rf[,ind_causal],2,sd)
		    }
		    sel_data[as.logical(y),ind_causal]<-round(t(t(sel_data[as.logical(y),ind_causal])+beta_per*ap_sd))
		    sel_data_rf[as.logical(y),ind_causal]<-round(t(t(sel_data_rf[as.logical(y),ind_causal])+beta_per*ap_sd_rf))
		    return(list(y=y,sel_data=sel_data,sel_data_rf=sel_data_rf))
		  },y=y,beta_per=beta_per,sel_data=sel_data,sel_data_rf=sel_data_rf)
	}
  } 
  
}

#here!!!!!!!!!!!!!!!



simulSet_by_methods<-function(list,X_P_ori,methods,totalreads,useMC=T){
#browser()
	if(useMC){
		#lapply_to_mclapply(list,function(list_argu){
		#	X_P_bf<-list_argu$sel_data
		#	X_P_bfC <- totalreads - apply(X_P_ori$nrfD,1,sum)
		#	ttr<- X_P_bfC+apply(X_P_bf,1,sum)
		#	lapply(methods,function(method){
		#		    if(str_detect(method,"TMAT")){
		#		      type_TMAT<-as.numeric(str_sub(method,5,-1))
		#		    }
		#		simulSet_by_methods_fixed(list_argu$y,list_argu$sel_data,list_argu$sel_data_rf,ttr=totalreads,type_TMAT,method)
		#	})
		#},n_core=n_core)
		lapply_to_mclapply(list,function(list_argu){
			X_P_bf<-list_argu$sel_data
			X_P_bfC <- totalreads - apply(X_P_ori$nrfD,1,sum)
			ttr<- X_P_bfC+apply(X_P_bf,1,sum)
			lapply(methods,function(method){
				    if(str_detect(method,"TMAT")){
				      type_TMAT<-as.numeric(str_sub(method,5,-1))
				    }
				simulSet_by_methods_fixed(list_argu$y,list_argu$sel_data,list_argu$sel_data_rf,ttr=totalreads,type_TMAT,method)
			})
		},n_core=n_core)
	}else{
		  lapply(list,function(list_argu){
		    X_P_bf<-list_argu$sel_data
		    X_P_bfC <- totalreads - apply(X_P_ori$nrfD,1,sum)
		    ttr<- X_P_bfC+apply(X_P_bf,1,sum)
		    #X_P_bf_rf <- list_argu$sel_data_rf
		    #X_P_bf_rfC <- total.reads.rf - apply(X_P_ori$rfD,1,sum)
		    #rfftable_bf<-cbind(X_P_bf_rf,X_P_bf_rfC)
		    #otu.tab.rff_togo <- Rarefy(rfftable_bf)$otu.tab.rff[,-dim(rfftable_bf)[2]]
		    #dataSet<-data.frame(matrix(X_P_bf3,ncol=dim(X_P_bf3)[2]))
		    #names(dataSet)<-colnames(X_P_bf3)
		    	lapply(methods,function(method){
				    if(str_detect(method,"TMAT")){
				      type_TMAT<-as.numeric(str_sub(method,5,-1))
				    }
				simulSet_by_methods_fixed(list_argu$y,list_argu$sel_data,list_argu$sel_data_rf,ttr=totalreads,type_TMAT,method)
			})
		  })
	}


}


statistic_Set<-function(t_ind,data,func,indi_tAnal,nperm=NULL,cov=NULL,total.reads=NULL,type=NULL,useMC=T,...){
	#browser()
		#print("t_ind")
		#print(t_ind)

	crit<-try(random_tree)

	
	if(useMC){
		lapply_to_mclapply(data[[t_ind]],function(tt_data,tindi_tAnal,...){
			if(crit==T){
				tindi_tAnal$subtree_tmp[[1]]<-getRandomTree(tindi_tAnal$subtree_tmp[[1]],factor_random)
			}
			func(tt_data,indi_tAnal=tindi_tAnal,indind=t_ind,...)
		},n_core=n_core,tindi_tAnal=indi_tAnal[[t_ind]],nperm=nperm,cov=cov,total.reads=total.reads,type=type,...)


		#chunk1<-1000
		#  len<-length(data[[t_ind]])
		#  chunk2<-ceiling(len/n_core)
		#  chunk<-min(chunk1,chunk2)
		#  J<-ceiling(len/chunk)
		#  outputSet2<-mclapply(1:J,function(j){
		#    #outputSet<-lapply(data[[t_ind]][(1+chunk*(j-1)):min(chunk*j,len)],func,indi_tAnal=indi_tAnal[[t_ind]],nperm=nperm,cov=cov,total.reads=total.reads,type=type,indind=t_ind)
		#    outputSet<-lapply(data[[t_ind]][(1+chunk*(j-1)):min(chunk*j,len)],function(tt_data,indi_tAnal,...){
		#				if(random_tree){
		#					indi_tAnal$subtree_tmp[[1]]<-getRandomTree(indi_tAnal$subtree_tmp[[1]],factor_random)
		#				}
		#				func(tt_data,indi_tAnal=indi_tAnal,...)
		#			},indi_tAnal=indi_tAnal[[t_ind]],nperm=nperm,cov=cov,total.reads=total.reads,type=type,indind=t_ind)
		#  },mc.cores=n_core)
		#print(outputSet2)
		#  outputSet3<-do.call("c",outputSet2)
		#  do.call("cbind",outputSet3)
	}else{
		
		 #sapply(data[[t_ind]],func,indi_tAnal=indi_tAnal[[t_ind]],nperm=nperm,cov=cov,total.reads=total.reads,type=type,indind=t_ind)
		 sapply(data[[t_ind]],function(tt_data,tindi_tAnal,...){
			if(crit==T){
				tindi_tAnal$subtree_tmp[[1]]<-getRandomTree(tindi_tAnal$subtree_tmp[[1]],factor_random)
			}
			func(tt_data,indi_tAnal=tindi_tAnal,indind=t_ind,...)
		},tindi_tAnal=indi_tAnal[[t_ind]],nperm=nperm,cov=cov,total.reads=total.reads,type=type,...)
	}
}

linearize_out<-function(list){
  outputlist<-list()
  for ( i in 1:length(outoutout)){
    if(class(outoutout[[i]])!="matrix"){
      for ( j in 1:length(outoutout[[i]])){
        outputlist<-c(outputlist,outoutout[[i]][j])
      }
    }else{
      outputlist<-c(outputlist,outoutout[i])
    }
  }
  return(outputlist)
}