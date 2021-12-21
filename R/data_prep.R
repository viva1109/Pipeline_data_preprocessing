make_pheno_data_antipro<-function(file,t_df_table,antipro=T,WhatToOne,readcountCut=10000){
 #browser()
  phenotable<-read.csv(file,stringsAsFactors = F)
  phenotable<-phenotable[nchar(phenotable[,1])!=0,]
  phenotable_matched<-phenotable[match(names(t_df_table)[-(1:3)],phenotable[,1]),]
  if (antipro){
	ind_nobiotic<-phenotable_matched$Antibiotic=="-"&phenotable_matched$Probiotic=="-"
  }else{
	ind_nobiotic<-rep(T, dim(phenotable)[1])
  } 
  phenos_txt<-phenotable_matched[ind_nobiotic,2]
  phenos<-ifelse(str_detect(phenos_txt,WhatToOne),1,0)
  df_table_formed_bf<-matrix(as.numeric(t(t_df_table[,-(1:3)])[ind_nobiotic,]),ncol=dim(t_df_table)[1])
  colnames(df_table_formed_bf)<-t_df_table[,1]
  rownames(df_table_formed_bf)<-names(t_df_table)[-(1:3)][ind_nobiotic]
  summed_df<-apply(df_table_formed_bf,1,sum)
  if(!is.null(readcountCut)){
	ind_out<-summed_df>=readcountCut
  }else{
	ind_out<-1:length(summed_df)
  }
  
  df_table_formed<-df_table_formed_bf[ind_out,]
  return(list(phenos=phenos[ind_out],data=df_table_formed))
}


gettable<-function(nooption_ori,sample_labels){
#browser()
  #   dim(nooption_ori)
  #   head(nooption_ori)
  mg2<-nooption_ori[,c(dim(nooption_ori)[2],1:(dim(nooption_ori)[2]-1))]
  #   head(mg2)
  #   dim(mg2)
  #   dim(nooption_ori)
  mg_title<-str_split(mg2$taxonomy,"; ")
  mg_title2<-sapply(mg_title,function(data){
    out<-""
    for (i in 6:1){
      if(!is.na(data[i])){
        if(str_detect(data[i],"D")){
          out<-paste0(data[i:7],collapse="_")
          break;
        }
      }
    }
    return(out)
  })
  
  mg2$g_s<-mg_title2
  mg3<-mg2[,c(2,1,dim(mg2)[2],3:(dim(mg2)[2]-1))]
  # head(mg3)
  a<-names(mg3)
  if(!is.null(sample_labels)){
  	sample_labels[,1]<-gsub("_","",sample_labels[,1])
	# head(sample_labels)
  
  out<-sample_labels[,2][match(a,sample_labels$ID)]
 
  }else{
  out<-a
  }
  
  out[1:3]<-c("taxon_code","taxonomy","genus_species")
  names(mg3)<-out
  return(mg3)
  # write.csv(mg3,out_where,row.names = F)
}

howmany_sps<-function(tmp_whichone,candidates,phylum_list_perphy){
  ind_perphy_out<-phylum_list_perphy!=candidates[tmp_whichone]
  sum(!ind_perphy_out)
}


setAnalysis<-function(n_sim,phenos_original,plus_minus=F,allOTUcausal=F){

t<-NA
for ( i in 1:length(propor_factor_list)){
  for (z in 1:length(beta_per_list)){
    for(z2 in 1:length(ww_list)){
      tmp_info<-unlist(df_genus_info[ww_list[z2],])
      
      t<-rbind(t,c(propor_factor_list[i],beta_per_list[z],ww_list[z2],node_quantile,plus_minus,tmp_info))
      #t<-rbind(t,cbind(propor_factor_list[i],beta_per_list[z],ww_list[z2],node_quantile,plus_minus,df_genus_info[ww_list[z2],]))
    }
  }
}
#browser()
t_table<-t[-1,]
t_table2<-data.frame(cbind(1:dim(t_table)[1],t_table),stringsAsFactors=F)
t_table2[,1:4]<-apply(t_table2[,1:4],2,as.numeric)
colnames(t_table2)<-c("order","propor_factor_list","beta_per_list","ww_list","node_quantile","plus_minus",colnames(df_genus_info))
  go_phylum_setting<-vector("list",dim(t_table2)[1])
  pernu_type1<-lapply(1:n_sim,function(data){shuffle(1:length(phenos_original))})
  t_table2$num_under_node<-rep(0, dim(t_table2)[1])
  for (kk in 1:dim(t_table2)[1]){
    ww<-t_table2$ww_list[kk] 
    propor_factor<-t_table2$propor_factor_list[kk]
    subtree_tmp<-gettheSubtree(ww,upperlevel_list_perphy,upperlevel_candidates,pruned.tree2,taxo_sorted_df_table=taxo_sorted_df_table)
    if(!is.rooted(subtree_tmp[[1]])){
       subtree_tmp[[1]]<-root(subtree_tmp[[1]], 1, r = TRUE)
	tmp_tee<-subtree_tmp[[1]]
	if(length(which(!(tmp_tee$edge[,1]<tmp_tee$edge[,2] | tmp_tee$edge[,2]<=Ntip(tmp_tee))))>0){
		subtree_tmp[[1]]<-fix_tree(tmp_tee)
	}
    }
    pruned.tree_perphy<-subtree_tmp[[1]]
    ind_gogo<-match(subtree_tmp[[2]]$code,colnames(filterd_table))
    whatnode<-sapply((Ntip(pruned.tree_perphy)+1):(Ntip(pruned.tree_perphy)+Nnode(pruned.tree_perphy)),function(data,pt){length(find_childtips(data,pt))},pruned.tree_perphy)
    # target_tips<-Ntip(pruned.tree_perphy)+which(whatnode>=(Ntip(pruned.tree_perphy)/4)&whatnode<=((Ntip(pruned.tree_perphy)/4)*3))
    
    t_table2$num_under_node[kk]<-t_table2$node_quantile[kk]
    t_table2$no_causal_mic[kk]<-"random"
    if(propor_factor==0){
	t_table2$no_causal_mic[kk]<-1
    }
    
    go_phylum_setting[[kk]]$subtree_tmp<-subtree_tmp
    go_phylum_setting[[kk]]$permu_indcausal_pheno<-lapply(1:n_sim,function(data,t_candi_factor,t_hmn_list){
      if(allOTUcausal){
	      target_tips_info<- Ntip(pruned.tree_perphy)+which.max(whatnode)
      }else{
	      target_tips_info<-get_target_tips(whatnode,t_table2$node_quantile[kk],pruned.tree_perphy)
      }
      ind_which<-target_tips_info[1]
	#print(ind_which)
	#print(pruned.tree_perphy)
      candi_factor<-find_childtips(ind_which,pruned.tree_perphy)
      len_candi_factor<-length(candi_factor)
      no_causal<-len_candi_factor*propor_factor
      if(propor_factor==0){
        hmn_list<-rep(1,n_sim)
      }else{
        no_causal<-max(no_causal,1)
        floor_no_causal<-floor(no_causal)
        hmn_list<-rep(floor_no_causal+1,n_sim)
        hmn_list[1:(n_sim*(1-(no_causal-floor_no_causal)))]<-floor_no_causal
      }
	#print("lala1")
	#print(candi_factor)
      ind_cs<-sample(candi_factor,sample(hmn_list,1))
	#print("lala2")
      #print(ind_cs)
      list(ind_causal=ind_cs,pheno_index=pernu_type1[[data]])
    })
    go_phylum_setting[[kk]]$t_table<-t_table2[kk,]
    go_phylum_setting[[kk]]$plus_minus<-plus_minus
    go_phylum_setting[[kk]]$ind_gogo<-ind_gogo
    # go_phylum_setting[[kk]]$causal_mic_ind<-candi_factor[ind_sampled_factor]
    if(kk%%50==0){
      cat(kk,"\n")
    }
    
  }
  return(list(go_phylum_setting=go_phylum_setting,t_table=t_table2))
}


setAnalysisALL_C<-function(n_sim,phenos_original,plus_minus=F){

t<-NA
for ( i in 1:length(propor_factor_list)){
  for (z in 1:length(beta_per_list)){
    for(z2 in 1:length(ww_list)){
      tmp_info<-unlist(df_genus_info[ww_list[z2],])
      
      t<-rbind(t,c(propor_factor_list[i],beta_per_list[z],ww_list[z2],node_quantile,plus_minus,tmp_info))
      #t<-rbind(t,cbind(propor_factor_list[i],beta_per_list[z],ww_list[z2],node_quantile,plus_minus,df_genus_info[ww_list[z2],]))
    }
  }
}
#browser()
t_table<-t[-1,]
t_table2<-data.frame(cbind(1:dim(t_table)[1],t_table),stringsAsFactors=F)
t_table2[,1:4]<-apply(t_table2[,1:4],2,as.numeric)
colnames(t_table2)<-c("order","propor_factor_list","beta_per_list","ww_list","node_quantile","plus_minus",colnames(df_genus_info))
  go_phylum_setting<-vector("list",dim(t_table2)[1])
  pernu_type1<-lapply(1:n_sim,function(data){shuffle(1:length(phenos_original))})
  t_table2$num_under_node<-rep(0, dim(t_table2)[1])
  for (kk in 1:dim(t_table2)[1]){
    ww<-t_table2$ww_list[kk] 
    propor_factor<-t_table2$propor_factor_list[kk]
    subtree_tmp<-gettheSubtree(ww,upperlevel_list_perphy,upperlevel_candidates,pruned.tree2,taxo_sorted_df_table=taxo_sorted_df_table)
    if(!is.rooted(subtree_tmp[[1]])){
       subtree_tmp[[1]]<-root(subtree_tmp[[1]], 1, r = TRUE)
	tmp_tee<-subtree_tmp[[1]]
	if(length(which(!(tmp_tee$edge[,1]<tmp_tee$edge[,2] | tmp_tee$edge[,2]<=Ntip(tmp_tee))))>0){
		subtree_tmp[[1]]<-fix_tree(tmp_tee)
	}
    }
    pruned.tree_perphy<-subtree_tmp[[1]]
    ind_gogo<-match(subtree_tmp[[2]]$code,colnames(filterd_table))
    whatnode<-sapply((Ntip(pruned.tree_perphy)+1):(Ntip(pruned.tree_perphy)+Nnode(pruned.tree_perphy)),function(data,pt){length(find_childtips(data,pt))},pruned.tree_perphy)
    # target_tips<-Ntip(pruned.tree_perphy)+which(whatnode>=(Ntip(pruned.tree_perphy)/4)&whatnode<=((Ntip(pruned.tree_perphy)/4)*3))
    
    t_table2$num_under_node[kk]<-t_table2$node_quantile[kk]
    t_table2$no_causal_mic[kk]<-"random"
    if(propor_factor==0){
	t_table2$no_causal_mic[kk]<-1
    }
    
    go_phylum_setting[[kk]]$subtree_tmp<-subtree_tmp
    go_phylum_setting[[kk]]$permu_indcausal_pheno<-lapply(1:n_sim,function(data,t_candi_factor,t_hmn_list){
      #target_tips_info<-get_target_tips(whatnode,t_table2$node_quantile[kk],pruned.tree_perphy)
      target_tips_info<- Ntip(pruned.tree_perphy)+which.max(whatnode)
      print(target_tips_info)
      ind_which<-target_tips_info[1]
	#print(ind_which)
	#print(pruned.tree_perphy)
      candi_factor<-find_childtips(ind_which,pruned.tree_perphy)
      print(candi_factor)
      len_candi_factor<-length(candi_factor)
      no_causal<-len_candi_factor*propor_factor
      if(propor_factor==0){
        hmn_list<-rep(1,n_sim)
      }else{
        no_causal<-max(no_causal,1)
        floor_no_causal<-floor(no_causal)
        hmn_list<-rep(floor_no_causal+1,n_sim)
        hmn_list[1:(n_sim*(1-(no_causal-floor_no_causal)))]<-floor_no_causal
      }
	#print("lala1")
	#print(candi_factor)
      ind_cs<-sample(candi_factor,sample(hmn_list,1))
	#print("lala2")
      #print(ind_cs)
      list(ind_causal=ind_cs,pheno_index=pernu_type1[[data]])
    })
    go_phylum_setting[[kk]]$t_table<-t_table2[kk,]
    go_phylum_setting[[kk]]$plus_minus<-plus_minus
    go_phylum_setting[[kk]]$ind_gogo<-ind_gogo
    # go_phylum_setting[[kk]]$causal_mic_ind<-candi_factor[ind_sampled_factor]
    if(kk%%50==0){
      cat(kk,"\n")
    }
    
  }
  return(list(go_phylum_setting=go_phylum_setting,t_table=t_table2))
}




setAnalysisH0<-function(n_sim,phenos_original,plus_minus=F){
propor_factor_list<-0
beta_per_list<-0
t<-NA
for ( i in 1:length(propor_factor_list)){
  for (z in 1:length(beta_per_list)){
    for(z2 in 1:length(ww_list)){
      tmp_info<-unlist(df_genus_info[ww_list[z2],])
      
      t<-rbind(t,c(propor_factor_list[i],beta_per_list[z],ww_list[z2],node_quantile,plus_minus,tmp_info))
    }
  }
}
#browser()
t_table<-t[-1,]
t_table2<-data.frame(cbind(1:dim(t_table)[1],t_table),stringsAsFactors=F)
t_table2[,1:4]<-apply(t_table2[,1:4],2,as.numeric)
colnames(t_table2)<-c("order","propor_factor_list","beta_per_list","ww_list","node_quantile","plus_minus",colnames(df_genus_info))
  go_phylum_setting<-vector("list",dim(t_table2)[1])
  pernu_type1<-lapply(1:n_sim,function(data){shuffle(1:length(phenos_original))})
  for (kk in 1:dim(t_table2)[1]){
    ww<-t_table2$ww_list[kk] 
    propor_factor<-t_table2$propor_factor_list[kk]
    #subtree_tmp<-gettheSubtree(ww,upperlevel_list_perphy,upperlevel_candidates,pruned.tree2,taxo_sorted_df_table=taxo_sorted_df_table)
    subtree_tmp<-gettheSubtree(ww,upperlevel_list_perphy,upperlevel_candidates,pruned.tree2,taxo_sorted_df_table=taxo_sorted_df_table)
    if(!is.rooted(subtree_tmp[[1]])){
       subtree_tmp[[1]]<-root(subtree_tmp[[1]], 1, r = TRUE)
    }
    pruned.tree_perphy<-subtree_tmp[[1]]
    ind_gogo<-match(subtree_tmp[[2]]$code,colnames(filterd_table))
    whatnode<-sapply((Ntip(pruned.tree_perphy)+1):(Ntip(pruned.tree_perphy)+Nnode(pruned.tree_perphy)),function(data,pt){length(find_childtips(data,pt))},pruned.tree_perphy)
    # target_tips<-Ntip(pruned.tree_perphy)+which(whatnode>=(Ntip(pruned.tree_perphy)/4)&whatnode<=((Ntip(pruned.tree_perphy)/4)*3))
    target_tips_info<-get_target_tips(whatnode,t_table2$node_quantile[kk],pruned.tree_perphy)
    t_table2$num_under_node[kk]<-target_tips_info[2]
    ind_which<-target_tips_info[1]
    candi_factor<-find_childtips(ind_which,pruned.tree_perphy)
    len_candi_factor<-length(candi_factor)
    
    no_causal<-len_candi_factor*propor_factor
    if(propor_factor==0){
      # ind_sampled_factor<-sample(1:len_candi_factor,1)
      t_table2$no_causal_mic[kk]<-1
      floor_no_causal<-floor(no_causal)
      hmn_list<-rep(1,n_sim)
    }else{
      # ind_sampled_factor<-sample(1:len_candi_factor,ceiling(len_candi_factor*propor_factor))
      no_causal<-max(no_causal,1)
      t_table2$no_causal_mic[kk]<-no_causal
      floor_no_causal<-floor(no_causal)
      hmn_list<-rep(floor_no_causal+1,n_sim)
      hmn_list[1:(n_sim*(1-(no_causal-floor_no_causal)))]<-floor_no_causal
    }
    go_phylum_setting[[kk]]$subtree_tmp<-subtree_tmp
    go_phylum_setting[[kk]]$permu_indcausal_pheno<-lapply(1:n_sim,function(data,t_len_candi_factor,t_hmn_list){
      list(ind_causal=sample(1:t_len_candi_factor,t_hmn_list[data]),pheno_index=pernu_type1[[data]])
    },len_candi_factor,hmn_list)
    go_phylum_setting[[kk]]$t_table<-t_table2[kk,]
    go_phylum_setting[[kk]]$plus_minus<-plus_minus
    go_phylum_setting[[kk]]$ind_gogo<-ind_gogo
    # go_phylum_setting[[kk]]$causal_mic_ind<-candi_factor[ind_sampled_factor]
    if(kk%%50==0){
      cat(kk,"\n")
    }
    
  }
  return(list(go_phylum_setting=go_phylum_setting,t_table=t_table2))
}