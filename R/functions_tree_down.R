find_childtips<-function(k,tree){
  fils <- NULL
  pere <- res <- k
  repeat {
    for (i in 1:length(pere)) fils <- c(fils, tree$edge[,2][tree$edge[, 1] == pere[i]])
    res <- c(res, fils)
    pere <- fils
    fils <- NULL
    if (length(pere) == 0) 
      break
  }
  return(res[res<=Ntip(tree) ])
}

fix_tree<-function(tree2){
  #browser()
  tmp_tree<-tree2
  emat<-matrix(nrow=dim(tree2$edge)[1],ncol=dim(tree2$edge)[2])
  emat<-tree2$edge
  considered_node<-c()
  for ( i in 1:dim(emat)[1]){
    # considered_node<-c(considered_node,tree2$edge[i,1])
    if(tree2$edge[i,2]>Ntip(tree2) &  tree2$edge[i,2]>max(tree2$edge[1:i,1])){
      # tvalue<-tree2$edge[i,2]
      emat[tree2$edge==tree2$edge[i,2]]<-max(tree2$edge[1:i,1])+1
      emat[tree2$edge!=tree2$edge[i,2]&tree2$edge>Ntip(tree2)&tree2$edge>max(tree2$edge[1:i,1])]<-tree2$edge[tree2$edge!=tree2$edge[i,2]&tree2$edge>Ntip(tree2)&tree2$edge>max(tree2$edge[1:i,1])]+1
    }
  }
  tmp_tree$edge<-emat
  return(tmp_tree)
}
  
call_whattodo<-function(k_up,tree,ancestor=""){
 browser()
  base<-1:Ntip(tree)
  temp_output<-NULL
  temp_output_comp<-list(NULL)
  k<-tree$edge[,2][tree$edge[,1]==k_up][1]
  ancestor<-paste0(ancestor,k,sep="-")
  # cat(paste0("\nfocusing node:",(ancestor)))
  # cat(" target:",find_childtips(k))
  
  if(k<=Ntip(tree)){
    # cat(" left_stop",k)
    tempno<-vector("logical",Ntip(tree))
    temp_output_comp[[1]]$group1<-base%in%k
    k2<-tree$edge[,2][tree$edge[,1]==k_up][2]
    # cat(paste0(" / contrast node:",(k2)))
    if(k2<=Ntip(tree)){
      # cat(" right_stop",k2)
      temp_output_comp[[1]]$group2<-base%in%k2
      temp_output_comp[[1]]$root<-ancestor
      temp_output<-append(temp_output,temp_output_comp)
    }else{
      # cat(" right_group",find_childtips(k2,tree))
      temp_output_comp[[1]]$group2<-base%in%find_childtips(k2,tree)
      temp_output_comp[[1]]$root<-ancestor
      temp_output<-append(temp_output,temp_output_comp)
      temp_output<-append(temp_output,call_whattodo(k2,tree,ancestor))
    }
  }else{
    # cat(" left_group",find_childtips(k,tree))
    temp_output_comp[[1]]$group1<-base%in%find_childtips(k,tree)
    k2<-tree$edge[,2][tree$edge[,1]==k_up][2]
    # cat(paste0(" / contrast node:",(k2)))
    if(k2<=Ntip(tree)){
      # cat(" right_stop",k2)
      temp_output_comp[[1]]$group2<-base%in%k2
      temp_output_comp[[1]]$root<-ancestor
      temp_output<-append(temp_output,temp_output_comp)
    }else{
      # cat(" right_group",find_childtips(k2,tree))
      temp_output_comp[[1]]$group2<-base%in%find_childtips(k2,tree)
      temp_output_comp[[1]]$root<-ancestor
      temp_output<-append(temp_output,temp_output_comp)
      temp_output<-append(temp_output,call_whattodo(k2,tree,ancestor))
      # call_whattodo(k2,tree,ancestor)
    }
    # call_whattodo(k,tree,ancestor)
    temp_output<-append(temp_output,call_whattodo(k,tree,ancestor))
  }
  return(temp_output)
}


call_whattodo_part<-function(k_up,tree){
  child_nodes<-tree$edge[,2][tree$edge[,1]==k_up]
  temp_output<-list()
  g1_vec<-vector("logical",Ntip(tree))
  g2_vec<-vector("logical",Ntip(tree))
  g1_vec[find_childtips(child_nodes[1],tree)]<-TRUE
  g2_vec[find_childtips(child_nodes[2],tree)]<-TRUE
  temp_output$group1<- g1_vec
  temp_output$group2<- g2_vec
  temp_output$root<-k_up
  return(temp_output)
}

call_whattodo_v2<-function(tree){
  lapply(1:Nnode(tree)+Ntip(tree),call_whattodo_part,tree)
}

findroot<-function(to,tree,ancester=to){
  ind_edge<-tree$edge[,2]==to
  find<-tree$edge[ind_edge,1]
  if(to==Ntip(tree)+1){
    return(Ntip(tree)+1)
  }
  for(i in find){
    ancester<-c(ancester,i);
    if(i==Ntip(tree)+1){
      return(c(ancester))
    }else{
      return(findroot(i,tree,ancester=ancester))
    }
  }
}

howlong<-function(from,to,tree,total_length=0,ancester=to){
  # browser()
  ind_edge<-tree$edge[,2]==to
  find<-tree$edge[ind_edge,1]
  total_length<-total_length+tree$edge.length[ind_edge]
  if(from==to){
    return(c(0,to))
  }
  for(i in find){
    ancester<-paste0(ancester,"-",i);
    if(i==from){
      return(c(total_length,ancester))
    }else{
      return(howlong(from,i,tree,total_length=total_length,ancester=ancester))
    }
  }
}

getdistance<-function(from,to,tree){
  if(from==to){
    return(list(0,0,0))
  }else{
    candi<-intersect(findroot(from,tree),findroot(to,tree))
    candi2<-candi[candi>Ntip(tree)]
    maxim<-max(candi2)
    list(howlong(maxim,to,tree),howlong(maxim,from,tree),as.numeric(howlong(maxim,to,tree)[1])+as.numeric(howlong(maxim,from,tree)[1]))
  }
}

# ready_nodes<-test_nodes

nodes_to_distnc<-function(ready_nodes, tree){
  dstnc_btw_testnodes<-apply(outer(ready_nodes,ready_nodes,"paste"),2,function(ttdd,tt_tree){
    la_out<-str_split(ttdd," ")
    return(
      sapply(la_out,function(data,t_tree){
      temp_out <-as.numeric(data)
      getdistance(temp_out[1],temp_out[2],t_tree)[[3]]
    },tt_tree)
    )
  },tree)
  return(list(ori=dstnc_btw_testnodes,exp=exp(-dstnc_btw_testnodes)))
}


gettheSubtree<-function(whichone,phylum_list_perphy,phylums_candidates,pruned.tree2,taxo_sorted_df_table){
  #browser()
  ind_perphy_out<-phylum_list_perphy!=phylums_candidates[whichone]
  pruned.tree_perphy<-drop.tip(pruned.tree2,pruned.tree2$tip.label[ind_perphy_out])
  outout_taxonomy_perphy<-list(full=taxo_sorted_df_table[match(pruned.tree_perphy$tip.label,taxo_sorted_df_table[,1]),2],code=taxo_sorted_df_table[match(pruned.tree_perphy$tip.label,taxo_sorted_df_table[,1]),1])
  list(pruned.tree_perphy,outout_taxonomy_perphy)
}


get_target_tips<-function(whatnode,whichone,pruned.tree_perphy){
  #browser()
  x<-whatnode
  dqt <- outer(x, quantile(x), function(x,y) sqrt( (x-y)^2) )
  list_q<-apply(dqt, 2, which.min)
  if(whichone=="q1"){
    c(Ntip(pruned.tree_perphy)+list_q[2],x[list_q[2]])
  }else if(whichone=="q2"){
    c(Ntip(pruned.tree_perphy)+list_q[3],x[list_q[3]])
  }else if(whichone=="q3"){
    c(Ntip(pruned.tree_perphy)+list_q[4],x[list_q[4]])
  }else if(whichone=="random"){
    ind_node<-sample(1:(Ntip(pruned.tree_perphy)-1),1)
    c(Ntip(pruned.tree_perphy)+ind_node,x[ind_node])
  }
}

splitTips<-function(cs){
  if(length(cs)==1){
    cs
  }else{
    tout<-sapply(1:length(cs),function(ind){
      sample(c(T,F),1)
    })
    while((prod(tout)==1 | sum(tout)==0)){
      tout<-sapply(1:length(cs),function(ind){
        sample(c(T,F),1)
      })
      #print(tout)
    }
    tempS_left<-cs[tout]
    tempS_right<-cs[!tout]
    list(splitTips(tempS_left),splitTips(tempS_right))
  }
}
buildEdge<-function(resultsTmp,sel_node,resultMat=matrix(NA,nrow=0,ncol=2)){
  # browser()
  edgeMat<-matrix(NA,nrow=2,ncol=2)
  if(length(resultsTmp)==2){
    
    
    leftbase<-resultsTmp[[1]]
    if(nrow(resultMat)==0){
      startnode<-sel_node
    }else{
      startnode<-max(max(resultMat[,2])  ,sel_node)
    }
    if(length(leftbase)!=1){
      edgeMat[1,]<-c(sel_node,startnode+1)
    }else{
      edgeMat[1,]<-c(sel_node,leftbase)
    }
    resultMat<-rbind(resultMat,edgeMat[1,])
    resultMat<-buildEdge(leftbase,edgeMat[1,2],resultMat)
    
    rightbase<-resultsTmp[[2]]
    if(nrow(resultMat)==0){
      startnode<-sel_node
    }else{
      startnode<-max(max(resultMat[,2])  ,sel_node)
    }
    if(length(rightbase)!=1){
      edgeMat[2,]<-c(sel_node,startnode+1)
    }else{
      edgeMat[2,]<-c(sel_node,rightbase)
    }
    
    resultMat<-rbind(resultMat,edgeMat[2,])
    resultMat<-buildEdge(rightbase,edgeMat[2,2],resultMat)
  }else{
    resultMat
  }
  return(resultMat)
}


getRandomTree<-function(tmpTree,factor=0.5){
# browser()
  Nt<-Ntip(tmpTree)
  Nn<-Nnode(tmpTree)
  node_seq<-1:Nnode(tmpTree)+Ntip(tmpTree)
  list_c_tips<-lapply(node_seq,find_childtips,tmpTree)
  no_childs<-sapply(list_c_tips,length)
  
  number_to_minimize<-no_childs-Nt*factor

  ind_pos<-which(number_to_minimize>=0)
  candidate<-which(number_to_minimize==min(number_to_minimize[number_to_minimize>=0]))
  if(length(candidate)==1){
    ind_target<-candidate
  }else{
    ind_target<-sample(candidate,1)
  }
  # ind_target
  
  # no_childs[ind_target]
  sel_node<-node_seq[ind_target]
  # sel_node
  # plot(tmpTree)

  
  cs<-find_childtips(sel_node,tmpTree)
  results<-splitTips(cs)
  re<-buildEdge(results,sel_node)
  # tmpTree$edge
  # tmpTree
  start<-min(which(tmpTree$edge[,1]==sel_node))
  end<-start+dim(re)[1]-1
  #all true or all false
  tmpTree$edge[start:end,]<-re
  tmpTree$edge.length[start:end]<-sample(tmpTree$edge.length[start:end])
  return(tmpTree)
}
