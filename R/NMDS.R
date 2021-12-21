NMDS_gogo<-function(mb_table,group,subgroup,filename="NMDS.pdf",width,height,circlesize,textSize,labSize){
  sol<-metaMDS(mb_table,distance = "bray", k = 2, trymax = 50)
  NMDS=data.frame(x=sol$point[,1],y=sol$point[,2],Group=group,Subgroup=subgroup)
  plot.new()
  ord<-ordiellipse(sol, group ,display = "sites", kind ="sd", conf = 0.95, label = T)
  dev.off()
  veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
  {
    theta <- (0:npoints) * 2 * pi/npoints
    Circle <- cbind(cos(theta), sin(theta))
    t(center + scale * t(Circle %*% chol(cov)))
  }
  
  #Generate ellipse points
  df_ell <- data.frame()
  for(g in levels(NMDS$Group)){
    if(g!="" && (g %in% names(ord))){
      
      df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$Group==g,],
                                                       veganCovEllipse(ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale)))
                                    ,Group=g))
    }
  }
  NMDS.mean=aggregate(NMDS[,1:2],list(group=NMDS$Group),mean)
  
  library(ggplot2)
  
  shape_values<-seq(1,11)
  
  p<-ggplot(data=NMDS,aes(x,y,colour=Group))
  p<-p+ annotate("text",x=NMDS.mean$x,y=NMDS.mean$y,label=NMDS.mean$group,size=textSize)
  p<-p+ geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2), size=circlesize, linetype=2)
  p<-p+geom_point(aes(shape=Subgroup))+scale_shape_manual(values=shape_values)+theme_bw() 
  png(filename,width=width,height=height)
  print(p)
  dev.off()
  print(p)
}