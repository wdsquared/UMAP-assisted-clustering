# COVID-19 UMAP-assisted clustering analysis pipeline - Model development and validation 

# David Greenwood - Univeristy of Birmingham 2022 - drg707@student.bham.ac.uk, wdsquared.projects@gmail.com

# Including modified code originally written by Thomas Taverner https://research.birmingham.ac.uk/en/persons/thomas-taverner

#----------------------------------------------------------------#

# MODEL DEVELOPMENT 

### Libraries 
{
  require(mice) # mice(), complete()
  require(FactoMineR) # FAMD()
  require(factoextra) # get_eigenvalue(), fviz_screeplot()
  require(intrinsicDimension) # maxLikGlobalDimEst()
  require(uwot) # UMAP()
  require(mclust) # mclust()
  require(cluster) # silhouette()
  require(forcats) # fct_reorder()
  require(DataCombine) # PercChange()
  require(confintr) # CI_proportion()
  require(sjstats) # bootstrap()
  
  require(dplyr) # piping %>%
  require(reshape2) # dcast()
  require(ggplot2) # plotting 
  
}

### Functions 
{
  
  # Mclust wrapper with silhoutte and BIC comparison
  GMM_comparison <- function(data_input = NULL, # Data for model fitting (e.g. UMAP embedding)
                             K_vec = NULL,
                             modelNames = "VVV",
                             initialization = list("hcpairs"),
                             seed_no = NULL, 
                             saved_seed = NULL,...){
    
    # Lists for models 
    # mod_list <- list()
    
    # Data frame for BIC
    BIC_df<- data.frame(matrix(ncol=2,nrow=length(K_vec)))
    colnames(BIC_df)<-c("K","BIC")
    
    # List for silhouette analysis 
    sil_list <- list()
    sil_df<- data.frame(matrix(ncol=2,nrow=length(K_vec)))
    colnames(sil_df)<-c("K","av.sil")
    
    # Run GMM looped 
    for(j in 1:length(K_vec)){
      
      print(paste0("K ",K_vec[j]))
      print(paste0("modelNames ",modelNames))
      print(paste0("initialization ", initialization))
      
      GMM.mod <- NULL
      
      # Set seed
      if(!is.null(saved_seed)){
        .Random.seed <- saved_seed
        print(paste0("Setting saved seed"))
      } else{
        if(!is.null(seed_no)){
          print(paste0("Setting seed number: ",seed_no))
          set.seed(seed_no)
        } else{
          print(paste0("Warning: no seed set"))
        }
      }
      
      # RUN GMM
      GMM.mod <- Mclust(data = data_input,
                        G = K_vec[j],
                        modelNames = modelNames,
                        initialization = initialization)
      
      # Store model
      #mod_list[[paste0("K_",K_vec[j])]] <-GMM.mod
      
      # Run silhouette  
      if(is.null(GMM.mod)){
        
        # Check if model failed to fit 
        print(paste0("Error with ",paste0("K_",K_vec[j])))
        next
        
      } else {
        
        # Sil analysis 
        data_dist<- dist(data_input) # Distance input
        si <- silhouette(GMM.mod$classification, data_dist)
        sil_df[j,"av.sil"] <- mean(si[,"sil_width"])
        sil_df[j,"K"] <- paste0("K_",K_vec[j])
        
        # BIC comparison
        BIC_tmp<-GMM.mod$bic
        BIC_df[j,"K"] <- paste0("K_",K_vec[j])
        BIC_df[j,"BIC"]<- BIC_tmp
        
      }
      
    }
  
    return(list("av_sil"=sil_df,
                "BIC"=BIC_df))
    
  }
  
  # Visualise BIC - modified factoextra::fviz_mclust
  fviz_mclust_bic_mod<-function (object, model.names = NULL, shape = 19, color = "model", 
                                 palette = NULL, legend = NULL, main = "Model selection", 
                                 xlab = "Number of components", ylab = "BIC", 
                                 ...) {
    if (!inherits(object, "Mclust")) 
      stop("An object of class Mclust is required.")
    best_model <- object$modelName
    number_of_cluster <- object$G
    x <- object$BIC
    n <- ncol(x)
    dnx <- dimnames(x)
    x <- matrix(as.vector(x), ncol = n)
    dimnames(x) <- dnx
    x <- as.data.frame(x, stringsAsFactors = TRUE)
    if (is.null(model.names)) 
      model.names <- dimnames(x)[[2]]
    x <- x[, model.names, drop = FALSE]
    x <- cbind.data.frame(cluster = rownames(x), x, stringsAsFactors = TRUE)
    x <- tidyr::gather_(x, key_col = "model", value_col = "BIC", 
                        gather_cols = colnames(x)[-1])
    x <- x[!is.na(x$BIC), , drop = FALSE]
    x$model <- factor(x$model, levels = dnx[[2]])
    x$cluster <- factor(x$cluster, levels = unique(x$cluster))
    if (ggpubr:::.is_col_palette(palette)) 
      palette <- ggpubr:::.get_pal(palette, k = length(model.names))
    ggline.opts <- list(data = x, x = "cluster", y = "BIC", 
                        group = "model", color = color, shape = shape, 
                        palette = palette, main = main, xlab = xlab, ylab = ylab, 
                        ...)
    p <- do.call(ggpubr::ggline, ggline.opts) #+ labs(subtitle = paste0("Best model: ", best_model, " | Optimal clusters: n = ", number_of_cluster)) + 
    #geom_vline(xintercept = number_of_cluster, linetype = 2, 
    #           color = "red") + theme(legend.title = element_blank())
    if (missing(legend)) 
      p + theme(legend.position = c(0.7, 0.2), legend.direction = "horizontal", 
                legend.key.height = unit(0.5, "line")) + guides(color = guide_legend(nrow = 5, 
                                                                                     byrow = TRUE))
    else p + theme(legend.position = legend)
  }
  
  # Compare binary outcome by cluster
  clust_N_event_summary = function(mod=NULL,
                                   outcome=NULL){
    
    stopifnot(all(outcome==0|outcome==1))
    
    mod_name = paste0("K_",mod$G)
    npts<- data.frame(table(mod$classification))
    npts$Var1 <- as.numeric(as.character(npts$Var1))
    
    # Data with outcome by classification 
    tmp<-data.frame("outcome"=outcome,
                    "cluster"=mod$classification)
    
    # Get N patients by cluster in each outcome group
    tmp2<-tmp %>%
      group_by(cluster,outcome) %>%
      summarise(n=n())
    
    # Get N patient by cluster in total 
    tmp3<-tmp %>%
      group_by(cluster) %>%
      summarise(n_total=n())
    
    # Join 
    tmp4<-left_join(tmp2,tmp3)
    
    # Calculate proprtion
    tmp4$prop = tmp4$n/tmp4$n_total
    
    # Calculate range of proportions
    tmp<-data.frame(tmp4,stringsAsFactors = F)
    
    # Filter to positive events (e.g. death )
    tmp<-tmp[tmp$outcome==1,] 
    
    # Median and IQR 
    median_tmp <- summary(tmp$prop)["Median"]
    lower_tmp <- summary(tmp$prop)["1st Qu."]
    upper_tmp <-summary(tmp$prop)["3rd Qu."]
    
    # Min / Max
    tmp2<-cbind(tmp[order(tmp$prop)[1],],group="min")
    tmp3<-cbind(tmp[rev(order(tmp$prop))[1],],group="max")
    tmp<-rbind(tmp2,tmp3)
    min_prop = tmp[tmp$group=="min","prop"]
    min_n =     tmp[tmp$group=="min","n_total"]
    max_prop =  tmp[tmp$group=="max","prop"]
    max_n =     tmp[tmp$group=="max","n_total"]
    diff_prop = max_prop-min_prop
    
    # Summary 
    summary_df<-data.frame("model"=mod_name,
                           "min_prop"=min_prop,
                           "min_n"=min_n,
                           "max_prop"=max_prop,
                           "max_n"=max_n,
                           "diff_prop"=diff_prop,
                           "median_prop"=as.numeric(as.character(paste0(median_tmp))),
                           "lower_quartile"=as.numeric(as.character(paste0(lower_tmp))),
                           "upper_quartile"=as.numeric(as.character(paste0(upper_tmp))),
                           stringsAsFactors = F)
    
    
    out <- list("model"=mod_name,
                "npts"=npts,
                "outcome"=data.frame(tmp4,stringsAsFactors = F),
                "min_max"=tmp,
                "summary"=summary_df)
    
    return(out)
  }
  
  # Silhouette wrapper function for complete data
  clust_sil<-function(data=NULL,
                      labels=NULL){
    
    si <- silhouette(labels, dist(data))
    si<-as.data.frame(cbind(cluster=si[,"cluster"],neighbor=si[,"neighbor"],sil_width=si[,"sil_width"]))  
    
    si2 <- 
      si %>%
      group_by(cluster) %>%
      summarise("N"=n())
    
    for(i in 1:length(unique(si2$cluster))){
      clust <- unique(si2$cluster)[i]
      tmp <- si[si$cluster==clust,]
      if(length(tmp$cluster)>=10){
        sd.sil <- se_fun(tmp$sil_width)
        tmp <- ci_mean(tmp$sil_width)
        si2[si2$cluster==clust,"Q_hat"]=tmp$estimate
        si2[si2$cluster==clust,"S_hat"]=sd.sil
        si2[si2$cluster==clust,"lower"]=tmp$interval[1]
        si2[si2$cluster==clust,"upper"]=tmp$interval[2]
      } else{
        sd.sil <- se_fun(tmp$sil_width)
        tmp <- mean(tmp$sil_width)
        si2[si2$cluster==clust,"Q_hat"]=tmp
        si2[si2$cluster==clust,"S_hat"]=sd.sil
      }
      
    }
    return(list("sil_summary"=si2,
                "sil_by_observation"=si))
  }
  
  # Silhouette wrapper function for multiple imputations
  mi_clust_sil<-function(data_list = NULL,labels = NULL){
    
    si3<-NULL
    
    for(i in 1:length(data_list)){
      si<-clust_sil(data=data_list[[i]],
                    labels=labels)
      si2<-si$sil_summary
      si2$imp<-paste0(i)
      if(is.null(si3)){
        si3<-si2
      } else {
        si3<-rbind(si3,si2)
      }  
    }
    
    si4<-data.frame()
    for(i in 1:length(table(si3$cluster))){
      cluster <- names(table(si3$cluster))[i]
      data_all <- si3[si3$cluster==cluster,]
      si4[i,"cluster"] <- cluster
      stopifnot(all(si3[si3$cluster==cluster,"N"]))
      si4[i,"N"] <- unique(si3[si3$cluster==cluster,"N"])
      rub <- rubinEstimate(estimate = data_all$Q_hat,se = data_all$S_hat)
      si4[i,"pool_Q_hat"] <-rub["rubinMean"]
      si4[i,"pool_lower"] <-rub["lower"]
      si4[i,"pool_upper"] <-rub["upper"]
    }
    si4
  }  

  # Utility functions
  se_fun <- function(x) sd(x)/sqrt(length(x))
  rubinEstimate <- function(estimate, se, K = length(estimate)){
    rubinMean <- mean(estimate)
    if(length(se) == 1){
      rubinVar <- se^2
    } else {
      rubinVar <- var(estimate)*(1 + 1/K) + mean(se^2)
    }
    lower <- rubinMean - 1.96*sqrt(rubinVar)
    upper <- rubinMean + 1.96*sqrt(rubinVar)
    c("rubinMean"=rubinMean,"lower"=lower,"upper"=upper)
  }
  getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  cont_tbl <- function(x,y, make_binary=T,labels=c("censored","death")){
    if(make_binary==T){
      levels(x)[which(levels(x)==labels[2])]<-1
      levels(x)[which(levels(x)==labels[1])]<-0
      levels(y)[which(levels(y)==labels[2])]<-1
      levels(y)[which(levels(y)==labels[1])]<-0
    }
    x<- as.numeric(as.character(x))
    y <- as.numeric(as.character(y))
    tbl <- rbind(c(sum(x), sum(abs(x-1))),
                 c(sum(y), sum(abs(y-1))))
    tbl
  } 
  
}   

### Development data
{
  
  setwd("/path_to/project_folder/")
  
  # Read in data 
  {
    # Incomplete data, wide format, with additional variables such as outcome 
    dev.dat<- read.csv("/path_to/project_folder/csv/development_data.csv")[,-1]
    
    # Imputed data 
    miDat.dev<- read.csv("/path_to/project_folder/Imputed/development_miData.RDS")
    
  }
  
  # Determine variables for UMAP dimenion reduction
  {
    # Variables for dimension reduction by UMAP
    umap_pred <- c("age","charlsondementiaicd10","cancer","CVD","systolicbp",
                   "diastolicbp","heartrate","temperature","resps","o2sats",
                   "Alt","CRP","Hion","HCO3","ur2","egfr","Hb","lymphocytes2",
                   "baseexc","lactate2","cough","fever","delirium",
                   "sex","bmi","neutlymph_ratio","pocfio2","frailty_cts")
    
    # Note categorical & continuous
    cat_pred <- c("charlsondementiaicd10","cancer","CVD","cough","fever","delirium","sex")
    cts_pred <- umap_pred[umap_pred%in%cat_pred==F]
    length(umap_pred)
    
    # Completed data, single imputation
    impt<-complete(miDat.dev,1)[,c(umap_pred)] 
    
    # Check either continuous or binary
    head(impt)
    str(impt)
    for(i in 1:length(umap_pred)){
      vec <- impt[,umap_pred[i]]
      name <- umap_pred[i]
      if(name %in%cat_pred){
        if(all(vec==1|vec==0)==F){
          print(paste0("issue with cat var: ",name))
        } else{
          print(paste0("cat var: ",name))
          #print(summary(vec))
        }
      }
      if(name %in%cts_pred){
        vec = as.numeric(vec)
        if(anyNA(vec)==T){
          print(paste0("issue with cts var: ",name))
        }else{
          print(paste0("cts var: ",name))
          #print(summary(vec))
        }
      }
      summary(impt)
    }
    
    # Make numeric
    impt <- data.frame((lapply(impt,function(x){as.numeric(as.character(x))})),stringsAsFactors = F)
    str(impt)
    
  }
}

### Selecting model parameters 
{
  # Estimate data complexity 
  {
    # FMDA and intrinsic dimensionality estimates used to select a range for D, the number of UMAP dimensions to embed the data
    
    # FAMD 
    {
      
      # Data 
      impt.tmp <- impt
      
      # Factorise categorical
      impt.tmp[,colnames(impt.tmp)%in%cat_pred]<-
        data.frame(lapply(impt.tmp[,colnames(impt.tmp)%in%cat_pred],function(x){
          as.factor(x)
        }))
      str(impt.tmp)
      
      # Apply FAMD
      res.famd <- FAMD (impt.tmp, ncp = ncol(impt.tmp), ind.sup = NULL, graph = TRUE)
      eig.val <- get_eigenvalue(res.famd)
      head(eig.val)  
      
      # Plot
      pl1<- fviz_screeplot(res.famd)+
        theme_bw()+labs(title="",x="FAMD dimensions",y="Explained variance (%)")+
        scale_y_continuous(limits=c(0,15))
      
      # Visualise
      print(pl1)
      
    }
    
    # Maximum likelihood global intrinsic dimension estimation
    {
      
      # Data
      impt.tmp <- impt
      str(impt.tmp)
      impt.tmp<-as.matrix(impt.tmp)
      
      # DF for results
      df<-data.frame(matrix(ncol=3,nrow=0))
      colnames(df)<-c("method","k","estimate")
      
      # Estimate dependant on number of distances used in calculation
      # repeat for range (long run time for larger K)
      tmp<-seq(from=30,to=6098,by = 120)
      
      # Maximum likelhood estimation 
      for(i in 1:length(tmp)){
        out<-maxLikGlobalDimEst(impt.tmp, k=tmp[i],
                                integral.approximation = 'Haro',
                                unbiased = TRUE,
                                neighborhood.aggregation = 'robust')
        df[i,"method"]="maxLikGlobalDimEst"
        df[i,"k"]=tmp[i]
        df[i,"estimate"]=out$dim.est
      }
      plot(df$estimate,df$k)
      
      
      pl2<-ggplot(df,aes(x=k,y=estimate))+
        geom_point()+scale_y_continuous(limits=c(1,10),breaks = seq(from=1,to=10,by=1))+
        geom_line()+labs(x="N. distances considered",y="Dimension estimate")+theme_bw()
      
      print(pl2)
      
      
    }
    
  }
  
  # Choose initial range for D
  print(pl1) # FAMD elbow point around 4-D or 5-D
  print(pl2) # Global intrinsic dimension estiamtes range from 2-D to 8-D, but most estimates < 6-D 
  D_range <-seq(from=2,to=6)
  
  # UMAP dimension reduction - repeated for each D 
  {
    
    # Run UMAP with each D
    umap_list<-list()
    for(j in 1:length(D_range)){
      print(paste0("UMAP dimension reduction from ",paste0(length(umap_pred)),"-D to ",D_range[j],"-D"))
      set.seed(10)
      umap_out <-uwot::umap(impt,
                            scale=T,
                            n_threads = 10,
                            n_neighbors = 40,
                            min_dist = 0.25,
                            learning_rate = 0.5,
                            init="normlaplacian",
                            ret_model = F,
                            n_components = D_range[j],
                            verbose=T)
      
      
      umap_list[[paste0("N_dim_",D_range[j])]]<-umap_out
    }
    remove(umap_out)
    
  }
  
  # GMM clustering - repeated for each K
  {
    
    # Choose K, the number of components in the model
    K_vec = c(2:50)
    
    # Choose constraint for GMM covariance (note, unconstrained covariance matrices will dramatically increase run time)
    mod_type = "VVV" # see ?mclust for more options
    
    # Run comparison 
    {
      # Data frame for average silhouette width 
      sil_out <- data.frame(matrix(ncol=3,nrow=0))
      colnames(sil_out) <- c("N_dim","K","av_sil")
      
      # BIC
      BIC_out <- data.frame(matrix(ncol=3,nrow=0))
      colnames(BIC_out) <- c("N_dim","K","BIC")
      
      # Run
      for(i in 1:length(umap_list)){
        
        # Get embedding
        tmp<-gsub("N_dim_","",names(umap_list)) # Get number of dimensons
        print(paste0(tmp[i])) 
        print(paste0(names(umap_list)[i])) # Name of embedding
        
        # Fit GMM models to this embedding with K components (clusters)
        
        out<-GMM_comparison(data_input = umap_list[[i]],
                            modelNames =   mod_type,
                            BIC=T,
                            seed_no = 10,
                            #saved_seed = seed_save,
                            K_vec=K_vec)
        
        # Note number of dimensions
        out$av_sil$N_dim<-tmp[i]
        
        # Note number of dimensions
        out$BIC$N_dim<-tmp[i]
        
        sil_out<-rbind(sil_out,out$av_sil)
        BIC_out<-rbind(BIC_out,out$BIC)
        
      }
      
      # Combine
      tmp<-dplyr::left_join(sil_out,BIC_out)
      tmp<-tmp[,c("N_dim","K","av.sil","BIC")]
      
      # Format 
      tmp$K_numeric<-as.numeric(gsub("K_","",tmp$K))
      tmp$N_dim<-as.numeric(tmp$N_dim)
      
      # Add NA for models that didnt converge 
      tmp2 <- data.frame("K_numeric"=as.numeric(rep(seq(from=2,to=max(K_vec)),
                                                    max(D_range)-1)),
                         "N_dim"=as.numeric(rep(seq(from=2,to=max(D_range)),
                                                each=max(K_vec)-1)))
      tmp<-dplyr::left_join(tmp2,tmp)
      
      sil_bic_out<-tmp 
      
    }
    
    # Plot comparison
    {
      tmp <- sil_bic_out
      #tmp <- tmp[complete.cases(tmp),]
      
      # Order factor by K 
      tmp$K_ordered <- forcats::fct_reorder(tmp$K, tmp$K_numeric)
      
      # Tidy name 
      levels(tmp$K_ordered) <- gsub("K_","K = ",levels(tmp$K_ordered))
      
      # N dim 
      tmp$N_dim <- factor(tmp$N_dim,levels=seq(from=2,to=max(D_range),by=1))
      
      # order by K, then number of dim 
      tmp<-tmp[order(tmp$K_numeric,
                     as.numeric(as.character(tmp$N_dim))),]
      
      # Calculate delta BIC (change between K clusters, by D dimensional embedding)
      tmp<-DataCombine::PercChange(tmp,
                                   Var="BIC",
                                   GroupVar = "N_dim",
                                   NewVar="Delta_BIC",
                                   slideBy = -1,
                                   type = "percent")
      sil_bic_out <- tmp
      
      # Silhoute plot
      pl3<-ggplot(sil_bic_out,aes(x=K_numeric,
                                  y=av.sil,
                                  group=N_dim,
                                  color=N_dim))+
        geom_line()+
        theme_bw()+
        scale_x_continuous(breaks=seq(from=0,to=max(sil_bic_out$K_numeric),by=5),limits=c(2,max(sil_bic_out$K_numeric)))+
        labs(y="Average sil. width",x="K components",colour="D dimensions")
      
      # Add any lines
      #pl3+
      #  geom_vline(linetype=3,xintercept = 3,color="blue")+
      #  geom_vline(linetype=2,xintercept = 9,color="blue")
      
      
      # plot BIC
      pl4<-ggplot(sil_bic_out, aes(x=K_numeric,
                                   y=Delta_BIC,
                                   group=N_dim,
                                   colour=N_dim))+
        geom_line()+
        theme_bw()+
        scale_x_continuous(breaks=seq(from=0,to=max(sil_bic_out$K_numeric),by=5),limits=c(2,max(sil_bic_out$K_numeric)))+
        labs(y="Delta BIC (%)",x="K components",colour="D dimensions")
      
      
    }
    
  }
  
  # Select maximum D 
  print(pl3) # Silhouette plot: highest average sil. width with 2-D but little difference with >3-D
  print(pl4) # Delta BIC plot: delta BIC varaince in increased substantially with >5-D
  D_value <- 4
  
  # GMM clustering - optimal BIC / max silhoutte 
  {
    # Optimal BIC
    {
      # Check if optimal reached within range of K values tested
      tmp<-sil_bic_out[sil_bic_out$N_dim==D_value,]
      opt_acheived <- any(tmp$Delta_BIC>0) 
      plot(tmp$K_numeric,tmp$BIC) # Check optimal, if not set opt_acheived = FALSE
      
      # Extract optimal BIC 
      if(opt_acheived==T){
        optimal_BIC_K<-min(tmp[tmp$Delta_BIC>0,"K_numeric"],na.rm = T)
        print(paste0("Optimal BIC reached at K = ",optimal_BIC_K))
        pl5<-ggplot(subset(sil_bic_out,sil_bic_out$N_dim==D_value), aes(x=K_numeric,
                                                                        y=BIC,
                                                                        group=N_dim))+
          geom_line()+
          theme_bw()+
          scale_x_continuous(breaks=seq(from=0,to=max(sil_bic_out$K_numeric),by=5),limits=c(2,max(sil_bic_out$K_numeric)))+
          labs(y="BIC",x="K components")
        print(pl5)
      } 
      
      # Refit if not reached
      if(opt_acheived==F) {
        
        # Fit GMM to D dimensional embedding and compare BIC with K = K_vec
        K_vec = 2:100
        
        # Base umap
        dat_umap <- umap_list[[names(umap_list)[grepl(D_value,names(umap_list))]]]
        
        # Compare BIC
        set.seed(10)
        #.Random.seed <- seed_save
        bic_tmp <-   mclustBIC(dat_umap,
                               G = K_vec,
                               modelNames = mod_type,
                               initialization = list("hcpairs"))
        set.seed(10)
        #.Random.seed <- seed_save
        tmp_mod <- Mclust(dat_umap, x = bic_tmp)
        optimal_BIC_K<-summary(tmp_mod)$G
        
        # Plot  
        {
          pl5 <- fviz_mclust_bic_mod(tmp_mod,
                                     main = "Model selection",
                                     xlab = "Number of components",
                                     ylab = "BIC")
          
          # Format plot
          trunc(min(tmp_mod$BIC[1:99],na.rm = T))
          #y_min = -71000
          trunc(max(tmp_mod$BIC[1:99],na.rm = T))
          #y_max = -47000
          pl5<-pl5+theme(axis.text.x = element_text(size = 9),
                         axis.text.y = element_text(angle = 90,size = 9,hjust = 0.9),
                         plot.title = element_blank(),
                         plot.subtitle=element_blank(),
                         legend.title = element_blank(),legend.position="none") +
            scale_x_discrete(breaks = seq(from=0,to=max(K_vec),by=10),position = "bottom")+
            #scale_y_continuous(breaks = seq(from=y_min,to=y_max,by=10000))+
            labs(x="K components")
          #pl5+
          #  geom_vline(linetype=3,xintercept = 49,color="red")+
          #  geom_vline(linetype=2,xintercept = 28,color="red")
          print(pl5)
        }
        
        
      }
    }
    
    # Max silhouete 
    {
      tmp<-sil_bic_out[sil_bic_out$N_dim==D_value,]
      max_sil_K <- tmp[tmp$av.sil==max(tmp$av.sil),"K_numeric"]
    }
    
  }
  
  # Select possible K values
  print(pl3) # By silhouette 
  max_sil_K # Maximum acheived at K=3
  manual_sil_K <- 9 # Secondary peak at K=9
  print(pl5) # By BIC
  optimal_BIC_K # optimal at K=49
  manual_BIC_K <- 29 # Plateu visible at K=29
  K_vec <- sort(c(max_sil_K=max_sil_K,
                  manual_sil_K=manual_sil_K,
                  optimal_BIC_K=optimal_BIC_K,
                  manual_BIC_K=manual_BIC_K))
  print(K_vec)
  
  # Compare models by outcome
  # Calculate event rate by cluster and calculate range
  {
    if(length(unique(K_vec))==1){
      print(paste0("Single K = ",unique(K_vec),", nothing to compare"))
    } else{
      print(paste0("Comparing proportion of events by K = ",paste0(K_vec)," clusters"))
    }
    
   
    
    # Base umap
    dat_umap <- umap_list[[names(umap_list)[grepl(D_value,names(umap_list))]]]
    
    # Fit for each selected K 
    mod_list<-list()
    for(i in 1:length(K_vec)){
      name <- paste0("K_",K_vec[i])
      #.Random.seed <- seed_save
      set.seed(10)
      mod_list[[name]] <- Mclust(dat_umap,
                                 G = K_vec[i],
                                 modelNames = mod_type,
                                 initialization = list("hcpairs"))
      
    }
    
    
    # Compare event rates by cluster for each model 
    clust_events <- list()
    for(i in 1:length(mod_list)){
      tmp<-mod_list[[i]]
      tmp2<-clust_N_event_summary(mod=tmp,
                                  outcome = dev.dat$death28)
      clust_events[[tmp2$model]]<-tmp2
    }
    
    # Format and plot 
    {
      # Bind summaries 
      clust_events.df <- clust_events[[1]]$summary
      for(i in 2:length(clust_events)){
        tmp<-clust_events[[i]]
        clust_events.df<-rbind(clust_events.df,tmp$summary)
      }
      
      # sort 
      clust_events.df$model_name <- gsub("K_","",clust_events.df$model)
      clust_events.df$model_name<-
        factor(clust_events.df$model_name,
               levels=sort(as.numeric(levels(factor(clust_events.df$model_name)))),
               ordered=T)
      
      str(clust_events.df)
      
      # X for annotation
      clust_events.df$x_annot_loc<- seq(from=1,to=nrow(clust_events.df))+0.4
      
      # 
      pl6<-ggplot(aes(x=model_name,
                      y=median_prop),
                  dat=clust_events.df) +
        geom_point()+
        geom_errorbar(aes(x=model_name, ymax=max_prop, ymin=min_prop),  na.rm=TRUE, width = 0.2)+
        theme_bw() +
        labs(x="K",y="Proportion of deaths")+
        theme(legend.position = "top")
      
      for(i in 1:length(unique(clust_events.df$model_name))){
        pl6<-pl6+ 
          annotate("text",x=clust_events.df$x_annot_loc[clust_events.df$model_name==unique(clust_events.df$model_name)[i]],
                   y=clust_events.df$max_prop[clust_events.df$model_name==unique(clust_events.df$model_name)[i]],
                   label=clust_events.df$max_n[clust_events.df$model_name==unique(clust_events.df$model_name)[i]],
                   size=2.7)+
          annotate("text",x=clust_events.df$x_annot_loc[clust_events.df$model_name==unique(clust_events.df$model_name)[i]],
                   y=clust_events.df$min_prop[clust_events.df$model_name==unique(clust_events.df$model_name)[i]],
                   label=clust_events.df$min_n[clust_events.df$model_name==unique(clust_events.df$model_name)[i]],
                   size=2.7)
      }
      
      print(pl6)
      
    }
    
  }
  
  # Select final K
  K_value = paste0(K_vec["manual_BIC_K"])
  D_value
  
  # Final selected parameters 
  print(paste0("Fitting GMM with ",K_value," components to a UMAP embedding with ", D_value, " dimensions"))
  
}

### Fit final models
{
  # Output: UMAP tranformation, fit GMM
  
  # UMAP transformation & save
  {
    # UMAP transform from orignal dimension to D dimensions 
    {
      suffix <- Sys.Date()
      
      # Run
      set.seed(10)
      umap_model <-uwot::umap(impt,scale=T,n_neighbors = 40,min_dist = 0.25,
                              learning_rate = 0.5,init="normlaplacian",ret_model = T,
                              n_components = D_value,verbose=T)
      
      rownames(umap_model$embedding)<-ext.dat[,"patientid"]
      
      
      # Create dir
      if(dir.exists("/path_to/project_folder/UMAP/")==F){
        dir.create("/path_to/project_folder/UMAP/")
        if(dir.exists("/path_to/project_folder/UMAP/model/")==F){
          dir.create("/path_to/project_folder/UMAP/model/")
        }
      }
      
      # Save transform 
      setwd("/path_to/project_folder/UMAP/model/")
      #uwot::save_uwot(umap_model,file = paste0("UMAP_",D_value,"D_model_export_",suffix,".rds"))
      
      # Load the transform
      umap_model<-uwot::load_uwot(file = paste0("UMAP_",D_value,"D_model_export_",suffix,".rds"))
      
      # Go back to main working folder
      setwd("/path_to/project_folder/")
    }
    
    # If D > 2, apply UMAP to visualise embedding in 2-D (not used for analysis, purely for visualistion)
    {
      if(D_value==2){
        print(paste0("2-D embedding chosen, no need to embedd to 2-D"))
      } else{
        print(paste0(D_value,"-D embedding chosen, embedding to 2-D"))
      }
      
      suffix <- Sys.Date()
      
      # Run
      set.seed(10)
      umap_model.visual <-uwot::umap(umap_model$embedding,
                                     scale=T,
                                     n_neighbors = 40,min_dist = 1,
                                     learning_rate = 0.5,
                                     init="normlaplacian",
                                     ret_model = T,
                                     n_components = 2,
                                     verbose=T,
                                     n_threads = 10)
      
      # Update patient id
      rownames(umap_model.visual$embedding)<-ext.dat[,"patientid"]
      
      # Save transform 
      setwd("/path_to/project_folder/UMAP/model/")
      #uwot::save_uwot(umap_model.visual,file = paste0("UMAP_",D_value,"D_to_2D_model_export_",suffix,".rds"))
      
      # Load the transform
      umap_model.visual<-uwot::load_uwot(file = paste0("UMAP_",D_value,"D_to_2D_model_export_",suffix,".rds"))
      
      # Go back to main working folder
      setwd("/path_to/project_folder/")
      
    }
    
  }
  
  # GMM model & save
  {
    suffix <- Sys.Date()
    
    # Base umap
    dat_umap <- umap_model$embedding
    
    # GMM mclust
    #.Random.seed <- seed_save
     set.seed(10) 
    GMM_model <- Mclust(umap_model$embedding,
                        G = K_value,
                        modelNames = mod_type,
                        initialization = list("hcpairs"))
    
    # Create dir
    if(dir.exists("/path_to/project_folder/GMM/")==F){
      dir.create("/path_to/project_folder/GMM/")
      if(dir.exists("/path_to/project_folder/GMM/model/")==F){
        dir.create("/path_to/project_folder/GMM/model/")
      }
    }
    
    # Save model
    #saveRDS(GMM_model,file = paste0("/path_to/project_folder/GMM/model/GMM_",D_value,"D_model_K",K_value,"_export_",suffix,".rds"))
    
  }
  
}

### Load fitted models from file 
{
  #umap_model<-uwot::load_uwot(file = "/path_to/project_folder/UMAP/model/UMAP_4D_model_export.rds")
  #umap_model.visual<-uwot::load_uwot(file = "/path_to/project_folder/UMAP/model/UMAP_4D_to_2D_model_export.rds")
  #GMM_model <- readRDS(paste0("/path_to/project_folder/GMM/model/GMM_4D_model_K29_export.rds"))
  
  K_value = 29
  D_value = 4
}

### Extract UMAP embedding and cluster labels 
{
  # Embedding
  embeddings <- list()
  embeddings[[paste0("dev_",D_value,"D")]]<-umap_model$embedding
  if(D_value!=2){
    embeddings[["dev_2d"]]<-umap_model.visual$embedding
  }
  
  # Classification 
  cluster_class<-list()
  cluster_class[["dev"]] <-GMM_model$classification
  
  # Store N 
  npts<- list()
  tmp<- data.frame(table(GMM_model$classification))
  tmp$Var1 <- as.numeric(as.character(tmp$Var1))
  tmp<-dplyr::left_join(data.frame("Var1"=seq(from=1,to=max(GMM_model$classification))),
                        tmp)
  npts[["dev"]] <- tmp 
  
  
}

### UMAP 2-D plot
{
  # By sex
  {
    if(D_value==2){
      embed <-data.frame(umap_model$embedding,stringsAsFactors = F)
    } else{
      embed <-data.frame(umap_model.visual$embedding,stringsAsFactors = F)
    }
    
    embed$group <- factor(dev.dat[,"sex"],levels=c(0,1),labels=c("Female","Male"))
    
    gg<-ggplot(embed,aes_string(x="X1",y="X2",colour="group"))+
      geom_point(size=0.3,alpha=0.7)+
      labs(x="UMAP1", y="UMAP2", color="Sex") +
      theme_bw() +
      theme(axis.line=element_blank(),axis.text.x=element_blank(),
            axis.text.y=element_blank(),axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank())+
      scale_color_manual(values=c("#00BFC4","#F8766D"))
    gg
  }
  
  # By cluster
  {
    # UMAP embeddign with group label and facet ID
    if(D_value==2){
      embed<-data.frame(umap_model$embedding,
                        "id"=rownames(umap_model$embedding),
                        "Group"=factor(cluster_class$dev))
    } else{
      embed<-data.frame(umap_model.visual$embedding,
                        "id"=rownames(umap_model.visual$embedding),
                        "Group"=factor(cluster_class$dev))
    }
    
    # Colours 
    colour_vec<-c("#32CD32","#F37FB8","#409388","#CE8BAE","#B23648",
                  "#ADD8E6","#D46E7E","#7E486B","#79AA7A","#FFEC2B",
                  "#8D5B96","#E41A1C","#00B4F0","#3A85A8","#488846",
                  "#BD6253","#46A169","#EB7AA9","#C4625D","#D8B62E",
                  "#d6c624","#77777C","#4F6FA1","#E1712E","#A65628",
                  "#B392A3","#E984B9","#F2E631","#999999")
    names(colour_vec) <- seq(from=1,to=length(colour_vec),by=1)
    
    
    # Plot (single)
    gg <-ggplot(embed,aes_string(x="X1",y="X2",colour="Group"))+
      geom_point(alpha=0.3,size=0.3,colour="grey",aes(group=id),data=embed[,!colnames(embed)%in%c("Group")])+
      geom_point(alpha=0.7,size=0.3)+
      xlab("UMAP1") + ylab("UMAP2") +
      theme_classic()+
      labs(color="Cluster")+
      theme(axis.line=element_blank(),axis.text.x=element_blank(),
            axis.text.y=element_blank(),axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())+
      guides(size = F)+
      scale_colour_manual(values=colour_vec)
    
    print(gg)
    
    # Plot (facet by cluster)
    gg1<-gg +
      facet_wrap(~Group)+
      theme(legend.position = "none")
    
    #print(gg1) # Very slow if N large 
  }
}

### Cluster 28 day mortality rates 
{
  #For each cluster (with n. pateints >10):
  # Calculate the proportion of deaths within 28 days
  
  # Data
  {
    # Development 
    dev_tmp <- data.frame("event"=factor(dev.dat$death28,levels=c(0,1),labels=c("censored","death")),
                          "cluster"=cluster_class$dev)
    
  }
  
  n_boot = 1000
  
  # Rank by proportion 
  {
    # Basic rank 
    test_prop<-dev_tmp %>%
      group_by(cluster,event) %>%
      summarise("n"=n()) %>% 
      mutate(prop = n / sum(n))
    tmp1<- test_prop[test_prop$event=="death",]
    tmp1<-tmp1[order(tmp1$prop),]
    tmp1 <- as.data.frame(tmp1)
    tmp1$rank <- c(1:length(tmp1$cluster))
    
    # Boostrap rank 
    rank_list<-list()
    boot_list <- bootstrap(dev_tmp,n=n_boot)
    for(b in 1:n_boot){
      dev_tmp.samp <- as.data.frame(boot_list$strap[[i]])
      
      test_prop2<-dev_tmp.samp %>%
        group_by(cluster,event) %>%
        summarise("n"=n()) %>% 
        mutate(prop = n / sum(n))
      tmp<- test_prop2[test_prop2$event=="death",]
      tmp<-tmp[order(tmp$prop),]
      tmp <- as.data.frame(tmp)
      tmp$rank <- c(1:length(tmp$cluster))
      tmp$samp <- rep(b,length(tmp$cluster))
      rank_list[[b]]<-tmp[,c("samp","cluster","rank")] 
    }
    tmp2<-melt(rank_list,id=c("cluster","rank"))
    tmp2$variable<-NULL
    tmp2$L1<-NULL
    tmp2<-dcast(tmp2,cluster~value,value.var = "rank")
    tmp2$median_rank <- apply(tmp2[,-1],1,median)
    tmp2 <- data.frame("cluster"=tmp2$cluster,"median_rank"=tmp2$median_rank)
    tmp2<-tmp2[order(tmp2$median_rank),]
    
    rank_df <- left_join(tmp1,tmp2)
    rank_df$cluster <- factor(rank_df$cluster,
                              levels = rank_df$cluster,
                              labels= rank_df$cluster,
                              ordered=T)
  }
  
  # Data frame 
  {
    # DF for proportion of deaths by cluster
    prop_df<-data.frame(matrix(nrow = 1,ncol=7))
    colnames(prop_df)<-c("cohort","cluster","cluster_n","cluster_n_event",
                         "cluster_prop_event","lower","upper")
    prop_list <- list()
  }
  
  # Check clusters have >10 observations
  dev_clust <- npts[["dev"]]
  dev_clust<-dev_clust[dev_clust$Freq>=10,"Var1"]
  
  # Proportions with bootstrapped CI
  for(i in 1:length(dev_clust)){
    
    # Set this cluster's data 
    clust <- dev_clust[i]
    dev_tmp.filt <-dev_tmp[which(dev_tmp$cluster==clust),]
    print(paste0("cluster: ",clust))
    
    # Check if uniform outcome (e.g. all death)
    dev_error = dim(table(as.numeric(dev_tmp.filt$event)))==1
    
    # If uniform outcome, skip bootstrap proportion
    if(dev_error==T){
      
      prop <- ifelse(unique(dev_tmp.filt$event)=="censored",0,1)
      
      prop_df.tmp <- prop_df
      prop_df.tmp$cohort = "Development"
      prop_df.tmp$cluster = clust
      prop_df.tmp$cluster_n = sum(table(dev_tmp.filt$event))
      prop_df.tmp$cluster_n_event = table(dev_tmp.filt$event)["death"]
      prop_df.tmp$cluster_prop_event = prop
      prop_df.tmp$lower = prop
      prop_df.tmp$upper = prop
      prop_df.tmp$cluster_prop_event_edit = prop
      prop_df.tmp$lower_edit = prop
      prop_df.tmp$upper_edit = prop
      
    } 
    
    # Proportion of events
    if(dev_error==F){ 
      
      # Calcaulte proportion & 95% CI
      dev_tmp.ci_prop <- ci_proportion(as.numeric(dev_tmp.filt$event)-1, type="bootstrap",R=n_boot)
      
      # store 
      prop_df.tmp <- prop_df
      prop_df.tmp$cohort = "Development"
      prop_df.tmp$cluster = clust
      prop_df.tmp$cluster_n = sum(table(dev_tmp.filt$event))
      prop_df.tmp$cluster_n_event = table(dev_tmp.filt$event)["death"]
      prop_df.tmp$cluster_prop_event = dev_tmp.ci_prop$estimate
      prop_df.tmp$lower = dev_tmp.ci_prop$interval[1]
      prop_df.tmp$upper = dev_tmp.ci_prop$interval[2]
    }
    
    prop_list[[i]] <- prop_df.tmp
  }
  
  # Combine to DF
  prop_df<-prop_list[[1]]
  for(i in 2:length(prop_list)){
    prop_df<-rbind(prop_df,prop_list[[i]])
  }
  
  # Output
  cluster_outcome.dev <- prop_df 
  
  # Add rank 
  cluster_outcome.dev$cluster <- factor(cluster_outcome.dev$cluster,levels=rank_df$cluster,ordered = T)
  cluster_outcome.dev<-cluster_outcome.dev %>% left_join(.,rank_df[,c("cluster","rank")]) %>% as.data.frame(.)
  
  # Check these agree
  cluster_outcome.dev[cluster_outcome.dev$cluster==1,]
  test_prop[test_prop$event=="death"&test_prop$cluster==1,]
  
  # Create dir
  if(dir.exists("/path_to/project_folder/outcome/")==F){
  #  dir.create("/path_to/project_folder/outcome/")
  }
  
  suffix <- Sys.Date()
  #write.csv(cluster_outcome.dev,file = paste0("/path_to/project_folder/outcome/cluster_outcome_dev_",suffix,".csv"))

}

### Cluster characterisation 
{
  
  # Input 
  {
    # Vector of categorical predictors for association tests
    # Continuous variables should first be categorised into ordinal factors - e.g. age group, urea low / high
    wrdcld_pred <- c(cat_pred)
    #wrdcld_pred <- c(cat_pred,"tempcat","heartratecat","systolicbpcat","diastolicbpcat","respscat")

    # Set P value threshold 
    pval_thresh <- 0.001
    
    # Complete data 
    impt<-complete(miDat.dev,1)[,wrdcld_pred] 
    
    # Edit labels 
    impt$sex <- ifelse(impt$sex==1,"Male","Female")
    impt$sex<-factor(impt$sex)
    
    # Vector showing cluster label per observation
    cluster<-factor(cluster_class$dev,
                    levels=seq(from=1,to=max(cluster_class$dev)),
                    ordered = T)
    
    # Vector of unique cluster labels
    K <- unique(cluster)
    
    # Add prefix to cluster labels
    levels(cluster) <- paste0("clust_",levels(cluster))
    
  }
  
  # Association tests
  {
    # Code written by Dr. Thomas Taverner (with minor modfication)
    # https://research.birmingham.ac.uk/en/persons/thomas-taverner
    
    # For details on thresholds see Chen, H., Cohen, P., & Chen, S. (2010). How Big is a Big Odds Ratio? Interpreting the Magnitudes of Odds Ratios in Epidemiological Studies. Communications in Statistics - Simulation and Computation, 39(4), 860-864. doi:10.1080/03610911003650383
    
    
    or_tabs <- list()
    pval_tabs <- list()
    for(var in colnames(impt)){
      print(var); flush.console()
      fac_levs <-levels(factor(impt[,var]))
      or_tmp <- pval_tmp <- matrix(NA, length(fac_levs), length(K))
      rownames(or_tmp) <- rownames(pval_tmp) <- fac_levs
      colnames(or_tmp) <- colnames(pval_tmp) <- paste0("clust_",K)
      for(facn in 1:length(fac_levs)){
        x <- 1*(impt[[var]]==fac_levs[facn])
        for(j in 1:length(K)){
          k <- K[j]
          k <- paste0("clust_", k)
          y <- 1*(cluster==k)
          ftest <- fisher.test(x, y)
          or_tmp[facn, j] <- as.numeric(ftest$estimate)
          pval_tmp[facn, j] <- ftest$p.value
        }
      }
      or_tabs[[var]] <- or_tmp
      pval_tabs[[var]] <- pval_tmp
    }
    
    effect_size_list <- sapply(colnames(impt), function(var){
      ortab <- or_tabs[[var]]
      pvaltab <- pval_tabs[[var]]
      newtab <- array("", dim(ortab))
      dimnames(newtab) <- dimnames(ortab)
      
      #0.518793793415168 1.24415459395877 1.90359895098359
      newtab[log(ortab) > 1.9 & newtab == "" & pvaltab < pval_thresh] <- "high"
      newtab[log(ortab) > 1.2 & newtab == "" & pvaltab < pval_thresh] <- "med"
      newtab[log(ortab) > 0.5 & newtab == "" & pvaltab < pval_thresh] <- "low"
      
      newtab
    }, simplify=FALSE)
    
    freq_list <- list()
    size <- c(high=">1.9", med="1.2-1.9", low="0.5-1.2")
    for(i in 1:length(K)){
      k <- K[i]
      clustname <- paste0("clust_", k)
      for(var in names(effect_size_list)){
        tab <- effect_size_list[[var]]
        if(length(rownames(tab)) == 2 && identical(rownames(tab), c("0", "1"))){
          tab <- tab[-1,,drop=FALSE]
          z <- tab[1, clustname, drop=FALSE]
          if(z != "") freq_list[[clustname]] <- c(freq_list[[clustname]],
                                                   setNames(size[z], var))
        } else {
          for(rn in rownames(tab)){
            z <- tab[rn, clustname]
            if(z != "") freq_list[[clustname]] <- c(freq_list[[clustname]],
                                                     setNames(size[z], paste0(var, "_", rn)))
          }
        }
      }
    }
    out_list <- freq_list
   
  } 
  
  # Export 
  {
    for(j in 1:length(out_list)){
      tmp = out_list[[j]]
      tmp <- data.frame(melt(tmp))
      tmp$value <- factor(tmp$value,levels = c("0.5-1.2",
                                               "1.2-1.9",
                                               ">1.9"),
                          labels = c("0.5-1.2",
                                     "1.2-1.9",
                                     ">1.9"),ordered = T)
      tmp$var <- rownames(tmp)
      #tmp <- tmp[rev(order(tmp$value)),]
      tmp <- tmp[order(tmp$var),]
      tmp$var <- NULL
      out_list[[j]] <- tmp
      names(out_list)[j]<-names(out_list)[j]
    }
    
    clust_order = gsub("clust_",replacement = "",names(out_list))
    clust_order = order(as.numeric(clust_order))
    
    out_list <- out_list[clust_order]
    names(out_list)<-gsub("clust_",replacement = "Cluster ",names(out_list))
    
    for(i in 1:length(out_list)){
      out_list[[i]]$var <- rownames(out_list[[i]])
    }
    
    # Create dir
    if(dir.exists("/path_to/project_folder/cluster_characterisation/")==F){
      dir.create("/path_to/project_folder/cluster_characterisation/")
    }
    
    # Save 
    #erer::write.list(out_list,file=paste0("/path_to/project_folder/cluster_characterisation/cluster_associations_dev_",suffix,".csv"),row.names=T)
   
    
    
    
    
    
  }
  
  # Continuous variable distributions 
  {
    # add dat and filter clusters 
    df <- data.frame(dev.dat,"cluster"=cluster_class$dev)

    # Order by outcome rank 
    df$cluster <- factor(df$cluster,levels = levels(cluster_outcome.dev$cluster),ordered=T)
    df$cluster<- droplevels(df$cluster)
    
    axis_rotate_theme<-
      theme(axis.text.x = element_text(angle = 90,vjust=0.5, hjust=1),
            axis.text.y = element_text(angle = 90,vjust=0.5, hjust=0.5))
    
    plot_list <- list()
    for(i in 1:length(c(cts_pred))){
      pred<-c(cts_pred)[i]
      whisker_plot<-ggplot(df,aes_string(x="cluster",
                                         y=paste0(pred)))+
        geom_boxplot(varwidth=T,outlier.shape = NA) +
        axis_rotate_theme + theme_bw()
      whisker_plot<-whisker_plot+labs(y=paste0(pred),x="4-D Cluster")
      plot_list[[pred]]<-whisker_plot
    }
    
    print(plot_list$age)
  }
  
  # Comorbidity heatmap 
  {
    # Data
    {
      comorbid_tmp <- c("charlsondementiaicd10", "cancer", "CVD", "COPD_SleepApnea_Asthma", "dm_WithAndWithoutComp")
      
      # add dat and filter clusters 
      df <- data.frame(dev.dat,"cluster"=cluster_class$dev)
      
      # Order by outcome rank 
      df$cluster <- factor(df$cluster,levels = levels(cluster_outcome.dev$cluster),ordered=T)
      df$cluster<- droplevels(df$cluster)
      
      # Na comorbid assume absent
      for(i in 1:length(comorbid_tmp)){
        df[is.na(df[,comorbid_tmp[i]]),comorbid_tmp[i]]<-0
      }
      
    }
    
    # Count proportions 
    {
      cluster_n <- length(levels(df$cluster))
      count_out <- matrix(nrow = cluster_n, ncol = length(comorbid_tmp))
      colnames(count_out) <- comorbid_tmp
      rownames(count_out) <- levels(df$cluster)
      for (i in 1:cluster_n) {
        idx <- df$cluster == levels(df$cluster)[i]
        for (j in 1:length(comorbid_tmp)) {
          vec = df[which(idx==T),which(colnames(df)==comorbid_tmp[j])]
          if(sum(vec)==0){
            count_out[i, j] <- 0
          } else {
            count_out[i, j] <- sum(vec)/length(vec)
          }
          colnames(count_out)[j]=comorbid_tmp[j]
        }
        rownames(count_out)[i]=levels(df$cluster)[i]
      }
      colnames(count_out)[colnames(count_out)=="dm_WithAndWithoutComp"]="Diabetes (any)"
      colnames(count_out)[colnames(count_out)=="COPD_SleepApnea_Asthma"]="COPD/Sleep apnoea/Asthma"
      
      # format
      prop1 <-
        count_out %>%
        data.table::as.data.table(keep.rownames = "cluster") %>%
        melt(id.vars = "cluster", value.name = "prop") %>%
        .[,cluster := factor(cluster, levels = levels(df$cluster),ordered=T)] %>%
        .[,variable := factor(variable, levels = sort(colnames(count_out)))]
      
      
      # Decile
      decile <- quantile(prop1$prop, probs = seq(.1, .9, by = .1))
      prop1$dec <- ntile(prop1$prop, 10)
      
      prop1$variable<-droplevels(prop1$variable)
      prop1$variable <- factor(prop1$variable,ordered=T)
    }
    
    # Plot 
    {
      pl <- prop1%>% 
        mutate(.,"cut"=cut(prop,breaks = c(-Inf,0.2,0.4,0.6,0.8,Inf)))%>%
        ggplot(aes(cluster,variable))+
        geom_tile(aes(fill = factor(cut,labels=c("0-19","20-39","40-59",
                                                 "60-79","80-100"))),
                  colour = "white")+
        geom_text(aes(label=dec,colour = factor(cut,labels=c("0-19","0-19","20-39",
                                                             "20-39","20-39")))) +
        scale_fill_brewer(palette = 3)+
        scale_colour_manual(values=c(RColorBrewer::brewer.pal(n = 5,"Blues")[5],"white"),guide = FALSE)+
        guides(fill=guide_legend(title="Proportion (%)"))+
        labs(x="Cluster",y="Comorbidity",title = "Development")+
        theme(axis.text.x = element_text(vjust = 0.5, hjust=0.5))+theme_bw()
      print(pl)
    }
  }
  
}
  
#----------------------------------------------------------------#

# VALIDATION

### Validation data
{
  # Read in data (wave 1)
  {
    # W1
    # Incomplete data, wide format, with additional variables such as outcome (.csv file)
    w1.dat<- read.csv("/path_to/project_folder/csv/wave_1_data.csv")
    
    # Imputed (mice::mids object,.RDS object)
    miDat.w1<-readRDS(file = "/path_to/project_folder/Imputed/wave1_miData.RDS")
  
    # Patient filter (vector with N observations of 1 (keep) or 0 (remove),.csv file)
    w1_eligible <- read.csv("/path_to/project_folder/csv/wave1_eligible_patients.csv")[,-1]
    
    # Update column names 
    impt.tmp <- complete(miDat.w1,include = T,action = "long")
    colnames(impt.tmp)<-gsub("_0$","",colnames(impt.tmp))  # Remove day name
    miDat.w1<-as.mids(impt.tmp)
    
  }  
  
  # Read in data (wave 2)
  {
    # W2
    # Incomplete data, wide format, with additional variables such as outcome (.csv file)
    w2.dat<- read.csv("/path_to/project_folder/csv/wave_2_data.csv")
    
    # Imputed (mice::mids object,.RDS object)
    miDat.w2<-readRDS(file = "/path_to/project_folder/Imputed/wave2_miData.RDS")
   
    # Patient filter (vector with N observations of 1 (keep) or 0 (remove),.csv file)
    w2_eligible <- read.csv("/path_to/project_folder/csv/wave2_eligible_patients.csv")[,-1]
    
    # Update column names 
    impt.tmp <- complete(miDat.w2,include = T,action = "long")
    colnames(impt.tmp)<-gsub("_0$","",colnames(impt.tmp))  # Remove day name
    miDat.w2<-as.mids(impt.tmp)
  }
  
  # Checks
  {
    # w1
    {
      impt.tmp<-complete(miDat.w1,1) # For cts dat
      # Filter to UMAP pred
      impt.tmp<-impt.tmp[,umap_pred]
      # Filter to eligible pts
      impt.tmp<-impt.tmp[w1_eligible==1,]
      # Check either continuous or binary
      head(impt.tmp)
      str(impt.tmp)
      for(i in 1:length(umap_pred)){
        vec <- impt.tmp[,umap_pred[i]]
        name <- umap_pred[i]
        if(name %in%cat_pred){
          if(all(vec==1|vec==0)==F){
            print(paste0("issue with cat var: ",name))
          } else{
            print(paste0("cat var: ",name))
            #print(summary(vec))
          }
        }
        if(name %in%cts_pred){
          vec = as.numeric(vec)
          if(anyNA(vec)==T){
            print(paste0("issue with cts var: ",name))
          }else{
            print(paste0("cts var: ",name))
            #print(summary(vec))
          }
        }
      }
      summary(impt.tmp)
    }
    
    # w2
    {
      impt.tmp<-complete(miDat.w2,1) # For cts dat
      # Filter to UMAP pred
      impt.tmp<-impt.tmp[,umap_pred]
      # Filter to eligible pts
      impt.tmp<-impt.tmp[w2_eligible==1,]
      # Check either continuous or binary
      head(impt.tmp)
      str(impt.tmp)
      
      for(i in 1:length(umap_pred)){
        vec <- impt.tmp[,umap_pred[i]]
        name <- umap_pred[i]
        if(name %in%cat_pred){
          if(all(vec==1|vec==0)==F){
            print(paste0("issue with cat var: ",name))
          } else{
            print(paste0("cat var: ",name))
            #print(summary(vec))
          }
        }
        if(name %in%cts_pred){
          vec = as.numeric(vec)
          if(anyNA(vec)==T){
            print(paste0("issue with cts var: ",name))
          }else{
            print(paste0("cts var: ",name)) 
            #print(summary(vec))
          } 
        }
      }
      summary(impt.tmp)
    }
    
    
  } 
  
}

### Transform and cluster new observations 
{
  # UMAP transformation onto embedding learned from the development cohort
  # Classification into clusters using fitted GMM  
  # Repeated for each imputation - mode classifacation used as final cluster classifcation

  # Wave 1
  {
    
    # List for UMAP embedding of each imputation 
    proj_list<-list() # This is for the embedding used in clustering 
    if(D_value!=2){
      print(paste0(D_value,"-D embedding for clustering, also embedding to 2-D for visualisation"))
      proj_list_2d<-list() # 2D 
    } 
    
    # Data frame for cluster classifcation of each individual per imputation
    mi_cluster_class <- data.frame(matrix(ncol=miDat.w1$m,
                                          nrow=length(w1_eligible)))
    
    # Classify observations based off each imputation 
    for(j in 1:miDat.w1$m){
      
      dat_name = paste0("imp_",j) 
      
      impt.tmp<- complete(miDat.w1,j)
      impt.tmp <- impt.tmp[,umap_pred]
      impt.tmp <- impt.tmp[w1_eligible==1,]
      
      # Make numeric 
      impt.tmp <- data.frame(apply(impt.tmp,2,function(x){as.numeric(as.character(x))}),stringsAsFactors = F)
      
      rownames(impt.tmp) <- paste0(as.character(w1.dat[w1_eligible==1,"Patient.ID"]))
      
      
      # Check predictors levels all match 
      #colnames(impt.tmp)<-gsub("_[0-9]","",colnames(impt.tmp))
      
      # Transform & cluster
      {
        set.seed(10)
        UMAP_t <- uwot::umap_transform(impt.tmp,umap_model)
        rownames(UMAP_t) <- rownames(impt.tmp)
        # Store transform 
        proj_list[[dat_name]] <- UMAP_t
        
        # Make cluster predictions
        pred <- predict.Mclust(newdata = UMAP_t, # This transformation 
                               object = GMM_model) # Fit model 
        names(pred$classification)<-rownames(UMAP_t)
        
        mi_cluster_class[,j]<-pred$classification
        rownames(mi_cluster_class)<-rownames(UMAP_t)
        
        #mi_cluster_prob[[j]]<-pred$z
        #rownames(mi_cluster_prob[[j]])<-rownames(UMAP_t)
      }
      
      # 2D transform 
      if(D_value!=2){
        
        set.seed(10)
        UMAP_t.2d <- uwot::umap_transform(UMAP_t,umap_model.visual)
        rownames(UMAP_t.2d)<-rownames(impt.tmp)
        
        # Store transform 
        proj_list_2d[[dat_name]] <- UMAP_t.2d
      } 
      
    }
    
    # Store UMAP embeddings
    embeddings[[paste0("w1_",D_value,"D_list")]]<-proj_list
    if(D_value!=2){
      # Store 2D embeddings  
      embeddings[[paste0("w1_2D_list")]]<-proj_list_2d
    }
    
    # Get cluster classification by taking mode across multiple imputations
    mi_cluster_class$mode_clust<-apply(mi_cluster_class,1,getmode)
    mi_cluster_class$n_diff<-apply(mi_cluster_class,1,function(x)length(unique(x)))
    summary(mi_cluster_class$n_diff)
    
    # Store 
    cluster_class[["w1"]] <- mi_cluster_class$mode_clust
    names(cluster_class[["w1"]]) <- rownames(mi_cluster_class)
    
    # Count per cluster
    tmp<- data.frame(table(mi_cluster_class$mode_clust))
    tmp$Var1 <- as.numeric(as.character(tmp$Var1))
    tmp<-dplyr::left_join(data.frame("Var1"=seq(from=1,to=max(GMM_model$classification))),
                          tmp)
    npts[["w1"]] <- tmp
    
    
  }  
  
  # Wave 2
  {
    
    # List for UMAP embedding of each imputation 
    proj_list<-list() # This is for the embedding used in clustering 
    if(D_value!=2){
      print(paste0(D_value,"-D embedding for clustering, also embedding to 2-D for visualisation"))
      proj_list_2d<-list() # 2D 
    } 
    
    # Data frame for cluster classifcation of each individual per imputation
    mi_cluster_class <- data.frame(matrix(ncol=miDat.w2$m,
                                          nrow=length(w2_eligible)))
    
    # Classify observations based off each imputation 
    for(j in 1:miDat.w2$m){
      
      dat_name = paste0("imp_",j) 
      
      impt.tmp<- complete(miDat.w2,j)
      impt.tmp <- impt.tmp[,umap_pred]
      impt.tmp <- impt.tmp[w2_eligible==1,]
      
      # Make numeric 
      impt.tmp <- data.frame(apply(impt.tmp,2,function(x){as.numeric(as.character(x))}),stringsAsFactors = F)
      
      rownames(impt.tmp) <- paste0(as.character(w2.dat[w2_eligible==1,"Patient.ID"]))
      
      
      # Check predictors levels all match 
      #colnames(impt.tmp)<-gsub("_[0-9]","",colnames(impt.tmp))
      
      # Transform & cluster
      {
        set.seed(10)
        UMAP_t <- uwot::umap_transform(impt.tmp,umap_model)
        rownames(UMAP_t) <- rownames(impt.tmp)
        # Store transform 
        proj_list[[dat_name]] <- UMAP_t
        
        # Make cluster predictions
        pred <- predict.Mclust(newdata = UMAP_t, # This transformation 
                               object = GMM_model) # Fit model 
        names(pred$classification)<-rownames(UMAP_t)
        
        mi_cluster_class[,j]<-pred$classification
        rownames(mi_cluster_class)<-rownames(UMAP_t)
        
        #mi_cluster_prob[[j]]<-pred$z
        #rownames(mi_cluster_prob[[j]])<-rownames(UMAP_t)
      }
      
      # 2D transform 
      if(D_value!=2){
        
        set.seed(10)
        UMAP_t.2d <- uwot::umap_transform(UMAP_t,umap_model.visual)
        rownames(UMAP_t.2d)<-rownames(impt.tmp)
        
        # Store transform 
        proj_list_2d[[dat_name]] <- UMAP_t.2d
      } 
      
    }
    
    # Store UMAP embeddings
    embeddings[[paste0("w2_",D_value,"D_list")]]<-proj_list
    if(D_value!=2){
      # Store 2D embeddings  
      embeddings[[paste0("w2_2D_list")]]<-proj_list_2d
    }
    
    # Get cluster classification by taking mode across multiple imputations
    mi_cluster_class$mode_clust<-apply(mi_cluster_class,1,getmode)
    mi_cluster_class$n_diff<-apply(mi_cluster_class,1,function(x)length(unique(x)))
    summary(mi_cluster_class$n_diff)
    
    # Store 
    cluster_class[["w2"]] <- mi_cluster_class$mode_clust
    names(cluster_class[["w2"]]) <- rownames(mi_cluster_class)
    
    # Count per cluster
    tmp<- data.frame(table(mi_cluster_class$mode_clust))
    tmp$Var1 <- as.numeric(as.character(tmp$Var1))
    tmp<-dplyr::left_join(data.frame("Var1"=seq(from=1,to=max(GMM_model$classification))),
                          tmp)
    npts[["w2"]] <- tmp
    
    
  } 
  
  # Join
  tmp<-cbind.data.frame(npts)
  tmp<-tmp[,c(1,2,4,6)]
  tmp
  npts.df <- tmp
}

### Cluster silhouette 
{
  plot_colours<- gg_color_hue(3)
  
  # Development
  {
    sil.dev <- clust_sil(data=embeddings$dev_4D,
                         labels=cluster_class$dev)
    sil.dev <- sil.dev$sil_summary
    sil.dev$S_hat <- NULL
    sil.dev$cohort <- "Development"
  }
  
  # Wave 1 
  {
    sil.w1 <- mi_clust_sil(data_list = embeddings$w1_4D_list,
                           labels = cluster_class$w1)
    sil.w1$cohort <- "Wave 1"
    colnames(sil.w1) <-gsub("pool_","",colnames(sil.w1))
  }
  
  # Wave 2
  {
    sil.w2 <- mi_clust_sil(data_list = embeddings$w2_4D_list,
                           labels = cluster_class$w2)
    sil.w2$cohort <- "Wave 2"
    colnames(sil.w2) <-gsub("pool_","",colnames(sil.w2))
  }
  
  # Bind
  sil_out <- rbind(sil.dev,sil.w1)
  sil_out <- rbind(sil_out,sil.w2)
  
  # Rank 
  {
    sil_rank <- sil_out %>%
      filter( cohort=="Development") %>%
      select(cluster,Q_hat)
    sil_rank <-sil_rank[order(sil_rank$Q_hat),]
    sil_rank$rank <- c(length(sil_rank$cluster):1)
    sil_rank$cluster <- factor(sil_rank$cluster,levels = sil_rank$cluster,
                               labels= sil_rank$cluster,
                               ordered=T)
    sil_rank[order(sil_rank$rank,decreasing = F),]
  }
  
  # Plot
  {
    # Plot proportions 
    tmp <- sil_out %>% filter(N>=10)
    tmp$cluster <- factor(tmp$cluster,levels=sil_rank$cluster,ordered = T)
    gg_sil<-   ggplot(aes(x=cluster, y=Q_hat),
                      dat=tmp[,]) +
      geom_hline(yintercept=0, linetype="11", colour="grey60") +
      geom_point(position=position_dodge(-0.8),
                 aes(colour=cohort)) +
      geom_pointrange(aes(x=cluster, ymax=upper, ymin=lower,
                          group=cohort, colour=cohort),
                      na.rm=TRUE,size=0.3,
                      position=position_dodge(width = -0.8))+
      scale_color_manual(values = plot_colours,labels=c("Development","Wave 1","Wave 2"))+
      coord_flip(clip = "off") +
      theme_classic() +labs(x="Cluster",colour="Cohort")+
      theme(legend.position = "top")+ylab("Average silhouette width")
    
    print(gg_sil)
  } 
  
}

### Cluster 28 day mortality rates & log ratios 
{
  # For each cluster (with n. pateints >10):
  # Calculate the proportion of deaths within 28 days
  # Calculate a log ratio between cohorts (1 K bootstrap resamples)
  
  cohort_colours<- gg_color_hue(3)
  
  # Data
  {
    # Event status and cluster label 
    
    # Development 
    dev_tmp <- data.frame("event"=factor(dev.dat$death28,levels=c(0,1),labels=c("censored","death")),
                          "cluster"=cluster_class$dev)
    
    # W1
    w1_tmp <- data.frame("event"=factor(w1.dat[w1_eligible==1,"death28"],levels=c(0,1),labels=c("censored","death")),
                         "cluster"=cluster_class$w1)
    
    # W2
    w2_tmp <- data.frame("event"=factor(w2.dat[w2_eligible==1,"death28"],levels=c(0,1),labels=c("censored","death")),
                         "cluster"=cluster_class$w2)
    
    # Aim:
    test_prop<-rbind(cbind(dev_tmp,cohort="dev"),
                     cbind(w1_tmp,cohort="w1"),
                     cbind(w2_tmp,cohort="w2")) %>%
      group_by(cohort,cluster,event) %>%
      summarise("n"=n()) %>% 
      mutate(prop = n / sum(n))
    test_prop[test_prop$event=="death",]
  }
  
  # Data frame 
  {
    # DF for proportion of deaths by cluster
    prop_df<-data.frame(matrix(nrow = 1,ncol=7))
    colnames(prop_df)<-c("cohort","cluster","cluster_n","cluster_n_event",
                         "cluster_prop_event","lower","upper")
    prop_df.blank<-prop_df
    
    # Empty df for difference   
    diff_df<-data.frame(matrix(nrow = 1,ncol=11))
    colnames(diff_df)<-c("cohort","cluster","cluster_n","odds_ratio","lower","upper","test",
                         "odds_ratio_w1w2","lower_w1w2","upper_w1w2","test_w1w2")
    diff_df.blank<-diff_df
    
    
  }
  
  # Clusters
  {
    # Get proportion for all clusters in development with >=10 observations 
    dev_clust <- npts.df[npts.df$dev.Freq>=10,"dev.Var1"]
    
    # Note clusters with <10 observations in validation cohorts 
    w1_clust <- npts.df[npts.df$w1.Freq<10|is.na(npts.df$w1.Freq),"dev.Var1"]
    w1_clust <- w1_clust[!is.na(w1_clust)]
    w2_clust <- npts.df[npts.df$w2.Freq<10|is.na(npts.df$w2.Freq),"dev.Var1"]
    w2_clust <- w2_clust[!is.na(w2_clust)]
    
  }
  
  n_boot = 1000
  
  # Repeat for each cluster
  for(i in 1:length(dev_clust)){
    
    # Set this cluster's data 
    {
      clust <- dev_clust[i]
      dev_tmp.filt <-dev_tmp[which(dev_tmp$cluster==clust),]
      w1_tmp.filt<-w1_tmp[which(w1_tmp$cluster==clust),]
      w2_tmp.filt<-w2_tmp[which(w2_tmp$cluster==clust),]
      print(paste0("cluster: ",clust))
      
    }
    
    # Check which can be compared
    {
      # Check error (low N or uniform outcome)
      w1_error = NULL
      w2_error = NULL
      
      # Low N w1
      if(clust%in%w1_clust){
        w1_error = T
        print(paste0("cluster: ",clust," <10 pts W1"))
      } else{
        w1_error = dim(table(as.numeric(w1_tmp.filt$event)))==1 
        # Single outcome w1
        if(w1_error==T){
          print(paste0("cluster: ",clust," >10 pts, single outcome, W1"))
        } else{
          print(paste0("cluster: ",clust," >10 pts, mixed outcome, W1"))
        }
      }
      
      # Low N w1
      if(clust%in%w2_clust){
        w2_error = T
        print(paste0("cluster: ",clust," <10 pts W2"))
      } else{
        w2_error = dim(table(as.numeric(w2_tmp.filt$event)))==1 
        # Single outcome w1
        if(w2_error==T){
          print(paste0("cluster: ",clust," >10 pts, single outcome, W2"))
        } else{
          print(paste0("cluster: ",clust," >10 pts, mixed outcome, W2"))
        }
      }
      
      # Either W1/W2
      if(w1_error==T|w2_error==T){
        w1w2_error=T
      } else{
        w1w2_error=F
      }
      
    }
    
    # W1 proportion, boostrap CI and logratio with development
    {
      
      # Skipping bootstrap if either low N or single outcome 
      if(w1_error==T){
        
        # If single outcome but >10
        if(clust%in%w1_clust==F){
          prop <- ifelse(unique(w1_tmp.filt$event)=="censored",0,1)
        } 
        
        # If <10 
        if(clust%in%w1_clust==T){
          prop <- NA
        } 
        
        w1.prop_df <- prop_df.blank
        w1.prop_df$cohort = "w1"
        w1.prop_df$cluster = clust
        w1.prop_df$cluster_n = sum(table(w1_tmp.filt$event))
        w1.prop_df$cluster_n_event = table(w1_tmp.filt$event)["death"]
        w1.prop_df$cluster_prop_event = prop
        w1.prop_df$lower = prop
        w1.prop_df$upper = prop
        
        
        w1.diff_df <- diff_df.blank
        w1.diff_df$cohort = "w1"
        w1.diff_df$cluster = clust
        w1.diff_df$cluster_n = sum(table(w1_tmp.filt$event))
        
      } 
      
      # Bootstrap proportion and log ratio 
      if(w1_error==F){
        
        # Bootstrap proportion
        w1_tmp.ci_prop <- ci_proportion(as.numeric(w1_tmp.filt$event)-1, type="bootstrap",R=n_boot)
        
        # store 
        w1.prop_df <- prop_df.blank
        w1.prop_df$cohort = "w1"
        w1.prop_df$cluster = clust
        w1.prop_df$cluster_n = sum(table(w1_tmp.filt$event))
        w1.prop_df$cluster_n_event = table(w1_tmp.filt$event)["death"]
        w1.prop_df$cluster_prop_event = w1_tmp.ci_prop$estimate
        w1.prop_df$lower = w1_tmp.ci_prop$interval[1]
        w1.prop_df$upper = w1_tmp.ci_prop$interval[2]
        
        # Comparison W1 to training 
        chi_tbl <- cont_tbl(x=w1_tmp.filt$event,y=dev_tmp.filt$event)
        w1_tmp.ftest <-fisher.test(chi_tbl,B=n_boot)
        
        # store 
        w1.diff_df <- diff_df.blank
        w1.diff_df$cohort = "w1"
        w1.diff_df$cluster = clust
        w1.diff_df$cluster_n = sum(table(w1_tmp.filt$event))
        w1.diff_df$odds_ratio = w1_tmp.ftest$estimate
        w1.diff_df$lower = w1_tmp.ftest$conf.int[1]
        w1.diff_df$upper = w1_tmp.ftest$conf.int[2]
        w1.diff_df$test = w1_tmp.ftest$p.value
        
      }
    } 
    
    # W2 proportion, boostrap CI and logratio with W1 & development
    {
      # Skipping bootstrap if either low N or single outcome 
      if(w2_error==T){
        
        # If single outcome but >10
        if(clust%in%w2_clust==F){
          prop <- ifelse(unique(w2_tmp.filt$event)=="censored",0,1)
        } 
        
        # If <10 
        if(clust%in%w2_clust==T){
          prop <- NA
        } 
        
        w2.prop_df <- prop_df.blank
        w2.prop_df$cohort = "w2"
        w2.prop_df$cluster = clust
        w2.prop_df$cluster_n = sum(table(w2_tmp.filt$event))
        w2.prop_df$cluster_n_event = table(w2_tmp.filt$event)["death"]
        w2.prop_df$cluster_prop_event = prop
        w2.prop_df$lower = prop
        w2.prop_df$upper = prop
        
        
        w2.diff_df <- diff_df.blank
        w2.diff_df$cohort = "w2"
        w2.diff_df$cluster = clust
        w2.diff_df$cluster_n = sum(table(w2_tmp.filt$event))
        
      } 
      
      # Bootstrap proportion and log ratio 
      if(w2_error==F){
        
        # Bootstrap proportion
        w2_tmp.ci_prop <- ci_proportion(as.numeric(w2_tmp.filt$event)-1, type="bootstrap",R=n_boot)
        
        # store 
        w2.prop_df <- prop_df.blank
        w2.prop_df$cohort = "w2"
        w2.prop_df$cluster = clust
        w2.prop_df$cluster_n = sum(table(w2_tmp.filt$event))
        w2.prop_df$cluster_n_event = table(w2_tmp.filt$event)["death"]
        w2.prop_df$cluster_prop_event = w2_tmp.ci_prop$estimate
        w2.prop_df$lower = w2_tmp.ci_prop$interval[1]
        w2.prop_df$upper = w2_tmp.ci_prop$interval[2]
        
        # Comparison w2 to training 
        chi_tbl <- cont_tbl(x=w2_tmp.filt$event,y=dev_tmp.filt$event)
        w2_tmp.ftest <-fisher.test(chi_tbl,B=n_boot)
        
        # store 
        w2.diff_df <- diff_df.blank
        w2.diff_df$cohort = "w2"
        w2.diff_df$cluster = clust
        w2.diff_df$cluster_n = sum(table(w2_tmp.filt$event))
        w2.diff_df$odds_ratio = w2_tmp.ftest$estimate
        w2.diff_df$lower = w2_tmp.ftest$conf.int[1]
        w2.diff_df$upper = w2_tmp.ftest$conf.int[2]
        w2.diff_df$test = w2_tmp.ftest$p.value
        
        if(w1w2_error==F){
          # Comparison w2 to training 
          chi_tbl <- cont_tbl(x=w2_tmp.filt$event,y=w1_tmp.filt$event)
          w2_tmp.w1_ftest <-fisher.test(chi_tbl,B=n_boot)
          w2.diff_df$odds_ratio_w1w2 = w2_tmp.w1_ftest$estimate
          w2.diff_df$lower_w1w2 = w2_tmp.w1_ftest$conf.int[1]
          w2.diff_df$upper_w1w2 = w2_tmp.w1_ftest$conf.int[2]
          w2.diff_df$test_w1w2 = w2_tmp.w1_ftest$p.value
        }
        
      } 
    }
    
    # Combine prop df 
    prop_df <- rbind(prop_df,w1.prop_df,w2.prop_df)
    prop_df<-prop_df[which(is.na(prop_df$cohort)==F),] # Na made by R bind, drop
    
    # Combine diff df
    diff_df<-rbind(diff_df,w1.diff_df,w2.diff_df)
    diff_df<-diff_df[which(is.na(diff_df$cohort)==F),] # Na made by R bind, drop
    
  }
  
  # Ouput: proporiton of deaths by cluster with CI
  {
    #
    cluster_outcome.all <- rbind(cluster_outcome.dev[,!colnames(cluster_outcome.dev)=="rank"],prop_df)
    cluster_outcome.all<-cluster_outcome.all %>% left_join(.,rank_df[,c("cluster","rank")]) %>% as.data.frame(.)
    
    # Add label
    cluster_outcome.all$label=paste(round(cluster_outcome.all$cluster_prop_event*100),"%",sep="")
    cluster_outcome.all$label[grep("NA%",cluster_outcome.all$label)]<-""
  }
  
  # Output: odds ratios by cluster between cohorts 
  {
    cluster_OR <- diff_df
    
    cluster_OR$cluster <- factor(cluster_OR$cluster,levels=rank_df$cluster,ordered = T)
    cluster_OR<-cluster_OR %>% left_join(.,rank_df[,c("cluster","rank")]) %>% as.data.frame(.)
    
    #cluster_OR$signif <- ifelse(is.na(cluster_OR$test),"",ifelse(cluster_OR$test<=0.05,"*",""))
    
  }
  
  
  print(head(cluster_outcome.all))
  print(head(cluster_OR))
  
}
