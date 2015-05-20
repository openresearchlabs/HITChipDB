
#' Phylogeneticenrichments
#'
#' Description: Computes the dependence of the higher phylogenetic
#' level groups of the given list x, which can be a vector/data frame of either
#' oligos, species, or level2 groups. Computes only the enrichments
#' if onlyEnrich=T (neglects the unexpected disappearences of the
#' groups in the given list). Outputs only the enrichments below the
#' given p-value threshold p.th.
#' Requires phylogeny.info=the first 4 columns from phylogenyprofile,
#' origlevel=the level from which the list is, and maplevel= the
#' level to which the mapping should be made and the enrichments of
#' which computed.
#' Outputs a four component list: p-values from Fisher's tests,
#' input tables for tests, the actual tests outputs, and the full
#' phylogenetic mapping information for the given list.
#' Changes: Version 1 computed two-tailed p-values also for single
#' occurrences on orig
#' v2: two-tailed p-values, but not for single occurrences
#' v3: one-tailed p-values (enrichments), if wanted, outputting only
#' enrichments under the given p-value
#'
#' Arguments:
#'   @param x x
#'   @param phylogeny.info phylogeny.info
#'   @param origlevel origlevel
#'   @param maplevel maplevel 
#'   @param onlyEnriched onlyEnriched
#'   @param p.th p-value threshold
#'
#' Returns:
#'   @return enrichments
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

Phylogeneticenrichments <- function(x, phylogeny.info, origlevel = colnames(phylogeny.info)[3], maplevel = "L1", onlyEnriched = T, p.th = 0.05)
{
  
  ## Convert to character vector
  x <- as.character(x)
    
  ## collect origlevel groups for all x items
  origlevel.split <- split(phylogeny.info,phylogeny.info[[origlevel]])
  x.origlevelGroups <- origlevel.split[x]
  
  ## vector (inX and inXind) showing which items from the origlevel of
  ## phylogeny.info are in the x list
  origlevel.ugroups <- names(origlevel.split)
  inXind <- match(x, origlevel.ugroups)
  inX <- rep(F,length(origlevel.ugroups))
  inX[inXind] <- T

  ## Collect the full phyloinfo
  phyloM <-  x.origlevelGroups[[1]]
  nolGroups <- length(x.origlevelGroups)
  if(nolGroups>1)
    for(i in 2:nolGroups){
      phyloM <- rbind(phyloM,x.origlevelGroups[[i]])
    }
  
  if(length(x)>2){
    
    ## compute enrichments of maplevel groups in x
    e <- c()
    pvals <- c()
    estimates <- c()
    tables <- list()
    maplevel.ugroups <- as.character(unique(phyloM[[maplevel]]))

    for(g in maplevel.ugroups){

      ## compute in which maplevel groups the given item occurs
      inMaplevelGroup <- unlist(lapply(origlevel.split,
                                         function(y){is.element(g,y[[maplevel]])
                                                   }))

      if(sum(inMaplevelGroup)>0){
        if(onlyEnriched)
          tmp <- try(fisher.test(inX, inMaplevelGroup, alternative="g"))
        else
          tmp <- try(fisher.test(inX, inMaplevelGroup, alternative="t"))

          e[g] <- list(tmp)
          pvals[g] <- tmp$p.value          
          tables[g] <- list(table(inX,inMaplevelGroup))
          estimates[g] <- tmp$estimate

      }

      pvals.adjusted <- p.adjust(pvals, method = "BH")
      inds <- (pvals.adjusted < p.th)
      e <- e[inds]
      pvals <- pvals.adjusted[inds]
      tables <- tables[inds]
      estimates <- estimates[inds]

    }

    ## Return enrichments i.e. Fisher's exact tests in table form

    if (length(pvals)>0)
       list(pval.adj=t(rbind(pvals,estimates)), 
           tables=tables, tests=e, phylogeny.info=phyloM) 
    else
       list(pvalues=pvals, tables=tables, tests=e, phylogeny.info=phyloM) 

  } else
    list(pvalues=1, tables=NULL, tests=NULL, phylogeny.info=phyloM) 
}

#' Description: BAGGED RDA Bootstrap solutions that follows the
#' Jack-knife estimation of PLS by Martens and Martens, 2000.  Solves
#' rotational invariance of latent space by orthogonal procrustes
#' rotations
#'
#' Arguments:
#'   @param X a matrix, samples on columns, variables (bacteria) on rows.
#'   @param Y vector with names(Y)=rownames(X), for example
#'   @param boot Number of bootstrap iterations
#'
#' Returns:
#'   @return List with items:
#'   	     loadings: bagged loadings
#' 	     scores: bagged scores
#' 	     significance: significances of X variables
#' 	     etc. TBA.	    
#'
#' @export
#' @importFrom vegan scores
#' @importFrom vegan procrustes
#'
#' @examples # NOT RUN data(peerj32); x <- as.matrix(peerj32$microbes)[1:20, 1:6]; y <- rnorm(nrow(x)); names(y) <- rownames(x); res <- Bagged.RDA(x, y , boot = 5)
#'
#' @references See citation("microbiome") 
#' @author Contact: Jarkko Salojarvi \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

Bagged.RDA <- function(X, Y, boot = 1000){

  ## Jarkko Salojärvi 7.8.2012
  ##  #17.8.2012 fixed problem with multiclass RDA  
  rda <- NULL

   if (is.numeric(boot)){
      class.split=split(names(Y),Y)

      boot=replicate(boot,unlist(sapply(class.split,function(x) sample(x,length(x),replace=T))),simplify=F)
   }
   nboot=length(boot)
   n.lev=length(levels(Y))
   TT=scores(rda(t(X)~Y),choices=1:max(n.lev-1,2),display="sites")
   nRDA=ncol(TT) 
   # rotate 
   rotateMat=function(M,TT,x){
      M.rot=procrustes(TT[x,],M)
      return(M.rot$Yrot)
   }
   nearest.centers=function(xx,cc){
     nC=nrow(cc)
     nS=nrow(xx)
     MM=as.matrix(dist(rbind(cc,xx)))[1:nC,(nC+1):(nC+nS)]
     if (nS>1)
       apply(MM,2,which.min)
     else
       which.min(MM)
   }
   # bootstrap loadings
   Tx=lapply(boot,function(x){ 
      nC=length(levels(Y[x]))
      M=rda(t(X[,x])~Y[x])
      # get scores
      TT.m=scores(M,choices=1:max(nC-1,2),display="sites")
      # bootstrap error
      testset=setdiff(colnames(X),x)
      err=NA
      if (length(testset)>0){
        Pr=predict(M,t(as.data.frame(X[,testset])),type="wa",model="CCA")
        centers=apply(TT.m,2,function(z) sapply(split(z,Y[x]),mean))
        if (nC==2)
           y.pred=apply(sapply(Pr,function(x) (x-centers[,1])^2),2,which.min)
        else
           y.pred=nearest.centers(Pr,centers) 
        err=mean(y.pred-as.numeric(Y[testset])!=0)
      }
      # procrustes rotation of scores
      TT.m=rotateMat(TT.m,TT,x)
      # solve loadings
      a=t(TT.m) %*% TT.m
      b=as.matrix(X[,x]-rowMeans(X[,x]))
      loadingsX=t(solve(a,t(b %*% TT.m)))
      list(loadingsX=loadingsX,err=err)
   })
   # significances   
   sig.prob=list()
   for (i in 1:nRDA){
     tmp=sapply(Tx,function(x) x$loadingsX[,i])
     sig.prob[[i]]=apply(tmp,1,function(x){ x1=sum(x>0)/length(x); x2=1-x1; min(x1,x2)})
    }
    names(sig.prob)=colnames(TT)
    # bagged estimates
    bagged.loadings=Tx[[1]]$loadingsX
    for (i in 2:nboot)
      bagged.loadings=bagged.loadings+Tx[[i]]$loadingsX
    bagged.loadings=bagged.loadings/nboot
    colnames(bagged.loadings)=colnames(TT)

    # solve scores
    a=t(bagged.loadings) %*% bagged.loadings
    b=as.matrix(X-rowMeans(X))
    bagged.scores=t(solve(a,t(bagged.loadings) %*% b))
    colnames(bagged.scores)=colnames(TT)
    # Group centers 
    Group.center=apply(bagged.scores,2,function(x) sapply(split(x,Y),mean))
    err.t=mean((nearest.centers(bagged.scores,Group.center)-as.numeric(Y))!=0) 

    # bagged error
    #err.t=mean((Y-prediction(X,Y,bagged.loadings,bagged.loadings.Y,bagged.proj,mean(Y),rowMeans(X)))^2)
    #random.pred=Y-rmultinom(length(Y),1,
    err.random=replicate(nboot,mean((as.numeric(Y)-sample(as.numeric(Y)))!=0))
    bagged.error=mean(sapply(Tx, function(x) x$err),na.rm=T)
    R=(bagged.error-err.t)/(mean(err.random)-err.t)
    R=max(min(R,1),0)
    w=.632/(1-.368*R)
    bagged.R2=(1-w)*bagged.error+w*err.t
    #bagged.R2=1-((1-w)*err.t+w*err.b)/mean(Y^2)
    can.cor.R=apply(X,1,function(x) cor(x,bagged.scores))^2
    Rsquare=rowSums(can.cor.R)/sum(diag(var(t(X))))
    names(Rsquare)=colnames(bagged.scores)
    Rsquare.variable=t(can.cor.R/apply(X,1,var))
    colnames(Rsquare.variable)=colnames(bagged.scores)
    list(loadings=bagged.loadings,scores=bagged.scores,significance=sig.prob,error=bagged.R2,group.centers=Group.center,bootstrapped=Tx,err.random=mean(err.random),err.significance=sum(err.random>bagged.R2)/nboot,R2=Rsquare,R2.variables=Rsquare.variable)
}


#' Description: BAGGED RDA feature selection
#'
#' Arguments:
#'   @param X a matrix, samples on columns, variables (bacteria) on rows.
#'   @param Y vector with names(Y)=rownames(X), for example
#'   @param sig.thresh signal p-value threshold, default 0.1
#'   @param nboot Number of bootstrap iterations
#'   @param verbose verbose
#'
#' Returns:
#'   @return List with items:
#'   	     loadings: bagged loadings
#' 	     scores: bagged scores
#' 	     significance: significances of X variables
#' 	     group.centers: group centers on latent space
#' 	     bootstrapped: bootstrapped loadings
#' 	     data: data set with non-significant components dropped out.
#'
#' @examples #
#' #  library(microbiome)
#' #  data(peerj32)
#' #  x <- t(peerj32$microbes)
#' #  y <- factor(peerj32$meta$time); names(y) <- rownames(peerj32$meta)
#' # Bag.res <- Bagged.RDA.Feature.Selection(x, y, sig.thresh=0.05, nboot=100)
#' # PlotBaggedRDA(Bag.res, y)
#'
#' @export
#'
#' @references See citation("microbiome") 
#' @author Contact: Jarkko Salojarvi \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities


Bagged.RDA.Feature.Selection <- function(X,Y,sig.thresh=0.1,nboot=1000, verbose = T){

  stop.run=F
  class.split=split(names(Y),Y)
  dropped=vector("character",nrow(X))
  X.all=X
  mean.err=rep(1,nrow(X))
  while(stop.run==F){
    boot=replicate(nboot,unlist(sapply(class.split,function(x) sample(x,length(x),replace=T))),simplify=F)
    Bag.res=Bagged.RDA(X,Y,boot=boot)
    min.prob=Bag.res$significance[[1]]
    if (length(levels(Y))>2){
      for (i in 1:nrow(X))
         min.prob[i]=min(sapply(Bag.res$significance,function(x) x[i]))
    }
    mean.err[nrow(X)]=Bag.res$error
    dropped[nrow(X)]=rownames(X)[which.max(min.prob)]
    if (verbose) {message(c(nrow(X), Bag.res$error))}
    if (nrow(X)>max(length(class.split),2))
      X=X[-which.max(min.prob),]
    else
      stop.run=T
  }
  dropped[1:length(class.split)]=rownames(X)[order(min.prob)[1:length(class.split)]]
  best.res=which.min(mean.err)

  Bag.res=Bagged.RDA(X.all[dropped[1:best.res],],Y,boot=boot)
  Bag.res$data=X.all[dropped[1:best.res],]
  Bag.res$Err.selection=mean.err
  Bag.res$dropped=dropped

  plot(mean.err,xlab="X dimension")
  points(best.res,mean.err[best.res],col="red")
  Bag.res
}




#' Description: Bagged RDA visualization
#'
#' Arguments:
#'   @param Bag.res Output from  Bagged.RDA.Feature.Selection function
#'   @param Y vector with names(Y)=rownames(X), for example
#'   @param which.bac TBA
#'   @param ptype TBA
#'   @param comp TBA
#'   @param cex.bac TBA
#'   @param plot.names Plot names
#'   @param group.cols TBA
#'   @param ... Other arguments to be passed
#'   
#'
#' Returns:
#'   @return TBA
#'
#' @examples #
#'   #library(microbiome)
#'   #data(peerj32)
#'   #x <- t(peerj32$microbes)
#'   #y <- factor(peerj32$meta$time); names(y) <- rownames(peerj32$meta)
#'   #Bag.res <- Bagged.RDA.Feature.Selection(x, y, sig.thresh=0.05, nboot=100)
#'   #PlotBaggedRDA(Bag.res, y)
#'
#' @export
#' @importFrom ade4 s.class
#'
#' @references See citation("microbiome") 
#' @author Contact: Jarkko Salojarvi \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

PlotBaggedRDA <- function(Bag.res, Y, which.bac = 1:nrow(Bag.res$loadings),
	           ptype="spider", comp=1:2, cex.bac=0.5, plot.names=T,
		   group.cols = as.numeric(unique(Y)),...){

  # NOTE LL20131004: plot.bagged.RDA could not be used as a name
  # unless a particular plot function is defined for the object
  # within the package. Otherwise there will be warnings in package build.

  scaled.loadings <- (Bag.res$loadings/max(abs(Bag.res$loadings)))[,comp]
  scaled.scores <- (Bag.res$scores/max(abs(Bag.res$scores)))[,comp]

  plot(rbind(scaled.scores,scaled.loadings),type="n",xlab=paste(names(Bag.res$R2)[1]," (",format(100*Bag.res$R2[1],digits=2),"%)",sep=""),ylab=paste(names(Bag.res$R2)[2]," (",format(100*Bag.res$R2[2],digits=2),"%)",sep=""))
  if (ptype=="spider")
    s.class(scaled.scores,factor(Y),grid=F,col=group.cols,cellipse=0.5,cpoint=0,add.plot=T)
  if (ptype=="hull"){

    ll=split(rownames(scaled.scores),Y)
    hulls=lapply(ll,function(ii) ii[chull(scaled.scores[ii,])])
    for (i in 1:length(hulls))
      polygon(scaled.scores[hulls[[i]],],border=group.cols[i])
  }   
  if (plot.names){
     text(scaled.scores,rownames(scaled.scores),cex=0.5,...)
  }else{
    points(scaled.scores,...)
  }
  text(scaled.loadings[which.bac,],rownames(scaled.loadings)[which.bac],cex=cex.bac)
}


probe.summarization <- function (probedata, taxonomy, levels, summarization.methods) {

  # Ensure we have all levels
  levels <- intersect(levels, colnames(taxonomy))

  # Oligo summarization
  finaldata <- list()
  finaldata$oligo <- probedata

  for (level in levels) {

    finaldata[[level]] <- list()

    for (method in summarization.methods) {

        message(paste(level, method))
	# For species/L1/L2 summarization use the filtered phylogeny: phylogeny.filtered!
    	summarized.log10 <- summarize.probesets(
					taxonomy = taxonomy,		
			    		  oligo.data = log10(probedata), 
      			       	          method = method, 
					   level = level)$summarized.matrix

        # Store the data in absolute scale					
        finaldata[[level]][[method]] <- 10^summarized.log10

    }
  }

  finaldata

}



#' Description: summarize.probesets.through.species
#'
#' Arguments:
#'   @param level summarization level
#'   @param taxonomy oligo - phylotype matching data.frame
#'   @param oligo.data preprocessed probes x samples data matrix in log10 domain
#'   @param method summarization method
#'   @param verbose verbose
#'
#' Returns:
#'   @return List with two elements: summarized.matrix (summarized data matrix in log10 scale) and probe.parameters (only used with rpa, probe-level parameter estimates)
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

summarize.probesets.through.species <- function (level, taxonomy, oligo.data, method, verbose = FALSE) {

  # Get species-level summary first
  message("Summarizing through species...")
  probeset.summaries <- summarize.probesets.species(taxonomy, oligo.data, method, verbose = verbose)

  species.matrix   <- 10^probeset.summaries$summarized.matrix
  probe.parameters <- probeset.summaries$probe.parameters 

  if (level == "species") {
    res <- list(summarized.matrix = log10(species.matrix), probe.parameters = probe.parameters)
    return(res)
  }

  # List all species for the given level (L0 / L1 / L2)")
  phylogroups <- levelmap(phylotypes = NULL, from = level, to = "species", taxonomy)

  # Remove specified phylogroups
  # phylogroups <- phylogroups[setdiff(names(phylogroups), rm.phylotypes[[level]])]

  summarized.matrix <- matrix(NA, nrow = length(phylogroups), ncol = ncol(oligo.data))
  rownames(summarized.matrix) <- sort(names(phylogroups))
  colnames(summarized.matrix) <- colnames(oligo.data)

  # Go through each phylogroup and summarize from species level
  for (pg in names(phylogroups)) {

    specs <- unique(phylogroups[[pg]])
    mat <- matrix(species.matrix[specs,], nrow = length(specs))

    if (method == "ave") { vec <- colMeans(mat) }
    if (method == "sum") { vec <- colSums(mat)  } 
    if (length(grep("rpa", method)) > 0) { vec <- colSums(mat) } # For RPA, use the sum for L1/L2

    summarized.matrix[pg, ] <- vec

  }

  list(summarized.matrix = log10(summarized.matrix), probe.parameters = probe.parameters)

}




#' Description: Probeset summarization with various methods.
#' 
#' Arguments:
#'   @param phylogeny.info oligo - phylotype matching data.frame
#'   @param oligo.data preprocessed probes x samples data matrix in log10 domain
#'   @param method summarization method
#'   @param verbose print intermediate messages
#'   @param rm.species Species to exclude
#'   @param probe.parameters Optional. If probe.parameters are given,
#'          the summarization is based on these and model parameters are not
#' 	    estimated. A list. One element for each probeset with the following probe vectors: 
#'	    affinities, variances
#' Returns:
#'   @return List with two elements: summarized.matrix (summarized data matrix in log10 scale) and probe.parameters (only used with rpa, probe-level parameter estimates)
#'
#' @export
#' @importFrom RPA d.update.fast 
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

summarize.probesets.species <- function (phylogeny.info, oligo.data, method, verbose = TRUE, rm.species = c("Victivallis vadensis"), probe.parameters = list()) {

  level <- "species"			    
  probesets <- retrieve.probesets(phylogeny.info, level = level)
  probesets <- probesets[setdiff(names(probesets), rm.species)]
  nPhylotypesPerOligo <- n.phylotypes.per.oligo(phylogeny.info, level) 

  # initialize
  summarized.matrix <- array(NA, dim = c(length(probesets), ncol(oligo.data)), 
  		    	      dimnames = list(names(probesets), colnames(oligo.data))) 

  probe.parameters <- list()
  if (method == "frpa") {
    if (verbose) {message("Loading pre-calculated preprocessing parameters")}
    rpa.hitchip.species.probe.parameters <- list()
    load(system.file("extdata/probe.parameters.rda", package = "HITChipDB"))
    probe.parameters <- rpa.hitchip.species.probe.parameters

    # Ensure we use only those parameters that are in the filtered phylogeny
    for (bac in names(probe.parameters)) {
      probe.parameters[[bac]] <- probe.parameters[[bac]][intersect(names(probe.parameters[[bac]]), probesets[[bac]])]
    }
  }

  for (set in names(probesets)) {

    # if (verbose) { message(set) }

    # Pick expression for particular probes
    probes <- probesets[[set]]

    # Pick probe data for the probeset: probes x samples
    # oligo.data assumed to be already in log10
    dat <- matrix(oligo.data[probes,], length(probes)) 
    rownames(dat) <- probes
    colnames(dat) <- colnames(oligo.data)

    if (method == "frpa") {

      # Summarize with pre-calculated variances
      vec <- d.update.fast(dat, probe.parameters[[set]])

    } else if (method == "rpa") {

      # RPA is calculated in log domain
      # Downweigh non-specific probes with priors with 10% of virtual data and
      # variances set according to number of matching probes
      # This will provide slight emphasis to downweigh potentially
      # cross-hybridizing probes
      res <- rpa.fit(dat, 
      	     		  alpha = 1 + 0.1*ncol(oligo.data)/2, 
			  beta  = 1 + 0.1*ncol(oligo.data)*nPhylotypesPerOligo[probes]^2)

      vec <- res$mu
      probe.parameters[[set]] <- res$tau2

    } else if (method == "ave") {

      vec <- log10(colMeans((10^dat), na.rm = T))

    } else if (method == "sum") {

      # Weight each probe by the inverse of the number of matching phylotypes
      # Then calculate sum -> less specific probes are downweighted
      # However, set the minimum signal to 0 in log10 scale (1 in original scale)!
      dat2 <- (10^dat) / nPhylotypesPerOligo[rownames(dat)]
      dat2[dat2 < 1] <- 1
      vec <- log10(colSums(dat2, na.rm = T))
      vec[which(vec == -Inf)] <- 0

    }
    
    summarized.matrix[set, ] <- vec 

  }

  list(summarized.matrix = summarized.matrix, probe.parameters = probe.parameters)
  
}


#' Description: summarize.probesets.directly
#'
#' Arguments:
#'   @param level summarization level
#'   @param phylogeny.info oligo - phylotype matching data.frame
#'   @param oligo.data preprocessed probes x samples data matrix in log10 domain
#'   @param method summarization method
#'   @param verbose verbose
#'
#' Returns:
#'   @return List with two elements: summarized.matrix (summarized data matrix in log10 scale) and probe.parameters (only used with rpa, probe-level parameter estimates)
#'
#' @export
#' @importFrom RPA d.update.fast rpa.fit
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

summarize.probesets.directly <- function (level, phylogeny.info, oligo.data, method, verbose = TRUE) {

  message("Direct summarization..")			     

  # Retrieve oligos for each taxonomic group
  probesets <- retrieve.probesets(phylogeny.info, level = level)

  # Number of phylotypes per oligo
  nPhylotypesPerOligo <- n.phylotypes.per.oligo(phylogeny.info, level) 

  # initialize
  summarized.matrix <- array(NA, dim = c(length(probesets),   ncol(oligo.data)), 
  		    	    dimnames = list(names(probesets), colnames(oligo.data))) 

  probe.parameters <- list()
  if (method == "frpa") {
    stop("fRPA not yet implement for direct summarization")
    if (verbose) {message("Loading pre-calculated preprocessing parameters")}
    rpa.hitchip.species.probe.parameters <- list()
    load(system.file("extdata/probe.parameters.rda", package = "HITChipDB"))
    probe.parameters <- rpa.hitchip.species.probe.parameters
  }

  for (set in names(probesets)) {
   
    # Pick data for particular probes
    probes <- probesets[[set]]

    # Pick probe data for the probeset: probes x samples
    # oligo.data assumed to be already in log10
    dat <- matrix(oligo.data[probes,], length(probes)) 
    rownames(dat) <- probes
    colnames(dat) <- colnames(oligo.data)

    if (method == "frpa") {

      message(method)

      # Summarize with pre-calculated variances
      vec <- d.update.fast(dat, probe.parameters[[set]])

    } else if (method == "rpa") {

      message(method)

      	    # RPA is calculated in log domain
     	    # Downweigh non-specific probes with priors with 10% of virtual data and
      	    # variances set according to number of matching probes
      	    # This will provide slight emphasis to downweigh potentially
      	    # cross-hybridizing probes
      	    res <- rpa.fit(dat, 
      	     		  alpha = 1 + 0.1*ncol(oligo.data)/2, 
			  beta  = 1 + 0.1*ncol(oligo.data)*nPhylotypesPerOligo[probes]^2)

      	    vec <- res$mu
      	    probe.parameters[[set]] <- res$tau2

   } else if (method == "rpa.with.affinities") {

      message(method)

     # Also include affinities in summarization

      	    res <- rpa.fit(dat, 
      	     		  alpha = 1 + 0.1*ncol(oligo.data)/2, 
			  beta  = 1 + 0.1*ncol(oligo.data)*nPhylotypesPerOligo[probes]^2, 
			  summarize.with.affinities = TRUE)

      	    vec <- res$mu
      	    probe.parameters[[set]] <- res$tau2

     } else if (method == "ave") {

       vec <- log10(colMeans((10^dat), na.rm = T))

     } else if (method == "sum") {

            # Weight each probe by the inverse of the number of matching phylotypes
      	    # Then calculate sum -> less specific probes are downweighted
      	    # However, set the minimum signal to 0 in log10 scale (1 in original scale)!
      	    dat2 <- (10^dat) / nPhylotypesPerOligo[rownames(dat)]
      	    dat2[dat2 < 1] <- 1
      	    vec <- log10(colSums(dat2, na.rm = T))
      	    vec[which(vec == -Inf)] <- 0

     } 

     summarized.matrix[set, ] <- vec 

  }

  list(summarized.matrix = summarized.matrix, probe.parameters = probe.parameters)  

}

#' Description: Generate toydata for the package
#'
#' Arguments:
#'   @param output.dir output directory name
#'
#' Returns:
#'   @return output file name
#'
#' @export
#' @importFrom microbiome read.profiling
#' @importFrom microbiome GetPhylogeny
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

GenerateSimulatedData <- function (output.dir) {

  data.directory <- system.file("extdata/", package = "microbiome")

  #phylogeny.info <- read.profiling(level = "phylogeny.info", data.dir = data.directory)
  phylogeny.info <- as.data.frame(GetPhylogeny("HITChip", "filtered"))

  oligo.matrix.nolog.simulated <- read.profiling(data.dir = data.directory, log10 = FALSE)
  N <- ncol(oligo.matrix.nolog.simulated)
  colnames(oligo.matrix.nolog.simulated) <- paste("Sample.", 1:N, sep = "")

  # Oligo summarization
  finaldata <- list()
  finaldata[["oligo"]] <- oligo.matrix.nolog.simulated
  levels <- c("species", "L2", "L1")
  for (level in levels) {
    finaldata[[level]] <- list()
    for (method in c("sum", "rpa")) {

        message(paste(level, method))
    	summarized.log10 <- summarize.probesets(phylogeny.info, log10(oligo.matrix.nolog.simulated), method = method, level = level)$summarized.matrix
      			       	          
        # Store the data in absolute scale					
        finaldata[[level]][[method]] <- 10^summarized.log10

    }
  }

  # Write summary matrices into the output directory
  outd <- WriteChipData(finaldata, output.dir, phylogeny.info)

  set.seed(344)
  metadata.simulated <- data.frame(list(
              sampleID = I(colnames(oligo.matrix.nolog.simulated)),
	      time = rep(1:4, 5),
              age = runif(N, 0, 100),
              bmi = runif(N, 20, 40),
	      subjectID = I(paste("subjectID", rep(1:4, 5), sep = "-")),
	      group = I(sample(paste("group", rep(1:4, 5), sep = "-"))),
              gender = I(sample(c("M", "F"), N, replace = TRUE)),
              diet = I(sample(c("Apricots", "Beverages", "Carrots"), N, replace = TRUE))))
  
  write.table(metadata.simulated, file = paste(output.dir, "/metadata.tab", sep = ""), quote = FALSE, sep = "\t")

  output.dir

}



#' Description: count
#'
#' For cross-hyb control
#'
#' Arguments:
#'   @param d TBA
#'
#' Returns:
#'   @return list
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Jarkko Salojarvi \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

count <- function(d){
   tabulate(d)
}



#' Description: summarize.oligos 
#'
#' Arguments:
#'   @param oligo.matrix oligo.matrix
#'   @param phylogeny.info phylogeny.info
#'   @param level taxonomic level
#'
#' Returns:
#'   @return list
#'
#' @export
#' @importFrom microbiome polish.phylogeny.info
#' @references See citation("microbiome") 
#' @author Contact: Jarkko Salojarvi \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

summarize.oligos <- function(oligo.matrix, phylogeny.info, level = "L2"){

  phylogeny.info <- polish.phylogeny.info(phylogeny.info)
		 
  oligo <- split(phylogeny.info$oligoID,phylogeny.info[,level])
  oligo <- lapply(oligo, function(x) unique(as.character(x)))

  Sim <- sapply(oligo,function(x){  
    if (length(x) > 1)
      return(colSums(oligo.matrix[x, ]))
    else
      return(oligo.matrix[x, ])
  })
  t(Sim)
}


#' Description: mixingMatrix
#'
#' For cross-hyb control
#'
#' Arguments:
#'   @param phylogeny.info phylogeny.info
#'   @param level taxonomic level
#'
#' Returns:
#'   @return oligos x phylotypes mixing matrix
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Jarkko Salojarvi \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

mixingMatrix <- function(phylogeny.info, level){

  M <- matrix(0,length(unique(phylogeny.info$oligoID)),length(unique(phylogeny.info[,level])),dimnames=list(sort(as.character(unique(phylogeny.info$oligoID))),sort(as.character(unique(phylogeny.info[,level])))))

  for (i in 1:nrow(phylogeny.info))
    M[as.character(phylogeny.info$oligoID[i]),as.character(phylogeny.info[i,level])]=1

  M <- apply(M, 2, function(x) x/sum(x))

  return(M)

}

#' Description: ngp: Non-negative probabilistic solution to probe mixing problem.
#'
#' Uses gamma prior to solve known cross-hybridisation problems. 
#' Perfect probe targets are given in phylogeny.info data frame. 
#'
#' Arguments:
#'   @param oligo.data oligo.data in absolute domain
#'   @param phylogeny.info oligo - phylotype mapping data frame
#'   @param level taxonomic level
#'   @param lambda - stregth of gamma prior. Default: 0.001.
#'   @param alpha - alpha parameter of gamma prior. Default: 1.
#'   @param beta - beta parameter of gamma prior. Default: 1.
#'
#' Returns:
#'   @return list
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Jarkko Salojarvi \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

ngp <- function(oligo.data, phylogeny.info, level, lambda=0.001, alpha=1,beta=1){

   # Avoid NOTE in buildcheck
   ginv <- NULL

   # oligos x phylotypes mixing matrix for the given level
   M <- mixingMatrix(phylogeny.info, level)
   coms <- intersect(rownames(oligo.data), rownames(M))
   M <- M[coms, ]  
   oligo.data <- oligo.data[coms, ]

   # starting guess using pseudoinverse
   W <- t(M)%*%M
   X <- t(oligo.data) %*% M
   A <- ginv(W)%*%t(X)
   A[which(A<0)] <- 0

   for (i in 1:ncol(oligo.data)){
      Acol=A[,i]
      for (j in 1:ncol(M)){  
         a=2*W[j,j]
         b=sum((W[j,-j]+W[-j,j])*Acol[-j])-2*X[i,j]-beta*lambda
         c=lambda*(alpha-1)
         if (a!=0){
           D=b^2-4*a*c
           if (D>0)
             A[j,i]=max((-b+sqrt(D)),(-b-sqrt(D)))/(2*a)
         }else
            A[j,i]=-c/b
     }
   }
   colnames(A) <- colnames(oligo.data)
   rownames(A) <- colnames(M)
   return(A)
}

#' Description: Deconvolution
#'
#' !OBSOLETE! Used nmf package (currently removed from CRAN) for solving cross-hyb problem.
#' Calls function ngp with standard arguments (SEE: ngp).
#'
#' Arguments:
#'   @param oligo.data oligo.data in absolute domain
#'   @param phylogeny.info oligo - phylotype mapping data frame
#'   @param level taxonomic level
#'   @param block.solution block.solution
#'   @param verbose verbose
#'   @param ... parameters to be passed to ngp function
#'
#'
#' Returns:
#'   @return list
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Jarkko Salojarvi \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

deconvolution.nonneg <- function(oligo.data, phylogeny.info, level, block.solution = T,verbose=F,...){
   ngp(oligo.data, phylogeny.info, level,...)
}

# --------------------------------------------------------------------



