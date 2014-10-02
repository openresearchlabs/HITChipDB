#' Core3D
#'
#' Description: Core visualization 3D
#'
#' Arguments:
#'  @param coreMat core matrix
#'  @param title title
#'  @param xlab X axis label
#'  @param cex.axis axis text size
#'
#' Returns:
#'  @return Used for its side effects
#'
#' @examples data(peerj32); 
#'          c3d <- Core3D(createCore(t(peerj32$microbes)))
#' @export 
#' 
#' @references 
#'   Jalanka-Tuovinen et al. Intestinal Microbiota in Healthy Adults: Temporal 
#'   Analysis Reveals Individual and Common Core and Relation to Intestinal 
#'   Symptoms. PLoS ONE 6(7):e23035, 2012.
#'   To cite the microbiome R package, see citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

Core3D <- function(coreMat, title = "Core microbiota", 
                   xlab = "Minimum Intensity", cex.axis = 0.7) {
    
    MinimumPrevalence <- as.numeric(colnames(coreMat))
    MinimumLogIntensity <- as.numeric(rownames(coreMat))
    tmp <- persp(MinimumLogIntensity, MinimumPrevalence, coreMat, 
                 theta = 60, phi = 5, 
        main = title, col = "light blue", axes = TRUE, ticktype = "detailed", 
        nticks = 9, 
        shade = 0.58, cex.axis = cex.axis, ylab = "Minimum Prevalence", 
        xlab = xlab, 
        zlab = "Core Size")
    
    return(NULL)
    
}

#' Description: core.which
#'
#' Arguments:
#'   @param data data matrix; phylotypes vs. samples
#'   @param intTr intTr
#'   @param prevalenceTr prevalenceTr
#'
#' Returns:
#'   @return TBA
#'
#' @references 
#'   A Salonen et al. The adult intestinal core microbiota is determined by 
#'   analysis depth and health status. Clinical Microbiology and Infection 
#'   18(S4):16 20, 2012. 
#'   To cite the microbiome R package, see citation('microbiome') 
#' @author Contact: Jarkko Salojarvi \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

core.which <- function(data, intTr, prevalenceTr) {
    d.bin <- data >= intTr
    prevalences <- rowSums(d.bin)
    nOTUs <- as.numeric(prevalences >= prevalenceTr)
    return(nOTUs)
}


#' bootstrap.microbecount
#'
#' Description: bootstrap.microbecount
#'
#' Arguments:
#'   @param D data
#'   @param Nsample sample size
#'   @param minprev minimum prevalence
#'   @param Nboot bootstrap sample size
#'   @param I.thr threshold
#'   @param ncore number of nodes for parallelization
#'
#' Returns:
#'   @return median microbe count in bootstrapped cores
#'
#' @examples \dontrun{data(peerj32); 
#' 	     tmp <- bootstrap.microbecount(t(peerj32$microbes),	Nboot = 5)}
#'
#' @export 
#' 
#' @references See citation("microbiome") 
#' @author Contact: Jarkko Salojarvi \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

bootstrap.microbecount <- function(D, Nsample = NULL, minprev = 1, Nboot = 1000, I.thr = 1.8, ncore = 1){

  if (is.null(Nsample)) {Nsample <- ncol(D)}

   boot <- replicate(Nboot,sample(ncol(D),Nsample,replace=T),simplify=F)

   # below:choose intensity such that there is at least one bacteria fulfilling prevalence criterion
   if (Nsample>1 && ncore > 1) {
     boot.which=mclapply(boot,function(x){ 
        Imax=max(apply(D[,x],1,min))
        Insty=runif(1,I.thr,Imax)
        sum(rowSums(D[,x]>Insty) >= minprev)
     }, mc.cores = ncore)
   } else if (Nsample>1 && ncore == 1) {
     boot.which=lapply(boot,function(x){ 
        Imax=max(apply(D[,x],1,min))
        Insty=runif(1,I.thr,Imax)
        sum(rowSums(D[,x]>Insty) >= minprev)
     }) 
   } else {
     boot.which=lapply(boot,function(x){ 
        Imax=max(D[,x])
        Insty=runif(1,I.thr,Imax)
        return(sum(D[,x]>=Insty))
     })
   }

   boot.prob <- as.matrix(as.data.frame(boot.which, check.names = F))
   t1 <- quantile(boot.prob, probs = c(0.05, 0.5, 0.95))
   #t1[2] <- mean(boot.prob)

   #print(t1)
   return(t1[2])
}

#' bootstrap.microbes
#'
#' Description: bootstrap.microbes
#'
#' Arguments:
#'   @param D data (phylotypes x samples)
#'   @param Nsample bootstrap sample size, default is the same size as data
#'   @param minPrev Lower limit for number of samples where microbe needs to exceed the intensity threshold for a 'present' call. 
#'   @param Nboot bootstrap sample size
#'   @param I.thr Lower limit for intensity threshold
#'   @param ncore number of nodes for parallelization - default 1
#'
#' Returns:
#'   @return data frame with microbes and their frequency of presence in the core
#'
#' @examples \dontrun{data(peerj32); 
#' 	     bs <- bootstrap.microbes(t(peerj32$microbes), Nboot = 5)}
#'
#' @export 
#' @import parallel
#' 
#' @references See citation("microbiome") 
#' @author Contact: Jarkko Salojarvi \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

bootstrap.microbes <- function(D, Nsample = NULL, minPrev = 2, Nboot = 1000, I.thr = 1.8, ncore = 1){

   if (is.null(Nsample)) {Nsample <- ncol(D)}

   boot <- replicate(Nboot, sample(ncol(D), Nsample, replace = T), simplify = F)

   # choose intensity such that there is at least one bacteria fulfilling prevalence criterion
   if (ncore > 1) {
     boot.which <- mclapply(boot, function(x){ 
       Prev = round(runif(1, minPrev, length(x)));
       Imax = max(apply(D[,x], 1, function(xx) quantile(xx, probs = (1 - Prev/length(x)))));
       Imax = max(I.thr, Imax); # Ensure Imax > I.thr, otherwise Insty gives NA's / LL 13.8.2012
       Insty = runif(1, I.thr, Imax);
       return(core.which(D[,x], Insty, Prev))
    }, mc.cores = ncore)
   } else {
     boot.which <- lapply(boot, function(x){ 
       Prev = round(runif(1, minPrev, length(x)));
       Imax = max(apply(D[,x], 1, function(xx) quantile(xx, probs = (1 - Prev/length(x)))));
       Imax = max(I.thr, Imax); # Ensure Imax > I.thr, otherwise Insty gives NA's / LL 13.8.2012
       Insty = runif(1, I.thr, Imax);
       return(core.which(D[,x], Insty, Prev))
    })
   }

   boot.prob <- rowSums(as.data.frame(boot.which))/Nboot

   df=data.frame(Microbe = rownames(D), Frequency = boot.prob)
   df=df[order(df$Frequency,decreasing=T),]

   mm=bootstrap.microbecount(D,Nsample = Nsample, minprev = minPrev, Nboot = Nboot, I.thr = I.thr, ncore = ncore)
   df$suggested.core=as.numeric(sapply(1:nrow(df),function(x) x<= mm))
   message("Bootstrapped median core size: ",mm)
   return(df)

}

#' plot_cumulative
#'
#' Description: plot_cumulative
#'
#' Arguments:
#'   @param d.sub d.sub
#'   @param writedir output directory
#'   @param fname output file name
#'   @param i.set i.set
#'   @param type plot type 
#'   @param ylim y axis limits
#'   @param phylogeny.info phylogeny.info matrix
#'
#' Returns:
#'   @return TBA
#'
#' @examples \dontrun{bs <- bootstrap.microbes(t(peerj32$microbes), Nboot = 5)
#' 	     	      plot_cumulative(bs, writedir = "./", "tmp-", 
#'		      			phylogeny.info = phylogeny.info)}
#'
#' @export 
#' 
#' @references 
#'   TBA.
#'   To cite the microbiome R package, see citation("microbiome") 
#' @author Contact: Jarkko Salojarvi \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

plot_cumulative <- function(d.sub, fname, i.set = NULL, type = "cumulative", ylim = NULL, phylogeny.info){

   PH.i <- unique(phylogeny.info[,1:2])
   PH.i <- PH.i[order(PH.i[,2]),]
   rownames(PH.i)=PH.i[,2]
   d.sub$Microbe <- PH.i[d.sub[,1],1]
   d.sub <- d.sub[order(d.sub[,2], decreasing = T),]

   if (is.null(i.set)) {
      i.set <- 1:length(levels(d.sub[,1]))
   }
   colmap <- colorRampPalette(c("Red", "Green","Blue"),space="rgb")(length(levels(d.sub[,1])))
   cnt <- 1;
   i.accept <- vector("logical",length(levels(d.sub[,1])))

   for (i in i.set){
     l.res <- as.numeric(d.sub[,1]==levels(d.sub[,1])[i])
     if (type=="cumulative"){

        out=cumsum(l.res)

     }
     if (type=="gsea"){

        l.res[l.res==0]=-1
        out=l.res
        for (j in 2:length(l.res))
           out[j]=max(l.res[j]+out[j-1],0)
     }
     t1=seq(max(d.sub[,2]),min(d.sub[,2]),-0.01)
     out=vector("numeric",length(t1))
     null.cum=matrix(NA,length(t1),3)
     for (j in 1:length(t1)){
       out[j]=sum(l.res*(d.sub[,2]>t1[j]))
       null.cum[j,]=quantile(replicate(1000,sum(sample(l.res,length(l.res))*(d.sub[,2]>=t1[j]))),probs=c(0.025,0.5,0.975))
     }
     yplot <- (out-null.cum[,2])/max(abs(out-null.cum[,2]))
     if (sum(out<null.cum[,1])>0 | sum(out>null.cum[,3])>0){
      if (cnt==1){
        if (is.null(ylim))
           plot(t1, yplot, type="l", main=paste(type,"prevalence of L1 taxa (",fname,")"), xlim=c(max(t1),min(t1)), col=colmap[i], ylab="proportion of total",xlab="Frequency")
        else
           plot(t1,yplot,type="l",main=paste(type,"prevalence of L1 taxa (",fname,")"),ylim=ylim,xlim=c(max(t1),min(t1)),col=colmap[i],ylab="proportion of total",xlab="Frequency")
      }else
         lines(t1,yplot,col=colmap[i])
      cnt=cnt+1;
      i.accept[i] <- TRUE
    }
   }
   legend(max(t1),1,levels(d.sub[,1])[which(i.accept==T)],fill=colmap[which(i.accept==TRUE)],cex=0.5)

   fname

}


