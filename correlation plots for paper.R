  #plot function for the paper - cowplot

# 17th July 2019 - Paul T

#-------------------------------------------------------------------------#
require(doBy)
library(tidyr)
require(tidyverse)
require(MASS)
require(stats)
library(matrixcalc)
library(ggpubr)

#-------------------------------------------------------------------------#

options(scipen = 999) #turn off scientific notation


#-------------------------------------------------------------------------#

gsca_sim_multigene_plot<-function(nsnp,nsnpeff,ngenes=1,nsub=120,ncases=100000,gpcorr=0.1,n2sim=3,cb=cb,cluster=cluster,phenCorr=phenCorr)
{
  gene_break<-floor(nsnp/ngenes)
  
  #set the correlation 
  gpcov<-gpcorr
  
  snpnames<-paste0('snp', 1:nsnp)
  maf<-runif(nsnp,.25,.5) #minor allele freq set to value from .25 to .5
  
  #Setup a dataframe to hold to population simulated data.
  mydata<-data.frame(matrix(nrow=ncases,ncol=(3+sum(nsnp))))
  
  
  myfilenames <- c('DatasetRand_N',
                   paste0('DatasetUncorr_',nsnpeff,'SNP_',10*gpcorr),
                   'Dataset4Block_N',
                   paste0('Dataset',ngenes,'Block_',nsnpeff,'SNP_',10*gpcorr))
  thisfile<-4 #User specified: Indicates which type of correlation structure in datafile (could use loop here)
  mydatafile<- myfilenames[thisfile]
  
  #SNP values are 0, 1 or 2, but we start with random normal numbers
  
  #-----------------------------------------------------------------
  #simulate the SNPS as correlated zscores, correlation is mycorr
  #-----------------------------------------------------------------
  #Where there is correlation, assume correlation depends on how close SNPs are, 
  #achieved by making size of correlation proportional to distance
  # (where sequence order is proxy for distance)
  
  h <-nsnp
  
  mymean<-rep(0,h) #vector of means to use in mvrnorm function
  mycorr<-0 #default is uncorrelated
  mycov2<-matrix(mycorr,nrow=h,ncol=h) #cov matrix for mvrnorm
  
  diag(mycov2)<-rep(1,h) #ones on diagonal of cov matrix
  
  MyRandCorrMat<-function(h=nsnp)
  {
    #R <- matrix(rbeta(h^2,2,2.5), ncol=h) 
    R <- matrix(rbeta(h^2,.23,1.23), ncol=h) 
    R <- (R * lower.tri(R)) + t(R * lower.tri(R)) 
    diag(R) <- 1
    return(R)
  }
  
  mycov2 <- MyRandCorrMat(nsnp)
  
  
  #-------------------------------------------------------------------------#
  
  mycov_PT<-matrix(mycorr,nrow=nsnp+n2sim,ncol=nsnp+n2sim)
  mycov_PT[1:nsnp,1:nsnp]<-mycov2
  m<-diag(n2sim)
  pheno_mat<-m[lower.tri(m)|upper.tri(m)]<-phenCorr
  mycov_PT[(nsnp+1):(nsnp+n2sim),(nsnp+1):(nsnp+n2sim)]<-pheno_mat
  
  #-------------------------------------------------------------------------#
  for(j in (nsnp+1):(nsnp+n2sim))
  {
    if(!nsnpeff==0){
      if(cluster==0){snpeff_pos<-sample(1:nsnp,nsnpeff,replace=F)} else {snpeff_pos<-c((nsnp-(nsnpeff-1)):nsnp)}
      
      mycov_PT[snpeff_pos,j] <- gpcov#*ifelse(rbinom(1,1,c(0.5,0.5))==1,-1,1)
      mycov_PT[j,snpeff_pos] <- t(mycov_PT[snpeff_pos,j])
    }
  }
  mycov<-mycov_PT
  
  #-------------------------------------------------------------------------#
  
  for(i in 1:nsnp)
  {
    for(j in 1:nsnp)
    {
      mycov[i,j] <- ifelse(rbinom(1,1,prob=0.5)==1,mycov[i,j],(-1*mycov[i,j]))
    }
  }
  
  for(i in (nsnp+1):(nsnp+n2sim))
  {
    for(j in 1:nsnp)
    {
      mycov[i,j] <- ifelse(rbinom(1,1,prob=0.5)==1,mycov[i,j],(-1*mycov[i,j]))
      mycov[j,i] <- ifelse(rbinom(1,1,prob=0.5)==1,mycov[j,i],(-1*mycov[j,i]))
    }
  }
  
  
  mycov[upper.tri(mycov)] <- t(mycov)[upper.tri(mycov)]
  
  diag(mycov)<-rep(1,dim(mycov)[1])
  
  #-------------------------------------------------------------------------#
  
  #then set diagonal values to 1 for all
  #diag(mycov)<-rep(1,n2sim)
  if(matrixcalc::is.positive.definite(mycov)==FALSE)
  {mycov<-Matrix::nearPD(mycov,keepDiag=TRUE)$mat}
  
  mymean<-rep(0,dim(mycov)[1])
  mydata=mvrnorm(n = ncases, mymean, mycov)
  
  mydata<-as.data.frame(mydata)
  
  colnames(mydata)[1:nsnp]<-snpnames
  colnames(mydata)[(nsnp+1):(nsnp+3)]<-c('NwdRepPheno','LangPheno','NeurodevPheno')
  
  #-------------------------------------------------------------------------#
  #Convert gene scores to integers: 0-2 for autosomal
  
  firstcol<-1
  lastcol<-nsnp
  p<-c(0,0,0) #initialise a vector to hold p values for different N alleles
  for (i in 1:nsnp){
    p[1]<-(1-maf[i])^2
    p[2]<-2*(1-maf[i])*maf[i]
    p[3]<-maf[i]^2
    
    #now use p-values to assign z-score cutoffs that convert to 0,1,2 or 3 minor alleles
    temp<-mydata[,i]
    w<-which(temp<qnorm(p[1]))
    mydata[w,i]<-0
    w<-which(temp>qnorm(p[1]))
    mydata[w,i]<-1
    w<-which(temp>qnorm(p[2]+p[1]))
    mydata[w,i]<-2
    
  }
  
  myr<-cor(mydata) #to check you have desired correlation structure, View(myr)
  
  
  #----------------------------------------------------------------------#

    cor.data<-as.matrix(myr)
    cor.data[lower.tri(cor.data)] <- NA
    
    cor_tri_N <- as.data.frame(cor.data) %>% 
      mutate(Var1 = factor(row.names(.), levels=row.names(.))) %>% 
      gather(key = Var2, value = value, -Var1, na.rm = TRUE, factor_key = TRUE) 
    
    
    g1<-ggplot(data = cor_tri_N, aes(Var2, Var1, fill = value)) + 
      geom_tile()+scale_fill_gradient2(low = "blue",mid="white", high = "red")+theme_bw()+ theme(axis.text.x = element_text(angle = 90, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.title = element_blank(),legend.title=element_blank(),text=element_text(size=14),title=element_text(size=10))+geom_hline(yintercept = (nsnp+0.5),colour="grey")+geom_vline(xintercept = (nsnp+0.5),colour="grey")+ggtitle("Random LD pattern")
    
 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    didat<-read.table('/Users/paulthompson/Dropbox/project SCT analysis/GSCA validation/SLICmergedPROcounttest.txt',header=TRUE,stringsAsFactors=FALSE)
    #Gene 1 - rs10255943 to rs2396765 (13 SNPs) foxp2
    #Gene 2 - rs10271363 to rs17170742, 23 SNPs cntnap2
    #Gene 3 - rs16957277 to rs16957385, 7 SNPs (edited)
    #Gene 4 - rs11149652 to rs4782962, 18 SNPs atp2c2
    #didat<-didat[,3:length(didat)] #strip off first 2 columns
    #didat<-didat[c(1:52,54:71,73:75,77:135),] #remove those with no phenotype data
    
    
    #for SNPs replace missing data with 1
    #didat[is.na(didat)]<-1 #just to get this to run I have replaced NA with 1
    
    #now read in 2nd dataset with more SNPs
    didat2<-read.table('/Users/paulthompson/Dropbox/project SCT analysis/GSCA validation/GSCA_SLIC_paul.txt',header=TRUE,stringsAsFactors=FALSE)
    
    myrow<-intersect(didat$ID,didat2$IND )
    
    alldat<-dplyr::filter(didat,ID%in%myrow)
    alldat2<-dplyr::filter(didat2,IND%in%myrow)
    #checked visually these are in the same order, so can combine with cbind
    
    alldat<-cbind(alldat[3:63],alldat2[5:292])
    
    #######################################################################
    refdat<-data.frame(SNP=colnames(alldat2[5:292]),gene=rep(NA,288))
    
    
    refdat$gene[1:99]<-"foxp2"
    refdat$gene[100:185]<-"cntnap2"
    refdat$gene[186:221]="ubr1"
    refdat$gene[222:288]="atp2c2"
    ########################################################################
    
    bigtab<-cor(alldat,use='pairwise.complete.obs')
    
   
      #Setup a dataframe to hold to population simulated data.
      mydata<-data.frame(matrix(nrow=ncases,ncol=(3+sum(nsnp))))
      
      
      myfilenames <- c('DatasetRand_N',
                       paste0('DatasetUncorr_',nsnpeff,'SNP_',10*gpcorr),
                       'Dataset4Block_N',
                       paste0('Dataset',ngenes,'Block_',nsnpeff,'SNP_',10*gpcorr))
      thisfile<-4 #User specified: Indicates which type of correlation structure in datafile (could use loop here)
      mydatafile<- myfilenames[thisfile]
      
      #SNP values are 0, 1 or 2, but we start with random normal numbers
      
      #-----------------------------------------------------------------
      #simulate the SNPS as correlated zscores, correlation is mycorr
      #-----------------------------------------------------------------
      #Where there is correlation, assume correlation depends on how close SNPs are, 
      #achieved by making size of correlation proportional to distance
      # (where sequence order is proxy for distance)
      
      h <-nsnp
      LDbase<-0
      
      mymean<-rep(0,h) #vector of means to use in mvrnorm function
      mycorr<-0 #default is uncorrelated
      mycov2<-matrix(mycorr,nrow=h,ncol=h) #cov matrix for mvrnorm
      
      diag(mycov2)<-rep(1,h) #ones on diagonal of cov matrix
      if(thisfile>2){ #only add correlation for conditions where thisfile is 3 or 4
        thisr<-mycov2<-cor(alldat[,sample(1:dim(all.data)[2],h)],use='pairwise.complete.obs')
      }
      
      #################################################
      
      mycov_PT<-matrix(mycorr,nrow=nsnp+n2sim,ncol=nsnp+n2sim)
      mycov_PT[1:nsnp,1:nsnp]<-mycov2
      m<-diag(n2sim)
      pheno_mat<-m[lower.tri(m)|upper.tri(m)]<-0.75
      mycov_PT[(nsnp+1):(nsnp+n2sim),(nsnp+1):(nsnp+n2sim)]<-pheno_mat
      
      #snpeff_pos<-sample(1:nsnp,nsnpeff,replace=F)
      
      for(j in (nsnp+1):(nsnp+n2sim))
      {
        if(!nsnpeff==0){
          if(cluster==0){snpeff_pos<-sample(1:nsnp,nsnpeff,replace=F)} else {snpeff_pos<-c((nsnp-(nsnpeff-1)):nsnp)}
          
          mycov_PT[snpeff_pos,j] <- gpcov#*ifelse(rbinom(1,1,c(0.5,0.5))==1,-1,1)
          mycov_PT[j,snpeff_pos] <- t(mycov_PT[snpeff_pos,j])
        }
      }
      
      mycov<-mycov_PT
      
      for(i in (nsnp+1):(nsnp+n2sim))
      {
        for(j in 1:nsnp)
        {
          mycov[i,j] <- ifelse(rbinom(1,1,prob=0.5)==1,mycov[i,j],(-1*mycov[i,j]))
          mycov[j,i] <- ifelse(rbinom(1,1,prob=0.5)==1,mycov[j,i],(-1*mycov[j,i]))
        }
      }
      
      mycov[upper.tri(mycov)] <- t(mycov)[upper.tri(mycov)]
      
      diag(mycov)<-rep(1,dim(mycov)[1])
      #################################################
      
      #then set diagonal values to 1 for all
      #diag(mycov)<-rep(1,n2sim)
      if(matrixcalc::is.positive.definite(mycov)==FALSE)
      {mycov<-Matrix::nearPD(mycov,keepDiag=TRUE)$mat}
      
      
      
      mymean<-rep(0,dim(mycov)[1])
      mydata=mvrnorm(n = ncases, mymean, mycov)
      
      mydata<-as.data.frame(mydata)
      
      colnames(mydata)[1:nsnp]<-snpnames
      colnames(mydata)[(nsnp+1):(nsnp+3)]<-c('NwdRepPheno','LangPheno','NeurodevPheno')
      
      
      
      #-------------------------------------------------------------------------
      #Convert gene scores to integers: 0-2 for autosomal
      
      firstcol<-1
      lastcol<-nsnp
      p<-c(0,0,0) #initialise a vector to hold p values for different N alleles
      for (i in 1:nsnp){
        p[1]<-(1-maf[i])^2
        p[2]<-2*(1-maf[i])*maf[i]
        p[3]<-maf[i]^2
        
        #now use p-values to assign z-score cutoffs that convert to 0,1,2 or 3 minor alleles
        temp<-mydata[,i]
        w<-which(temp<qnorm(p[1]))
        mydata[w,i]<-0
        w<-which(temp>qnorm(p[1]))
        mydata[w,i]<-1
        w<-which(temp>qnorm(p[2]+p[1]))
        mydata[w,i]<-2
        
      }
      
      myr<-cor(mydata) #to check you have desired correlation structure, View(myr)
      
      
      #----------------------------------------------------------------------

        
        #cor.data<-as.matrix(abs(myr))
        cor.data<-as.matrix(myr)
        cor.data[lower.tri(cor.data)] <- NA
        
        cor_tri_N <- as.data.frame(cor.data) %>% 
          mutate(Var1 = factor(row.names(.), levels=row.names(.))) %>% 
          gather(key = Var2, value = value, -Var1, na.rm = TRUE, factor_key = TRUE) 
        
        g2<-ggplot(data = cor_tri_N, aes(Var2, Var1, fill = value)) + 
          geom_tile()+scale_fill_gradient2(low = "blue", mid="white", high = "red")+theme_bw()+ theme(axis.text.x = element_text(angle = 90, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.title = element_blank(),legend.title=element_blank(),text=element_text(size=14),title=element_text(size=10))+geom_hline(yintercept = (nsnp+0.5),colour="grey")+geom_vline(xintercept = (nsnp+0.5),colour="grey")+ggtitle("SLIC LD pattern")
        
       
    g3<-ggarrange(g1,g2,common.legend = TRUE,legend="top")
    
print(g3)
    ggsave(paste0("/Users/paulthompson/Dropbox/project SCT analysis/GSCA validation/July2019/Corr_plots/combined_plot2.tiff"),dpi = 400,width=10,height=6,units="in")
    
}


combos <- expand.grid(nsnpeff=c(5), nsnp=c(20), ngenes=1, nsub=100,cluster=0, gpcorr=c(0.15),phenocorr=c(0.75))

for(cb in 1:length(combos[,1]))
{

  gsca_sim_multigene_plot(nsnp=combos[cb,2],nsnpeff=combos[cb,1],ngenes=combos[cb,3],nsub=combos[cb,4],ncases=10000,gpcorr=combos[cb,6],n2sim=3,cb=cb,cluster=combos[cb,5],phenCorr = combos[cb,7])
  
}


