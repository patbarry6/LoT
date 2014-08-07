### This is code for Garvin et al. (in review) Cross validation paper
### describing a leave out ten procedure for mixed stock analysis
### The code was written by Garvin for his PhD thesis and substantially
### altered by Patrick Barry as a side project. 
### The code was written for PCs but should be able to run on 
### both linux and mac OS as there are not any packages or system commands
### implemented. 
### There are three major pieces to the code:
### Piece one makes many baseline.bse files
### Piece two makes many mixture files
### Piece three makes a control file to run with the program BAYES
### Because Bayes is extremely menu driven it is not possible to run 
### it in batches from MSDOS command prompt. 

### In order to execute the code you will need to have your data 
### formatted according to the example given:
### The first column is the Individual #
### The second column is population of the individual
### Then each column after is a single allele from a locus
### diploid loci will occupy two columns and should come before haploid loci
### Haploid loci should occupy two columns with the second column 
### coded with 888 denoting that it is haploid.


### To make the control file, a Stock.ID tab delim .txt file is needed
### to specify the priors for the dirichlet and other stock
### specific options for the program Bayes.

BaseMix_csv<-function(x,Prop=10,loci.diploid,loci.haploid,print_pop.sum = TRUE)
{
  #This is a function to produce baseline and mixture files
  #from an input file. These files can then be formatted
  #for analysis with Bayes
  
  #read in data file
  Y<-read.csv(x,na.string="?")
  
  #Recode all the 888 for haploid loci with na
  Y[Y==888]<-NA
  
  #How many populations do we have?
  Pops <- length(unique(Y$Population))
  #How many individuals per pop?
  Pop.sum<-matrix(data=NA,ncol=4,nrow=Pops)
  Pop.sum[,1]<-seq(from=1, to = Pops,by=1)
  for (q in 1:Pops){
    Pop.sum[q,2]<-nrow(Y[Y$Population==q,,])
  }
  
  #We can add to our vector how many individuals to be taken that are 1/10th the population
  #and the leftovers to be added back
  for(w in 1:Pops){
    Pop.sum[w,3]<-Pop.sum[w,2]%/%10 
    Pop.sum[w,4]<- Pop.sum[w,2]%%10
  }
  colnames(Pop.sum)<-c("Pop","n","N","LO") 
  
  if (print_pop.sum == TRUE) {write.table(Pop.sum,quote=FALSE,file="Pop.sum.txt")
                              cat("Population Summary Printed to Pop.sum\n", sep="")
  } 
  
  ###How many unique alleles for each locus?
  Allele.sum<- matrix(data=NA,nrow=loci.diploid+loci.haploid,ncol=2) #matrix to hold information to reference later
  colnames(Allele.sum)<- c("Locus", "Max.Num")
  
  for (a.s in 1:nrow(Allele.sum)){ # fill in loci names
    Allele.sum[a.s,1]<- colnames(Y[a.s*2+1]) 
    all.alleles<-list((as.vector(na.omit(suppressWarnings(as.numeric(rownames(table(Y[,(a.s*2+1):(a.s*2+2)]))))))),(as.vector(na.omit(suppressWarnings(as.numeric(colnames(table(Y[,(a.s*2+1):(a.s*2+2)]))))))))
    length(all.alleles)
    
    u.all<-(unique(Y[,a.s*2+1]))
    u.all<-append(u.all,unique(Y[,a.s*2+2]))
    u.all<-unique(u.all)
    u.all<-sort(u.all)
    Allele.sum[a.s,2] <- length(u.all)
    assign(paste("Locus",a.s,sep=""),u.all)# make a vector for each with a list of the alleles found
  }
  
  write.table(Allele.sum,file="LociSummary.csv",quote=F,sep=",", row.names=F,col.names=F)
 
### Create indexing in the Y file
  Y$Index<-NA # Create index vector 
  index<- vector()
  for (v in 1:Pops){
    index[(sum(Pop.sum[Pop.sum[,1] < v,2])+1):(as.numeric((sum(Pop.sum[Pop.sum[,1] < v,2])+Pop.sum[v,3]*10)))]<- rep(1:10,rep((Pop.sum[v,3]),10))
    
    if (as.logical(Pop.sum[v,4]>0)){
      index[(as.numeric((sum(Pop.sum[Pop.sum[,1] < v,2])+Pop.sum[v,3]*10))+1):(as.numeric((sum(Pop.sum[Pop.sum[,1] < v,2])+Pop.sum[v,3]*10)+Pop.sum[v,4]))]<- rep(11,Pop.sum[v,4])
    }
  }
  Y$Index<- index
  
  
  #add Index shuffle values to matrix
  Y$Index.shuff<-NA
  Index.shuff.temp<- vector(mode="numeric", length=sum(Pop.sum[,2]))
  #Shuffle the index within each population
  for (u in 1:Pops){
    P.temp<- Y[Y$Population==u,]### This gives me the data I want to use
    Index.shuff<-P.temp$Index
    Index.shuff<- sample(Index.shuff,length(Index.shuff),replace=FALSE)
    #Write the shuffled values back to Y
    Index.shuff.temp[(sum(Pop.sum[Pop.sum[,1] < u,2])+1):(as.numeric((sum(Pop.sum[Pop.sum[,1] < u,2])+Pop.sum[u,2])))]<- Index.shuff
  }
  Y$Index.shuff<- Index.shuff.temp
  
  
  #Recombine all 101 x 10 mixtures to create 10 mixed stocks:      
  for (i in 1:10){
    assign(paste('Mix',i,sep=""),Y[Y$Index==i,])
    assign(paste('Baseline',i,sep=""),Y[Y$Index!=i,])
  }
  
  #Clean up file by getting rid of shuffling col. 
  for (i in 1:10){
    assign(paste('Baseline',i,sep=""),eval(parse(text=paste('Baseline',i,'[,1:(ncol(',(paste('Baseline',i,sep="")),')-2)]',sep=""))))
  }
  
  #Export the files:
  for (i in 1:10){
    write.table(eval(parse(text=paste('Mix',i,sep=""))),quote=FALSE,file=paste('Mix',i,'.csv',sep=""),row.names = FALSE, sep=",")
  }
  
  for (i in 1:10){
    write.table(eval(parse(text=paste("Baseline",i,sep=""))),quote=FALSE,file=paste("Baseline",i,".csv",sep=""),row.names = FALSE, sep=",")
  }
  cat('Baseline and Mixture files written to working directory')
 
}

############################################################################
#RtoBAYESBASE
Bayes_format<-function(x,cores=2,loci.haploid,loci.diploid,npops){
  # x is the original baseline to test
  # cores = number of cores on the computer
  # L is number of loci
  # requires that .csv files produced by BaseMix_csv() 
  # be located in the working directory
  
  #RtoBAYESBASE
  ####Read baseline file
  library(doSNOW)
  library(foreach)
  cl<-makeCluster(cores)
  registerDoSNOW(cl)
  
  Y<-read.csv(x,na.string="?")
  L<-loci.haploid + loci.diploid
  Loci.names<-paste("Locus",seq(from=1,to=L,by=1),sep="")
  
  
  ###How many unique alleles for each locus?
  Allele.sum<- matrix(data=NA,nrow=L,ncol=2) #matrix to hold information to reference later
  colnames(Allele.sum)<- c("Locus", "Max.Num")
  
  for (a.s in 1:nrow(Allele.sum)){ # fill in loci names
    Allele.sum[a.s,1]<- colnames(Y[a.s*2+1]) 
    all.alleles<-list((as.vector(na.omit(suppressWarnings(as.numeric(rownames(table(Y[,(a.s*2+1):(a.s*2+2)]))))))),(as.vector(na.omit(suppressWarnings(as.numeric(colnames(table(Y[,(a.s*2+1):(a.s*2+2)]))))))))
    length(all.alleles)
    
    u.all<-(unique(Y[,a.s*2+1]))
    u.all<-append(u.all,unique(Y[,a.s*2+2]))
    u.all<-unique(u.all)
    u.all<-sort(u.all)
    Allele.sum[a.s,2] <- length(u.all)
    assign(paste("Locus",a.s,sep=""),u.all)# make a vector for each with a list of the alleles found
  }
  
  cat('Please be patient, this may take a moment\n')
  
  foreach (o = 1:10,.export=Loci.names) %dopar% {
    Base.temp<-read.csv(paste("Baseline",o,".csv", sep="")) 
    #Calculate the number of populations
    P<-unique(Base.temp$Population)
    
    C<-ncol(Base.temp)
    
    # Create a matrix to index the alleles at each locus
    maxalleles= as.numeric(max(Allele.sum[,2]))
    Base <- matrix(data=NA,ncol=maxalleles+3, nrow=L*npops)#### rows here probably have to extend much further!
    colnames(Base)<- c("Pop","Locus","N",rep(1:max(Allele.sum[,2])))
    Base[,1]<- rep(1:npops,rep(L,npops))
    Base[,2]<-rep(1:L,npops)
    
    for(k in 1:npops){
      Pop.k <- Base.temp[Base.temp$Population==k,]
      
      for(i in 1:L){
        allele.temp <- eval(parse(text=paste("Locus",i,sep="")))# which alleles are we dealing with  
        ##add up all alleles
        a.vec<-vector(mode="numeric",length=length(allele.temp))
        for (a in 1:length(allele.temp)){
          if (i<= loci.diploid){
            a.vec[a]<-sum(Pop.k[,((2*i+1):(2*i+2))]== eval(parse(text=paste("Locus",i,"[a]",sep=""))),na.rm=T)
          } else {
            a.vec[a]<-sum(Pop.k[,(2*i+1)]== eval(parse(text=paste("Locus",i,"[a]",sep=""))),na.rm=T)
          }#if
          if (i<= loci.diploid){
            Base[(k-1)*L+i,4:(3+length(allele.temp))]<-a.vec # need to be clever about the 1 in the base
            Base[(k-1)*L+i,3]<- sum(Base[(k-1)*L+i,4:(3+length(allele.temp))],na.rm=T)
          } else {
            Base[(k-1)*L+i,4:(3+length(allele.temp))]<-a.vec # need to be clever about the 1 in the base
            Base[(k-1)*L+i,3]<- sum(Base[(k-1)*L+i,4:(3+length(allele.temp))],na.rm=T)*2
          }#if
        }#over allele
      }#over loci
    }#over pops
    
    #write the Base to a file
    write.table(Base,quote=FALSE,file=paste("Baseline",o,".txt",sep=""),row.names = FALSE,col.names=FALSE, sep="\t", na="")
  }
  ##### format baseline #######################################################
  require(stringr)
  foreach (o = 1:10,.packages=c('stringr')) %dopar% { 
    file.create(paste("Baseline",o,".bse",sep=""))
    line.temp<-readLines(paste("Baseline",o,".txt",sep=""))
    line.temp<-str_trim(line.temp)
    line.temp<-str_split(line.temp,"\t")
    max.ln<-npops*L
    for (l.n in 1:max.ln){
      line.temp2<-line.temp[[l.n]]
      line.temp2<-gsub(",","",paste(formatC(line.temp2, width=6,flag=" "),collapse=","))
      lapply(line.temp2, write, paste("Baseline",o,".bse",sep=""), append=TRUE)
    }}
  
  cat('Baseline.bse files are now saved to the working directory\n')
  #####################################################################################
  #### Lets format all the Mixture files now
  foreach (m = 1:10,.export=Loci.names) %dopar% {
    mix.temp<- read.csv(paste("Mix",m,".csv", sep=""))
    
    mix <- matrix(data=NA, ncol=sum(as.numeric(Allele.sum[,2]))+(L-1), nrow=nrow(mix.temp))
    
    #count alleles at each locus
    for (m.l in 1:L){
      allele.temp2 <- eval(parse(text=paste("Locus",m.l,sep="")))
      
      mix.a.mtrx<-matrix(data=NA,nrow=nrow(mix.temp),ncol=length(allele.temp2))
      
      if (m.l<= loci.diploid){ # count over locus
        for (ind.m in 1:nrow(mix.temp)){ # count for each individual in mix
          for (a.m in 1:length(allele.temp2)){ # count over all alleles at a locus
            mix.a.mtrx[ind.m,a.m]<-sum(mix.temp[ind.m, ((m.l*2+1):(m.l*2+2))]==allele.temp2[a.m],na.rm=T)
          }}
        if (m.l<2){# this section puts the results in the right place
          mix[,1:length(allele.temp2)]<- mix.a.mtrx
        } else {
          mix[,((sum(as.numeric(Allele.sum[((m.l-1):1),2]))+1)+m.l-1):(((sum(as.numeric(Allele.sum[((m.l-1):1),2])))+(length(allele.temp2-1)))+m.l-1)]<-mix.a.mtrx 
        }# if it is the first locus it needs to be placed in the right spot.
      } else { # if it is haploid
        for (ind.m in 1:nrow(mix.temp)){
          for (a.m in 1:length(allele.temp2)){
            mix.a.mtrx[ind.m,a.m]<-sum(mix.temp[ind.m, (m.l*2+1)]==allele.temp2[a.m],na.rm=T)
          }}
        mix[,((sum(as.numeric(Allele.sum[((m.l-1):1),2])))+m.l):(((sum(as.numeric(Allele.sum[((m.l-1):1),2])))+(length(allele.temp2-1)))+m.l-1)]<-mix.a.mtrx  
        gsub("NA"," ",mix) 
      }
    }
    # write the mixture to a file
    write.table(mix,quote=FALSE,na=" ",file=paste("Mixture",m,".mix", sep=""),row.names = FALSE, sep="", col.names=FALSE)
    
  }
  cat ('Mixture files are now saved to working directory. Have fun running Bayes.exe')
}
############################################################################
Ctl_file<-function(MCMC.num=1000,TI.SP=10,TI.BAF=10,TI.SAI=100,BAYES.options="T T F T F T F",
                   Stock.ID='StockID.txt',Bayes.title='BayesTitle',Pops=74,
                   loci.haploid,loci.diploid){ 
# Let's make a control file
  # we have 10 mixtures and 10 baselines!
  
  #User input for making the control file
  #MCMC.num<-How many MCMC samples?
  #TI.SP<-thinning interval for stock proportion
  #TI.BAF<-thinning interval baseline allele freq
  #TI.SAI<-thinning interval stock assignment for individual.
  #BAYES.options<-options for output of program SEE bayes manual ex."T T F T F T F"
  #Stock.ID<- StockID block in the control file. No way to make this generic.
  #Bayes.title<- Name for the analysis
  Allele.sum<-read.csv('LociSummary.csv',header=F)
  write.table(Bayes.title,quote=FALSE,file=paste(Bayes.title,'.ctl',sep=""),row.names = FALSE, sep="", col.names=FALSE) #write name to file
  write.table('Baseline1.bse',quote=FALSE,file=paste(Bayes.title,'.ctl',sep=""),append=TRUE,row.names = FALSE, sep="", col.names=FALSE)#baseline file
  write.table('Mixture1.mix',quote=FALSE,file=paste(Bayes.title,'.ctl',sep=""),append=TRUE,row.names = FALSE, sep="", col.names=FALSE) #mixture file
  #output names
  Out.names<-c('Summary.SUM','MCMC_StockProp.BOT','MCMC_AlleleFreq.FRQ','Binary_Restart.B01','Stock_IndAssign.CLS','MCMC_SGRP.RGN')
  write.table(Out.names,quote=FALSE,file=paste(Bayes.title,'.ctl',sep=""),append=TRUE,row.names = FALSE, sep="", col.names=FALSE) #output names
  #MCMC
  write.table(format(MCMC.num, scientific=F),quote=FALSE,file=paste(Bayes.title,'.ctl',sep=""),append=TRUE,row.names = FALSE, sep="", col.names=FALSE) #MCMC
  # number of stocks - get from input file. 
  write.table(Pops,quote=FALSE,file=paste(Bayes.title,'.ctl',sep=""),append=TRUE,row.names = FALSE, sep="", col.names=FALSE) #stocks 
  # number of characters 
  Characters<-loci.diploid+loci.haploid
  write.table(Characters,quote=FALSE,file=paste(Bayes.title,'.ctl',sep=""),append=TRUE,row.names = FALSE, sep="", col.names=FALSE) #loci 
  # 3 random seeds
  Seeds<-sample(1:200000,3)
  write.table(Seeds,quote=FALSE,file=paste(Bayes.title,'.ctl',sep=""),append=TRUE,row.names = FALSE, sep="", col.names=FALSE) #Seeds 
  #Thinning intervals
  write.table(c(TI.SP,TI.BAF,TI.SAI),quote=FALSE,file=paste(Bayes.title,'.ctl',sep=""),append=TRUE,row.names = FALSE, sep="", col.names=FALSE) #Thinning intervals 
  # Fortran Format for the mixture file
  Fortran.Mix<-paste('(',(paste(Allele.sum[,2],'I1,','1X',sep="",collapse=',')),')',sep="")
  write.table(Fortran.Mix,quote=FALSE,file=paste(Bayes.title,'.ctl',sep=""),append=TRUE,row.names = FALSE, sep="", col.names=FALSE) #Fortran mixture
  # Fortran Format of Baseline
  # column numbers in baseline I6)
  FtF_bse<-paste('(',(as.numeric(max(Allele.sum[,2]))+3),'I','6)',sep="")
  write.table(FtF_bse,quote=FALSE,file=paste(Bayes.title,'.ctl',sep=""),append=TRUE,row.names = FALSE, sep="", col.names=FALSE)
  
  #Options for output from BAYES
  write.table(BAYES.options,quote=FALSE,file=paste(Bayes.title,'.ctl',sep=""),append=TRUE,row.names = FALSE, sep="", col.names=FALSE) #Bayes Options
  
  # character descriptions
  Ch.Des<- matrix(data=NA,ncol=4,nrow=Characters)
  Ch.Des[,1]<-seq(from=1,to=Characters,by=1)
  Ch.Des[,2]<-Allele.sum[,2]
  Ch.Des[(1:loci.diploid),3]<-rep('T',time=loci.diploid)
  Ch.Des[((loci.diploid+1):(loci.diploid+loci.haploid)),3]<-rep('F',time=loci.haploid) 
  Ch.Des[,4]<-Allele.sum[,1]
  Ch.D<-paste(formatC(Ch.Des[,1], width = 3, flag = "+"),formatC(Ch.Des[,2], width = 3, flag = "+"),(paste(formatC(Ch.Des[,3], width = 3, flag = "+"),formatC(Ch.Des[,4], width = 3, flag = "-"),sep=" ")),sep="")
  write.table(Ch.D,quote=FALSE,file=paste(Bayes.title,'.ctl',sep=""),append=TRUE,row.names = FALSE, sep="", col.names=FALSE) #Bayes Options
  
  #Stock ID block
  # would this be easiest to just make the user make a stock id block?
  # read in the StockID file
  Stk.ID <- readLines(con=Stock.ID)
  #split into strings
  Stk.ID<- str_split(Stk.ID,pattern='\t')
  #what is the max length of the population name?
  C5<-max(nchar(sapply(Stk.ID, "[", 4)))
  Stk.ID.w<-paste(formatC(sapply(Stk.ID, "[", 1), width = 3, flag = "+"),formatC(sapply(Stk.ID, "[", 2), width = 3, flag = "+"),formatC(gsub( "^0+" , "" , sapply(Stk.ID, "[", 3) ), width = 8, flag = "+"),formatC(sapply(Stk.ID, "[", 4), width = C5+1, flag = "+"),formatC(gsub( "^0+" , "" , sapply(Stk.ID, "[", 5) ), width = 8, flag = "+"),sep="")
  write.table(Stk.ID.w,quote=FALSE,file=paste(Bayes.title,'.ctl',sep=""),append=TRUE,row.names = FALSE, sep="", col.names=FALSE) #Bayes Options
  cat (paste((paste(Bayes.title,'.ctl',sep="")),'is saved in working directory',sep=" "))
}