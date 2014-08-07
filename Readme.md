This is code to implement the MSA cross validation method outlined in 
Garvin et al. (in review)
Code was initially written by Garvin for his PhD thesis and substantially
altered by Patrick Barry as a side project. 
The code was written on a PC box, but checked on a Mac running
OSX. 

The code is written as three R functions that should be run in order.
The three functions are:

BaseMix_csv: This function produces csv files for mixture and baseline files.
The arguments that can be passed are:
x - the input file name
Prop - the proportion you want to use, default is 10
loci.diploid - the number of diploid loci
loci.haploid - the number of haploid loci

Bayes_format: this formats the .csv files for use in Bayes. The arguments that 
can be passed are:
x - the original input file used to create the .csv files
cores - the number of cores the computer has
npops - the number of populations in the baseline
loci.diploid - the number of diploid loci
loci.haploid - the number of haploid loci

Ctl_file: this creates a control file for Bayes. The arguments that it can pass are:
MCMC.num - number of MCMC samples to take
TI.SP - The thinning interval for the stock proportions
TI.BAF - The thinning interval for the baseline allele frequencies
TI.SAI - The thinning interval for the 
BAYES.options - The options for bayes.exe that need to be passed  an example is 'T T T T T T T'
Stock.ID - The stock ID block for Bayes. This is a tab delim text file that defines the stock
ID block with things like the dirichlet prior. 
Pops - The number of pops in the baseline
loci.diploid - the number of diploid loci
loci.haploid - the number of haploid loci

To construct the baseline, mixture, and control files the LoT.R file should be opened in
an R editor such as Rstudio and executed in sequence. To create the .csv files for the example file
'Baseline74.csv' use the command
BaseMix_csv('Baseline74.csv',loci.diploid=20,loci.haploid=1,print_pop.sum=T)
Then to create the baseline.bse and mixture.mix files use the command
Bayes_format('Baseline74.csv',cores=4,loci.haploid=1,loci.diploid=20,npops=74)
And finally to create the control file use the command 
Ctl_file(MCMC.num=1000,TI.SP=10,TI.BAF=10,TI.SAI=100,BAYES.options="T T F T F T F",Stock.ID='StockID.txt',Bayes.title='BayesTitle',Pops=74,loci.haploid=1,loci.diploid=20)

Because Bayes is extremely menu driven it is not possible to run 
it in batches from MSDOS command prompt. So the user must run each trial seperately.

In order to execute the code you will need to have your data 
formatted according to the example given:
The first column is the Individual
The second column is population of the individual
Then each column after is a single allele from a locus
diploid loci will occupy two columns and should come before haploid loci
Haploid loci should occupy two columns with the second column 
coded with 888 denoting that it is haploid.

To make the control file, a Stock.ID tab delim .txt file is needed
to specify the priors for the dirichlet and other stock
specific options for the program Bayes.
