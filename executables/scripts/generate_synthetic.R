# In this version:
# - generates data according to the region based Epiclomal with missing data
# - generates the clones independently and also following a phylogeny tree
# - missing and error are probabilities, for each position, make it erroneous or missing with that probability
# - the missing step is applied after the error step
# - only 2 states: methylated (1) and unmethylated (0)

# simulate synthetic data

#======================
# libraries
#======================
#.libPaths(c("/extscratch/shahlab/dalai/R/x86_64-pc-linux-gnu-library-centos5/3.2", "/clusterapp/software/linux-x86_64-centos5/R-3.2.3/lib64/R/library"))


library("argparse")
suppressMessages(library("Matrix"))
suppressMessages(library("MCMCpack")) # for generating values from a Dirichlet distribution

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add a help option

parser$add_argument("--read_size", type="character",default="1_0", help="mean and std dev for the number of CpGs per read when generating data considering a read-based approach, default is no read based approach, that is 1_0")

parser$add_argument("--num_samples", type="integer", default=1, help="Number of different samples, default is 1. More than 1 was not fully implemented and tested.")

parser$add_argument("--num_loci", type="integer", default=100, help="Number of loci")
parser$add_argument("--num_clones", type="integer", default=3, help="Number of clones")

parser$add_argument("--num_cells", type="character", default=10, help="Number of cells per sample, starting with sample 1, then sample 2, etc")

parser$add_argument("--clone_prevalence", type="character", default="0.2_0.5_0.3", help="Probability of clone prevalences for all samples, separated by _, they have to sum up to 1. For example 0.2_0.5_0.3 for 3 clones generates clone 1 with pb. 0.2, clone 2 with pb 0.5 and clone 3 with pb. 0.3")
parser$add_argument("--error_probability", type="character", default="0.01_0.01", help="Error probability for each methylation state, for example 0.01_0.01")
# error_probability="0.01_0.01" corresponds to P(Y=1|G=0)=0.01 and P(Y=0|G=1)=0.01, respectively
# If we assume more than 2 states it would be better if error_probability is written as matrix

### the missing probability could be in principle different among samples, keeping the same for now
parser$add_argument("--missing_probability", type="character", default=0.1, help="Missing data probability, it could be one probability for all cells or different ones")

### Genotype probabilities: p(G_krl = s) = mu_krs
parser$add_argument("--genotype_prob", type="character", default="dirichlet", help="dirichlet (genotype probabilities are draws from a dirichlet distribution) or 0.5_fixed (all genotype probabilities fixed to 0.5)")

### Generate them using a Dirichlet distribution
parser$add_argument("--dirichlet_param_genotype_prob", type="character", default="1_1", help="Dirichlet parameters to draw genotype probabilities for each region r and clone k, for example 1_1")
parser$add_argument("--percent_regions_dirichlet_param", type="character", default="0.9", help="The percentage of the regions that have methylation drawn from the above distribution. For example if I want 90 percent of the regions to be hypermethylated and 10 percent to be hypomethylated, I will set dirichlet_param_genotype_prob=99_1 and percent_regions_dirichlet_param=0.9")


parser$add_argument("--num_regions", type="double", default=5, help="Number of regions")

parser$add_argument("--region_size_type", type="character", default="multinomial_equal", help="uniform, multinomial_equal or multinomial_nonequal. Fixed generated from uniform (from 1 to nloci), multinomial_equal (with prob 1/nregions) or multinomial_nonequal (currently this has hard-coded probabilities)")

parser$add_argument("--output_dir", type="character", default="output", help="Entire or just the beginning of the output directory file ")
parser$add_argument("--given_dir_complete", type="integer", default=0, help="If this is 0, it creates a long output dir name with the input parameters, if it is 1, the output dir is output_dir ")

parser$add_argument("--plot_data", type="character", default=0, help="If this is 1, use the visualization software to plot the data")

parser$add_argument("--prop_add_var", type="character", default="0_0.5", help="Proportion (0-1) of non-flipped regions to have their methylation calls changed with prob = 0.5, default is to NO cell to cell variability")

parser$add_argument("--bulk_depth", type="integer", default=60, help="Number of cells that will be used to generate bulk methylation levels. If zero no bulk data will be saved.")


### arguments to generate following a phylogenetic tree
parser$add_argument("--phylogenetic_generation", type="integer", default=1, help="1 or 0. If this is 1, use phylogenetic tree to generate the clones, if this is 0, the clones are independent.")
parser$add_argument("--num_reg_flip", type="double", default=1, help="Number of regions to flip at each branch of the tree")

parser$add_argument("--prop_cpg_flip", type="double", default=1, help="Proportion of CpGs to be flipped inside a region")

# writes: 6 files
# 1: all input parameters
# 2. region coordinates
# 3: true matrix of genotypes
# 4: Z's (true clone memberships)
# 5: complete data ( all X's already with error)
# 6: incomplete data (data with some X's missing)

parser$add_argument("--seed", help="The variability seed. You can set the seed for reproducibility. If this is not set, it will generate a random one and save it in the output params.yaml file")
parser$add_argument("--verbose", type="integer", default=1, help="Set to 1 if you want details of the data generation")
parser$add_argument("--saveall", type="integer", default=1, help="Set to 1 if you want the save all the data (with errors but without missing observations)")

args <- parser$parse_args()

print(args)

#####################
# other parameters ##
#####################
# these parameters should be the same for all simulated data sets
# therefore, we need to set another seed here for them and keep that the same across all simulated data sets from a particular scenario

# Setting a seed to get fixed genotype probabilities and fixed region sizes
#set.seed(args$fixing_seed)

S = 2 # number of states, we will consider only two for now in all simulations

# disable the scientific notation so that the output dir has the name 100000, and not 1e+05
options(scipen=999)

# Make the output directory have all the input parameters, only if given_dir_complete is 0, see below
# adding the regions arguments in the output directory
output_dir <- paste0(args$output_dir ,"_readsize",args$read_size,"_loci", args$num_loci, "_clones", args$num_clones,
                     "_cells", args$num_cells, "_prev", args$clone_prevalence, "_errpb", args$error_probability,
                     "_mispb", args$missing_probability, "_gpb", args$genotype_prob, "_dirpar", args$dirichlet_param_genotype_prob,
                     "_nregs", args$num_regions, "_regionsize-", args$region_size_type, "_rnonequal-", args$region_nonequal)

# if seed is not provided, then length(args$seed) is 0
if (length(args$seed) > 0) {
  print ("Setting seed")
  set.seed(args$seed)
  output_dir <- paste0(output_dir, "_seed", args$seed)
} else {
  print ("Setting seed to random")
  set.seed(Sys.time())
}

#==============================================
### Generating the genotype probabilites
### when genotype_prob="dirichlet" we generate different mu_krs's for each k and r
### Important note: when R=1 and region_sizes[1] = num_loci we have the basic-GeMM as in that case we assumed P(G_km = s) = mu_k (probability does not depend on the locus)
if (args$genotype_prob=="0.5_fixed") {
  mu_array <- array(0.5,c(S,args$num_regions,args$num_clones))
}
if (args$genotype_prob=="dirichlet") {
  dirichlet_param <- as.double(unlist(strsplit(args$dirichlet_param_genotype_prob, split="_")))
  mu_array <- array(NA,c(S,args$num_regions,args$num_clones))
  # Note: if phylogeny is used, only the first row from mu_array is used. For no phylogeny, the entire mu_array is used.
  for (k in 1:args$num_clones) {
    # MA: The following command was executed before we introduced the argument  percent_regions_dirichlet_param
    # mu_array[,,k] <- t(round(rdirichlet(args$num_regions, alpha=c(dirichlet_param[1],dirichlet_param[2])),2))
    # MA: adding the percent_regions_dirichlet_param argument here
    # first generating num_regions random numbers
    random_methylation <- runif (args$num_regions, 0, 1)
    for (r in 1:args$num_regions) {
      if (random_methylation[r] <= as.double(args$percent_regions_dirichlet_param))  {
        mu_array[,r,k] <- t(round(rdirichlet(1, alpha=c(dirichlet_param[1],dirichlet_param[2])),2))
      } else {
        mu_array[,r,k] <- t(round(rdirichlet(1, alpha=c(dirichlet_param[2],dirichlet_param[1])),2))
      }
    }
  }
}

if(args$verbose) {
  print("mu_array")
  print(mu_array)
}


#==============================================
### Generating the regions
M <- args$num_loci

if (args$num_regions == 1){
  region_sizes <- M
  reg_coord<- as.matrix(t(c(1,region_sizes)))
  #print(reg_coord)
  #print("number of regions equal to 1")
}

if (args$num_regions > 1) {
  print("number of regions greater than 1")
  
  if (args$region_size_type=="multinomial_equal") {
    region_sizes <- rmultinom(1,size=M,p=rep(1/(args$num_regions),args$num_regions))
    reg_coord <- NULL
    for(r in 1:args$num_regions){
      if(r==1){
        reg_coord <- rbind(reg_coord,c(1 , region_sizes[r]))} else{
          reg_coord <- rbind(reg_coord,c(sum(region_sizes[1:(r-1)]) + 1 , sum(region_sizes[1:r]) ) ) }
    }
  } else if (args$region_size_type=="uniform"){
    breaks_unif <- c(0,sort(round(runif(args$num_regions-1,min=0,max=M))),M)
    while(sum(duplicated(breaks_unif))!=0){
      print("generating regions again")
      breaks_unif <- c(0,sort(round(runif(args$num_regions-1,min=0,max=M))),M)
    }
    region_sizes <- diff(breaks_unif)
    reg_coord <- NULL
    for(r in 1:args$num_regions){
      reg_coord <- rbind(reg_coord,c(breaks_unif[r]+1,breaks_unif[r+1]))
    }
  } else if(args$region_size_type=="multinomial_nonequal") {  # NOT USED
    #print("nonuniform")
    p_reg <- rep(c(0.1,0.2,0.5,0.1,0.1),args$num_regions/5)
    region_sizes <- rmultinom(1,size=M,p=p_reg) ### unbalanced sizes
    
    ### when R is big we will need a dirichlet if you want some control on which regions to be big or small
    ### Example for R=500
    #p_reg <- t(rdirichlet(1, alpha=c(rep(1,100),rep(4,100),rep(50,100),rep(8,100),rep(6,100))))
    
    reg_coord <- NULL
    for(r in 1:R){
      if(r==1){
        reg_coord <- rbind(reg_coord,c(1 , region_sizes[r]))
      }
      else{
        reg_coord <- rbind(reg_coord,c(sum(region_sizes[1:(r-1)]) + 1 , sum(region_sizes[1:r]) ) )
      }
    }
  } else {
    stop("region_size_type has to be multinomial_equal, uniform or multinomial_nonequal")
  }
}


# BEGIN FUNCTIONS
# ==========================
# function that flips all the positions_to_flip in vector from 0 to 1 or from 1 to 0
flip <- function (vector, positions_to_flip) {
  for (locus in positions_to_flip) {
    if (vector[locus] == 0) {
      vector [locus] <- 1
    } else if (vector[locus] == 1) {
      vector [locus] <- 0
    } else {
      print ("ERROR: value of vector is not 0 or 1!!")
    }
  }
  return (c(vector))
}

# parent_child_flip <- function(parent,region_coord,region_flip,prop_flip){
#   for(i in 1:length(region_flip)){
#   if(prop_flip == 1){
#     pos_flip <- region_coord[region_flip[i],1]:region_coord[region_flip[i],2]
#     child <- flip(vector = parent, positions_to_flip = pos_flip)
#   }
#   
#   if(prop_flip != 1){
#     pos_flip <- sort(sample(region_coord[region_flip[i],1]:region_coord[region_flip[i],2],size=rbinom(1,size=(region_coord[region_flip[i],2]-region_coord[region_flip[i],1]),prob=prop_flip)))
#     child <- flip(vector = parent, positions_to_flip = pos_flip)
#   }
#   
#   }
#   return(child)
# }


parent_child_flip <- function(parent,region_coord,region_flip,prop_flip){
  
  if(prop_flip == 1){
    pos_flip <- NULL
    for(i in 1:length(region_flip)){
      pos_flip <- c(pos_flip,region_coord[region_flip[i],1]:region_coord[region_flip[i],2])
    }
    child <- flip(vector = parent, positions_to_flip = pos_flip)
    
  }
  
  if(prop_flip != 1){
    pos_flip <- NULL
    for(i in 1:length(region_flip)){
      pos_flip <- c(pos_flip,sort(sample(region_coord[region_flip[i],1]:region_coord[region_flip[i],2],size=rbinom(1,size=(region_coord[region_flip[i],2]-region_coord[region_flip[i],1]),prob=prop_flip))))
    }
    child <- flip(vector = parent, positions_to_flip = pos_flip)
  }
  
  return(child)
}


#function that generates the observed data for each cell
xobs.function <- function(z,epsilon,genotype_matrix){
  g <- genotype_matrix[z,] ### that's Step 2 - assigning the cell a true genotype
  x <-rep(NA,length(g))
  
  ### one way of generating the x's - that's Step 3
  #ptm = proc.time()
  x1 <- g[g==1]
  x1[sample(1:sum(g==1),size=rbinom(n=1,size=sum(g==1),prob=epsilon[2]))] <- 0 ### epsilon[2] = P(Y=0|G=1)
  x[g==1] <- x1
  
  x0 <- g[g==0]
  x0[sample(1:sum(g==0),size=rbinom(n=1,size=sum(g==0),prob=epsilon[1]))] <- 1 ### epsilon[1] = P(Y=1|G=0)
  x[g==0] <- x0
  #(f1 <- proc.time() - ptm)
  
  ### another of way of doing it - sampling each obs from a bernoulli distribution
  ### it is a bit slower than above when considering M=1e+6
  ### If we consider more than 2 states it is easier to do this way sampling from a multinomial distribution
  #ptm = proc.time()
  #x1 <- rbinom(n=sum(g==1),size=1,prob=(1-epsilon[2,1])) ### prob of success --> P(Y=1|G=1) == 1-epsilon[2,1]
  #x[g==1] <- x1
  
  #x0 <- rbinom(n=sum(g==0),size=1,prob=epsilon[1,2]) ### prob of success --> P(Y=1|G=0) == epsilon[1,2]
  #x[g==0] <- x0
  #(f2 <- proc.time() - ptm)
  return(x)
}


genotype_reg <- function(r,region_sizes,mu_k){
  mu_kr <- as.matrix(mu_k)[,r]
  g_rk <- sample(c(0,1),size=region_sizes[r],prob=mu_kr,replace=TRUE)
  return(g_rk)
}

write_data_file <- function (data_matrix, output_file, index="cell_id") {
  # first write the header into the file
  # write the \n separately, otherwise it adds an extra tab
  cat(sapply(c(index,1:ncol(data_matrix)), toString), file=output_file, sep="\t")
  cat("\n", file=output_file, append=TRUE)
  rownames(data_matrix) <- c(1:nrow(data_matrix))
  write.table (data_matrix, gzfile(output_file) , sep="\t", col.names=FALSE, quote=FALSE, append=TRUE)
}

# END FUNCTIONS
# =============================



# if argument given_dir_complete is 1, then just use the given output_dir
if (args$given_dir_complete)
  output_dir<-args$output_dir

print(paste0("Output dir is: ", output_dir))

dir.create(output_dir, showWarnings = TRUE)
# add a data subdirectory
output_dir <- file.path(output_dir, "data")
dir.create(output_dir, showWarnings = TRUE)


### WE NEED TO SAVE mu_array or at least the seed that generates it
save(mu_array,file=file.path(output_dir, "mu_array.Rdata"))


### SAVE region coordinates
# ===================================================
reg_coord<-reg_coord-1      # I want the coordinates to start from 0
if (args$verbose)   {
  print ("Region coordinates")
  print (reg_coord)
}
regions_file <- file.path(output_dir, "regions_file.tsv.gz")
reg_file <- gzfile(regions_file)
# write_data_file(as.matrix(reg_coord)-1, reg_file, index="region_id")
# write_data_file(reg_coord, reg_file, index="region_id")
# adding start and end to the header, and I want the regions ids to start from 0 -- MAYBE WE SHOULD DO THIS FOR THE OTHER IDS??
tmp <- cbind(1:nrow(reg_coord)-1,reg_coord)
colnames(tmp) <- c("region_id","start","end")
write.table (tmp, reg_file, sep="\t", row.names=FALSE, quote=FALSE)
rm(tmp)

# Measuring the time with Rprof. For more details see http://stackoverflow.com/questions/6262203/measuring-function-execution-time-in-r
#Rprof ( tf <- file.path(output_dir, "log.log"),  memory.profiling = TRUE )

# Create the genotype matrix, that is, the vectors G_k's for k=1,...,K
# ==========================
print ("GENERATING EPIGENOTYPES")

R <- args$num_regions
K <- args$num_clones

########################################################
### Generating independent epigenotypes, no phylogeny ##
########################################################
# It uses the probabilities in mu_array
# It includes regions
# Remember note above: R=1 and region_size = num_loci corresponds to the basic-Epiclomal model

if (args$phylogenetic_generation == 0){
  genotype_matrix <- NULL
  for(k in 1:K){
    
    if(k==1){
      g_k <- as.vector(unlist(sapply(1:R,genotype_reg,region_sizes=region_sizes,mu_k=mu_array[,,k]))) ### the S element of mu_array[,r,k] is the probability of sucess
      genotype_matrix <- rbind(genotype_matrix,g_k)
      #if (args$verbose)   {
      #     print ("First clone of epigenotype matrix")
      #    print (genotype_matrix) }
      
    }else{
      g_k <- as.vector(unlist(sapply(1:R,genotype_reg,region_sizes=region_sizes,mu_k=mu_array[,,k])))
      while (sum(duplicated(rbind(genotype_matrix,g_k))==TRUE) != 0) { ### making sure there all genotype vectors are different from each other
        #if (args$verbose) {
        #  print ("Clone is already there, regenerate")
        #}
        #c <- c+1
        g_k <- as.vector(unlist(sapply(1:R,genotype_reg,region_sizes=region_sizes,mu_k=mu_array[,,k])))
      }
      genotype_matrix <- rbind(genotype_matrix,g_k)
    }
  }
  
  #final <- proc.time() - ptm
  if (args$verbose) {
    print ("Final epigenotype matrix: ")
    print (genotype_matrix)
  }
  
  
  # NOW write the matrix with each clone genotype into a file
  # ========================================
  geno_file <- file.path(output_dir, "true_clone_epigenotypes.tsv.gz")
  write_data_file(genotype_matrix, geno_file)
  
}

##############################################
### Generating epigenotypes using phylogeny ##
##############################################

if (args$phylogenetic_generation == 1){
  
  if (args$verbose) {
    print ("Generating epigenotypes using phylogeny")
  }
  
  
  if (args$num_regions == 1){
    
    if (args$verbose) {
      print ("Number of regions > 1, it does not work for phylogenies")
    }
    
    stop()
  }
  
  
  genotype_matrix <- NULL
  flipped_regions <- NULL
  r <- args$num_reg_flip ### number of regions to flip at each step of the phylogeny tree
  
  print("Changing one or more regions at a time")  
  
  for(k in 1:K){
    
    print(k)
    
    if(k==1){
      g_k <- as.vector(unlist(sapply(1:R,genotype_reg,region_sizes=region_sizes,mu_k=mu_array[,,k]))) ### the S element of mu_array[,r,k] is the probability of sucess
      genotype_matrix <- rbind(genotype_matrix,g_k)
      #if (args$verbose)   {
      #     print ("First clone of epigenotype matrix")
      #    print (genotype_matrix) }
      
    }
    
    if(R < r){
      print("Error: number of regions to flip greater than total number of regions")
      stop()
    }
    
    if(k==2){
      
      if( (R - r*(k-1)) < 0 ){
        print("Error: no more regions left to flip")
        stop()
      }
      
      region_flip <- sort(sample(1:R,r)) ### selection one region to flip per lineage
      print("regions to flip 2nd cluster")
      print(region_flip)
      
      flipped_regions <-  c(flipped_regions,region_flip)
      
      g_k <- parent_child_flip(parent=genotype_matrix[k-1,],region_coord=(reg_coord+1),region_flip = region_flip,prop_flip=args$prop_cpg_flip)
      
      genotype_matrix <- rbind(genotype_matrix,g_k)
      
    }
    
    if(k > 2){
      
      if( (R - r*(k-1)) < 0 ){
        print("Error: no more regions left to flip")
        stop()
      }
      
      ## picking a parent epigenotype
      
      parent_k <- sample(1:(k-1),size=1)
      print("parent epigenotype")
      print(parent_k)
      
      region_flip <- sort(sample(1:R,r))
      print("region to flip")
      print(region_flip)
      
      g_k <- parent_child_flip(parent=genotype_matrix[parent_k,],region_coord=(reg_coord+1),region_flip = region_flip,prop_flip=args$prop_cpg_flip)
      
      while (sum(duplicated(rbind(genotype_matrix,g_k))==TRUE) != 0) {
        ### making sure there all genotype vectors are different from each other
        
        print ("Clone is already there, regenerate")
        
        parent_k <- sample(1:(k-1),size=1)
        print("parent epigenotype")
        print(parent_k)
        
        region_flip <- sort(sample(1:R,1))
        
        print("region to flip")
        print(region_flip)
        
        g_k <- parent_child_flip(parent=genotype_matrix[parent_k,],region_coord=(reg_coord+1),region_flip = region_flip,prop_flip=args$prop_cpg_flip)
      }
      # here region_flip is surely accepted as unique
      flipped_regions <- c(flipped_regions, region_flip)
      
      genotype_matrix <- rbind(genotype_matrix,g_k)
    }
  }
  
  print ("All the flipped regions are")
  print (flipped_regions)
  
  # Write the flipped regions to a file
  write(flipped_regions, file=file.path(output_dir, "flipped_regions.tsv"))
  
  if (args$verbose) {
    print ("Final epigenotype matrix: ")
    print (genotype_matrix)
  }
  
  # print(genotype_matrix)
  
  # NOW write the matrix with each clone genotype into a file
  # ========================================
  geno_file <- file.path(output_dir, "true_clone_epigenotypes.tsv.gz")
  write_data_file(genotype_matrix, geno_file)
  
}

# Now generate the cells
# =============================
print ("GENERATING CELLS")

### Overall, generating the cell data involves 4 steps
### Step 1: generate the the clone membership of each cell, that is, the Z_i's
### Step 2: assign to each cell to its true genotype from the genotype_matrix according to the value of Z
### Step 3: add error to the true genotype generating the observed data for each cell
### Step 4: remove some observations from each cell to allow for missing data

###########
### Step 1: generating the clone membership of each cell, that is, the Z_i's
###########

### We assume that P(Z_i=k) = pi_k where pi_k is the clonal prevalence of clone k

clone_prev <- as.double(unlist(strsplit(args$clone_prevalence, split="_"))) ### this is the vector with the pi_k's

num_cells <- as.double(unlist(strsplit(args$num_cells, split="_")))

### this allows the number of cells to be different between samples
Z <- NULL
for(s in 1:args$num_samples){
  #print(s)
  Z <- c(Z,sample(1:K,size=num_cells[s],prob=clone_prev[((s*(K-1)+1) - (K-s)):(s*K)],replace = TRUE))
}

#print(Z)
#print(length(Z))

sample_id <- NULL
for (s in 1:args$num_samples){
  sample_id <- c(sample_id, rep(s,num_cells[s]))
}

cell_id <- NULL
for (s in 1:args$num_samples){
  cell_id <- c(cell_id, seq(from=1,to=num_cells[s],by=1))
}

# MA 23 Sep 2017: for now I am removing _sample_id from the final cell ID because epiclomal doesn't know how to deal with that yet. I will add it later when I add samples in epiclomal.
#cell_id_sample_id <- paste0(cell_id,"_",sample_id)
cell_id_sample_id <- cell_id

#print(cell_id_sample_id)

tmp <- cbind(cell_id_sample_id,Z)
colnames(tmp) <- c("cell_id","epigenotype_id")
clone_filename <- file.path(output_dir, "true_clone_membership.tsv.gz")
clone_file <- gzfile(clone_filename)
write.table (tmp, clone_file, sep="\t", row.names=FALSE, quote=FALSE)
rm(tmp)

##################
### Steps 2 and 3: generating the observed data for each cell (X_nm or X_nrl)
##################
### generating a matrix where each row is the true genotype of a cell
#data_true <- (t(sapply(Z,function(x){return(genotype_matrix[x,])}))) ### my genotype_matrix is K by M, so data_true is N by M
### we do not need to save data_true as it can be easily calculated using Z and genotype_matrix
### we can combine Steps 2 and 3 together and simply generate the observed data

### adding error to data_true
### this will depend on the epsilon_st's
### P(Y_nm = t | G_km = s) = epsilon_st
### the epsilon's are written in a matrix, for example,
### epsilon <- rbind(c(0.99,.01),c(0.02,0.98))
### the rows of epsilon should sum to one
### apply(epsilon,1,sum)

## If we assume that the rows of epsilon are the same (same error for all states in S)
## we can then simply flip the true genotype values with probability epsilon[1,2] without considering the true genotype states

### xobs.function is the function that generates the vector of observed genotypes for each cell
### epsilon should be a 2 by 2 matrix where the first row gives you P(Y=0|G=0) and P(Y=1|G=0) and the second row P(Y=0|G=1) and P(Y=1|G=1)
### so far we can only accommodate S=2 states

error_prob <- as.double(unlist(strsplit(args$error_probability, split="_")))
mean_read_size = as.double(unlist(strsplit(args$read_size, split="_")))[1]
sd_read_size = as.double(unlist(strsplit(args$read_size, split="_")))[2]

#Rprof(tmp_prof <- tempfile(),line.profiling=TRUE)

if (args$saveall) {
  cat(sapply(c("cell_id",1:ncol(genotype_matrix)), toString), file= file.path(output_dir, "data_complete.tsv"), sep="\t")
  cat("\n", file=file.path(output_dir, "data_complete.tsv"), append=TRUE)
}

cat(sapply(c("cell_id",1:ncol(genotype_matrix)), toString), file= file.path(output_dir, "data_incomplete.tsv"), sep="\t")
cat("\n", file=file.path(output_dir, "data_incomplete.tsv"), append=TRUE)

########################################
### adding cell to cell variability  ###
########################################

prop_add_var <- as.double(unlist(strsplit(args$prop_add_var, split="_")))

print(prop_add_var)

print(prop_add_var[1])

if( prop_add_var[1] != 0){
  
  print("adding cell to cell variability")
  regions_c <- 1:R
  number_non_flipped <- length(regions_c[-flipped_regions])
  #print(number_non_flipped)
  number_non_flipped_var <- max(1,(prop_add_var[1]*number_non_flipped))
  #print(number_non_flipped_var)
  regions_non_flipped <- sort(sample(regions_c[-flipped_regions],size=number_non_flipped_var))
  #print(regions_non_flipped)
  
}

for (n in 1:length(cell_id_sample_id)){
  
  cell_n <- xobs.function(z=Z[n],epsilon=error_prob,genotype_matrix=genotype_matrix)
  
  ##########################
  ## saving complete data ##
  ##########################
  
  if (args$saveall) {
    write.table(t(as.matrix(c(cell_id_sample_id[n],cell_n))), file.path(output_dir, "data_complete.tsv") , sep="\t",row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE)
  }
  
  ###########################
  ## including missingness ##
  ###########################
  
  if(grepl("_",x=args$missing_probability) == TRUE){
    unif_range <- as.double(unlist(strsplit(args$missing_probability, split="_")))
    miss_prob <- round(runif(1,min=unif_range[1],max=unif_range[2]),2)
  }else{
    miss_prob <- as.numeric(args$missing_probability)
  }
  
  ideal_obs <- rbinom(1,size=length(cell_n),prob=(1-miss_prob))
  
  obs_indices <- sample(1:length(cell_n),size=ideal_obs,replace=F)
  
  vect_data_new <- rep(NA,length(cell_n))
  
  total_obs <- 0
  i <- 0
  
  #  if (args$verbose) {
  #    print(ideal_obs)
  #  }
  
  while ( total_obs < ideal_obs) {
    
    read_length <- round(rnorm(1,mean=mean_read_size,sd=sd_read_size))
    
    i <- i+1
    
    if( (ideal_obs - total_obs) > read_length){
      
      mid <- round(read_length/2)
      
      if((obs_indices[i] - mid) >= 1){
        vect_data_new[(obs_indices[i] - mid):(obs_indices[i]+((read_length-mid)-1))] <- cell_n[(obs_indices[i] - mid):(obs_indices[i]+((read_length-mid)-1))]
      }else{
        vect_data_new[((obs_indices[i] - mid)-(obs_indices[i] - mid)+1):((obs_indices[i])+((read_length-mid)-1))] <- cell_n[((obs_indices[i] - mid)-(obs_indices[i] - mid)+1):((obs_indices[i])+((read_length-mid)-1))]
      }
      
      vect_data_new <- vect_data_new[1:length(cell_n)]
      
    }else{
      vect_data_new[(obs_indices[i]:(obs_indices[i]+ ((ideal_obs - total_obs)-1) ))] <- cell_n[(obs_indices[i]:(obs_indices[i]+  ((ideal_obs - total_obs)-1) ))]
      vect_data_new <- vect_data_new[1:length(cell_n)]
    }
    
    total_obs <- sum(!is.na(vect_data_new))
    
  }
  
  
  if( prop_add_var[1] != 0){
    
    for(c in regions_non_flipped){
      
      change <- rbinom(1,1,prob=prop_add_var[2])
      
      if( change == 1 ){
        
        coord_tmp <- reg_coord[c,] + 1
        
        vect_tmp <- vect_data_new[coord_tmp[1]:coord_tmp[2]]
        
        vect_tmp[!is.na(vect_tmp)] <- 0
        
        vect_tmp[vect_data_new[coord_tmp[1]:coord_tmp[2]] == 0] <- 1
        
        vect_data_new[coord_tmp[1]:coord_tmp[2]] <- vect_tmp
        
        rm(vect_tmp)
      }
      
    }
    
  }
  
  
  vect_data_new[is.na(vect_data_new)] <- ""
  
  
  write.table(t(as.matrix(c(cell_id_sample_id[n],vect_data_new))), file.path(output_dir, "data_incomplete.tsv") , sep="\t",row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE)
  
}

if (args$saveall) {
  system(paste0("gzip --force ", file.path(output_dir, "data_complete.tsv")))
}

system(paste0("gzip --force ", file.path(output_dir, "data_incomplete.tsv") ))

#print("SUMMARY Rprof for generating cells")
#Rprof()
#print(summaryRprof(tmp_prof,lines="both"))


if( args$bulk_depth != 0 ){
  
  print("GENERATING BULK DATA")
  
  ## Rprof(tmp_prof_bulk <- tempfile(),line.profiling=TRUE)
  
  bulk_data <- matrix(0,nrow=dim(genotype_matrix)[2],ncol=3)
  # MA 19 Oct 2017: reversing meth_reads and unmetth_read, putting unmeth first because that is how we do in epiclomal
  colnames(bulk_data) <- c("position","unmeth_reads","meth_reads")
  
  bulk_data[,1] <- 1:dim(genotype_matrix)[2]
  
  for(s in 1:args$num_samples){
    
    Z <- sample(1:K,size=args$bulk_depth,prob=clone_prev[((s*(K-1)+1) - (K-s)):(s*K)],replace = TRUE)
    
    for(n in 1:length(Z)){
      
      cell_n <- xobs.function(z=Z[n],epsilon=error_prob,genotype_matrix=genotype_matrix)
      
      bulk_data[,3] <- bulk_data[,3] + cell_n
      
      ######################################
      ## saving complete data per sample ###
      ######################################
      #write.table(t(as.matrix(cell_n)), file.path(output_dir, "bulk_cell_data_complete","_sample_",s,".tsv"),sep="\t",row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE)
      
      rm(cell_n)
      
    }
    
    bulk_data[,2] <- args$bulk_depth - bulk_data[,3]
    
    #print(reg_coord)
    
    #print(mu_array)
    
    if (args$verbose){
      print(head(bulk_data))
    }
    
    write.table(bulk_data, file.path(output_dir, paste0("bulk_data_sample_",s,".tsv")),sep="\t",row.names=FALSE, col.names=TRUE, quote=FALSE, append=FALSE)
    
    rm(Z)
    
    #system(paste0("gzip --force ", file.path(output_dir, "bulk_cell_data_complete","_sample_",s,".tsv")))
    
    system(paste0("gzip --force ", file.path(output_dir, paste0("bulk_data_sample_",s,".tsv"))))
    
  }
  
  
  
  ## print("SUMMARY Rprof for bulk")
  ## Rprof()
  ## print(summaryRprof(tmp_prof_bulk,lines="both"))
}


print("End of data generation")

if (args$plot_data == 0) {
  print("NOT PLOTTING GENERATED DATA")
}


if (args$plot_data == 1) {
  
  
  print("PLOTTING GENERATED DATA")
  
  suppressMessages(library("REpiclomal"))
  library(pheatmap)
  
  meth_file <- file.path(output_dir, "data_incomplete.tsv.gz")
  # visualization is a function in library REpiclomal
  visualization(outdir=output_dir,
                input_CpG_data_file=meth_file,
                input_regions_file=regions_file,
                input_CN_data_file=NULL,
                name="data",
                inferred_clusters_file=NULL,
                true_clusters_file=paste0(clone_filename),
                order_by_true=1,
                regions_to_plot=file.path(output_dir, "flipped_regions.tsv"))
  
}
