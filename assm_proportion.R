library(optparse)

option_list <- list(
  make_option(c("-m","--mat"), action="store", type="character", default=NULL, help="Input the OTU table, with columns samples and rows OTUs, with last column is taxonomy"),
  make_option(c("-t","--tree"), action="store", type="character", default=NULL, help="Input otu trees. newick format"),
  make_option(c("-o", "--outdir"), action="store", default="./", type="character", help="The output dirctory, default is ./"),
  make_option(c("-r", "--reps"), action="store", default=999, type="double", help="Number of randomizations, default is 999"),
  make_option(c("-p", "--top"), action="store", default=1000, type="integer", help="Choosing top number OTUs for calculation, default is 1000")
)

opt <- parse_args(OptionParser(usage="%prog [options] file\n", option_list=option_list))

package_list <- c("picante","ggplot2","ecodist")
for(pk in package_list){
  if(!suppressWarnings(suppressMessages(require(pk, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    stop("WARNING: Please install ", pk, "!\n")
  }else{
    cat("YES:", pk, "was succesfully installed and loaded!\n")
  }
}
cat('####################START:', date(),'####################\n')

if(is.na(file.info(opt$outdir)$isdir)) dir.create(opt$outdir, recursive = TRUE)

# opt$mat <- '~/Desktop/assmbly_proportion/out.table.xls'
# opt$tree <- '/Users/wangpeng/Desktop/assmbly_proportion/rep_set.tre'

otu <-read.table(opt$mat, comment.char = '', 
                 sep = "\t", header = T, row.names = 1)

otu[,dim(otu)[2]] <- NULL

tran_outs <- as.data.frame(t(otu))


top <- ifelse(opt$top >= ncol(tran_outs), ncol(tran_outs), opt$top)

tran_outs <- tran_outs[ ,1:top]
all(tran_outs < 0.9)
phylo = read.tree(opt$tree);

match.phylo.otu = match.phylo.data(phylo, t(tran_outs));

# plot.phylo(match.phylo.otu$phy,typ="phylogram"); # a quick plot

beta.mntd.weighted = as.matrix(comdistnt(t(match.phylo.otu$data),cophenetic(match.phylo.otu$phy),abundance.weighted=T));

write.csv(beta.mntd.weighted, paste0(opt$outdir, '/betaMNTD_weighted.csv'),quote=F);

identical(colnames(match.phylo.otu$data),colnames(beta.mntd.weighted)); # just a check, should be TRUE
identical(colnames(match.phylo.otu$data),rownames(beta.mntd.weighted)); # just a check, should be TRUE

# calculate randomized betaMNTD

beta.reps = opt$reps; 

rand.weighted.bMNTD.comp = array(c(-999),dim=c(ncol(match.phylo.otu$data),ncol(match.phylo.otu$data),beta.reps));

for (rep in 1:beta.reps) {
  
  rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(match.phylo.otu$data),taxaShuffle(cophenetic(match.phylo.otu$phy)),abundance.weighted=T,exclude.conspecifics = F));
  
  print(c(date(),rep));
  
}

weighted.bNTI = matrix(c(NA),nrow=ncol(match.phylo.otu$data),ncol=ncol(match.phylo.otu$data));

for (columns in 1:(ncol(match.phylo.otu$data)-1)) {
  for (rows in (columns+1):ncol(match.phylo.otu$data)) {
    
    rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
    weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals);
    rm("rand.vals");
    
  };
};

rownames(weighted.bNTI) = colnames(match.phylo.otu$data);
colnames(weighted.bNTI) = colnames(match.phylo.otu$data);
write.csv(weighted.bNTI, paste0(opt$outdir, "/weighted_bNTI.csv"),quote=F);


pdf(paste0(opt$outdir, "/weighted_beta-NTI_Histogram.pdf"))
val.temp <- as.vector(as.matrix(weighted.bNTI))

val <- val.temp[!is.na(val.temp)]
ym = max(hist(val, breaks = 100, plot = F)$count)
xv <- sprintf("%.3f", mean(val, na.rm = T))
hist(val, breaks = 100, xlab = 'beta-NTI distribution', main = '')
abline(v = xv, lty=2, lwd=2,col="red")
text(x = xv, y=ym, labels = xv, col = 'blue')
dev.off()


#RC_bray
raup_crick_abundance = function(spXsite, plot_names_in_col1=TRUE, classic_metric=FALSE, split_ties=TRUE, reps=9999, set_all_species_equal=FALSE, as.distance.matrix=TRUE, report_similarity=FALSE){
  
  ##expects a species by site matrix for spXsite, with row names for plots, or optionally plots named in column 1.  By default calculates a modification of the Raup-Crick metric (standardizing the metric to range from -1 to 1 instead of 0 to 1). Specifying classic_metric=TRUE instead calculates the original Raup-Crick metric that ranges from 0 to 1. The option split_ties (defaults to TRUE) adds half of the number of null observations that are equal to the observed number of shared species to the calculation- this is highly recommended.  The argument report_similarity defaults to FALSE so the function reports a dissimilarity (which is appropriate as a measure of beta diversity).  Setting report_similarity=TRUE returns a measure of similarity, as Raup and Crick originally specified.  If ties are split (as we recommend) the dissimilarity (default) and similarity (set report_similarity=TRUE) calculations can be flipped by multiplying by -1 (for our modification, which ranges from -1 to 1) or by subtracting the metric from 1 (for the classic metric which ranges from 0 to 1). If ties are not split (and there are ties between the observed and expected shared number of species) this conversion will not work. The argument reps specifies the number of randomizations (a minimum of 999 is recommended- default is 9999).  set_all_species_equal weights all species equally in the null model instead of weighting species by frequency of occupancy.
  
  
  ##Note that the choice of how many plots (rows) to include has a real impact on the metric, as species and their occurrence frequencies across the set of plots is used to determine gamma and the frequency with which each species is drawn from the null model
  
  
  ##this section moves plot names in column 1 (if specified as being present) into the row names of the matrix and drops the column of names
  if(plot_names_in_col1){
    row.names(spXsite)<-spXsite[,1]
    spXsite<-spXsite[,-1]
  }
  
  
  ## count number of sites and total species richness across all plots (gamma)
  n_sites<-nrow(spXsite)
  gamma<-ncol(spXsite)
  
  ##build a site by site matrix for the results, with the names of the sites in the row and col names:
  results<-matrix(data=NA, nrow=n_sites, ncol=n_sites, dimnames=list(row.names(spXsite), row.names(spXsite)))
  
  ##make the spXsite matrix into a new, pres/abs. matrix:
  ceiling(spXsite/max(spXsite))->spXsite.inc
  
  ##create an occurrence vector- used to give more weight to widely distributed species in the null model:
  occur<-apply(spXsite.inc, MARGIN=2, FUN=sum)
  
  ##create an abundance vector- used to give more weight to abundant species in the second step of the null model:
  abundance<-apply(spXsite, MARGIN=2, FUN=sum)
  
  ##make_null:
  
  ##looping over each pairwise community combination:
  
  for(null.one in 1:(nrow(spXsite)-1)){
    for(null.two in (null.one+1):nrow(spXsite)){
      
      null_bray_curtis<-NULL
      for(i in 1:reps){
        
        ##two empty null communities of size gamma:
        com1<-rep(0,gamma)
        com2<-rep(0,gamma)
        
        ##add observed number of species to com1, weighting by species occurrence frequencies:
        com1[sample(1:gamma, sum(spXsite.inc[null.one,]), replace=FALSE, prob=occur)]<-1
        com1.samp.sp = sample(which(com1>0),(sum(spXsite[null.one,])-sum(com1)),replace=TRUE,prob=abundance[which(com1>0)]);
        com1.samp.sp = cbind(com1.samp.sp,1); # head(com1.samp.sp);
        com1.sp.counts = as.data.frame(tapply(com1.samp.sp[,2],com1.samp.sp[,1],FUN=sum)); colnames(com1.sp.counts) = 'counts'; # head(com1.sp.counts);
        com1.sp.counts$sp = as.numeric(rownames(com1.sp.counts)); # head(com1.sp.counts);
        com1[com1.sp.counts$sp] = com1[com1.sp.counts$sp] + com1.sp.counts$counts; # com1;
        #sum(com1) - sum(spXsite[null.one,]); ## this should be zero if everything work properly
        rm('com1.samp.sp','com1.sp.counts');
        
        ##same for com2:
        com2[sample(1:gamma, sum(spXsite.inc[null.two,]), replace=FALSE, prob=occur)]<-1
        com2.samp.sp = sample(which(com2>0),(sum(spXsite[null.two,])-sum(com2)),replace=TRUE,prob=abundance[which(com2>0)]);
        com2.samp.sp = cbind(com2.samp.sp,1); # head(com2.samp.sp);
        com2.sp.counts = as.data.frame(tapply(com2.samp.sp[,2],com2.samp.sp[,1],FUN=sum)); colnames(com2.sp.counts) = 'counts'; # head(com2.sp.counts);
        com2.sp.counts$sp = as.numeric(rownames(com2.sp.counts)); # head(com2.sp.counts);
        com2[com2.sp.counts$sp] = com2[com2.sp.counts$sp] + com2.sp.counts$counts; # com2;
        # sum(com2) - sum(spXsite[null.two,]); ## this should be zero if everything work properly
        rm('com2.samp.sp','com2.sp.counts');
        
        null.spXsite = rbind(com1,com2); # null.spXsite;
        
        ##calculate null bray curtis
        null_bray_curtis[i] = distance(null.spXsite,method='bray-curtis');
        
      }; # end reps loop
      
      ## empirically observed bray curtis
      obs.bray = distance(spXsite[c(null.one,null.two),],method='bray-curtis');
      
      ##how many null observations is the observed value tied with?
      num_exact_matching_in_null = sum(null_bray_curtis==obs.bray);
      
      ##how many null values are smaller than the observed *dissimilarity*?
      num_less_than_in_null = sum(null_bray_curtis<obs.bray);
      
      rc = (num_less_than_in_null )/reps; # rc;
      
      if(split_ties){
        
        rc = ((num_less_than_in_null +(num_exact_matching_in_null)/2)/reps)
      };
      
      
      if(!classic_metric){
        
        ##our modification of raup crick standardizes the metric to range from -1 to 1 instead of 0 to 1
        
        rc = (rc-.5)*2
      };
      
      results[null.two,null.one] = round(rc,digits=2); ##store the metric in the results matrix
      
      print(c(null.one,null.two,date()));
      
    }; ## end null.two loop
    
  }; ## end null.one loop
  
  if(as.distance.matrix){ ## return as distance matrix if so desired
    results<-as.dist(results)
  }
  
  return(results)
  
}

RC.bray <- raup_crick_abundance(tran_outs,reps=opt$reps,plot_names_in_col1=F)
RC.mat <- as.matrix(RC.bray)
write.csv(RC.mat,paste0(opt$outdir, '/RC_bray.csv'),quote=F)

total_pairs <- dim(weighted.bNTI)[1] * dim(weighted.bNTI)[2] / 2

selection <- sum(weighted.bNTI > 2, na.rm = T) / total_pairs

dispersal_limitation <- sum(RC.mat > 0.95 & weighted.bNTI < 2, na.rm = T) / total_pairs

homogenizing_dispersal <- sum(RC.mat < -0.95 & weighted.bNTI < 2, na.rm = T) / total_pairs

drift <- sum(abs(RC.mat) < 0.95 & weighted.bNTI < 2, na.rm = T) / total_pairs


cat("selection: ", selection, "\n")
cat("dispersal_limitation: ", dispersal_limitation, "\n")
cat("homogenizing_dispersal: ", homogenizing_dispersal, "\n")
cat("drift: ", drift, "\n")

res <- data.frame(Type = c("selection", "disp._lim.", "homo.disp.", "drift"),
                  proportion = c(selection, dispersal_limitation, homogenizing_dispersal, drift))

write.csv(res, file = paste0(opt$outdir, '/commu_assm_prop.csv'))
lb <- paste0(res$Type, " (", sprintf("%.2f", res$proportion*100), '%)')
pdf(paste0(opt$outdir, "/pie_chart.pdf"))
pie(res$proportion, labels = lb)
dev.off()


