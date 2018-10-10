# Written by DMW, 11/7/17
# Edited by JA, 5/2/18

args <- commandArgs(trailingOnly=TRUE)
  
# Check if all required arguments are supplied
if (length(args)<5) {
  stop("Missing 1 or more required arguments. Required: \n(1) <path to file with 'true' (observed) results>\n(2) <path to file with burden shift p-values from all annotation categories (rows) and permutations (columns)> \n(3) <path to table listing category names and category set membership ,as 0s and 1s in subsequent columns>\n(4) <number of permutations in results file (n columns - 1)>\n(5) <maximum p-value threshold for counting (2-sided binomial test)>\n(6) <name tag for output files>",call.=FALSE)
}

# Assign arguments to objects
burdenResFile <- args[1]
#
burdenShiftFile <- args[2]
#
catsFile <- args[3]
#
nPerms <- as.numeric(args[4])
#
maxP <- as.numeric(args[5])
if (maxP>1) {stop(paste("P-value threshold ",maxP," is greater than 1. Provide a p<1, e.g. 0.05."))}
#
if (length(args)==6) {
  outtag <- args[6]
} else {
  outtag <- "out"
}
#
rm(args)

# Load required packages
require(readr)
require(ggplot2)
# Function for counting categories
countCats <- function(pvals, pvalThresh) {
  nCase <- length(pvals[pvals>0 & pvals<=pvalThresh])
  nControl <- length(pvals[pvals<0 & abs(pvals)<=pvalThresh])
  return(c("case"=nCase,"control"=nControl))
}

# Load observed results file (true data)
burdenRes <- read.delim(burdenResFile, sep="\t", header=T, stringsAsFactors=F)

# Load burden shift results file (permuted data)
burdenShift <- read_delim(burdenShiftFile, delim="\t", col_names=T, col_types=paste("c",paste(rep("d",times=nPerms),collapse=""),sep=""))

# Load file noting category sets to process
catsets <- read.delim(catsFile, sep="\t", header=T, stringsAsFactors=F)

# Initialize results objects
permCaseTab <- permControlTab <- c()
pdf(file=paste("plotDistr_p",maxP,"_",outtag,".pdf",sep=""), height=3, width=4)

# Loop through categories specified in the catsets object
for (i in 1:ncol(catsets)) {
  
  # Define the name of the set, and the categories within it
  if (i==1) {
    setName <- "All"
    testCats <- catsets$Annotation_combo
  } else {
    setName <- colnames(catsets)[i]
    testCats <- catsets$Annotation_combo[catsets[,i]==1]
  }
  print(setName)
  print(length(testCats))
  #
  # Subset true results to only the categories in the set
  burdenResTrim <- data.frame(burdenRes[burdenRes$Annotation_combo %in% testCats,])
  #
  # Find case and control counts from this subset of categories
  nObsCase <- nrow(burdenResTrim[burdenResTrim$Binom_p<=maxP & burdenResTrim$Adjusted_relative_risk>1,])
  nObsControl <- nrow(burdenResTrim[burdenResTrim$Binom_p<=maxP & burdenResTrim$Adjusted_relative_risk<1,])
  #
  # Subset burden shift results to only the categories in the set
  burdenShiftTrim <- data.frame(burdenShift[burdenShift$Annotation_combo %in% testCats,])
  #
  # Count up the number of case and control categories passing maxP
  permCounts <- data.frame(t(apply(burdenShiftTrim[,2:ncol(burdenShiftTrim)], 2, countCats, pvalThresh=maxP)))
  #
  # Compare the observed counts to the permuted counts to calculate shift p-values
  nPermCase <- nrow(permCounts[permCounts$case>=nObsCase,])
  pCase <- nPermCase/nrow(permCounts)
  nPermControl <- nrow(permCounts[permCounts$control>=nObsControl,])
  pControl <- nPermControl/nrow(permCounts)
  #
  # Plot the results
  ggpermCounts <- data.frame("N_signif_tests"=c(permCounts$case,permCounts$control))
  myfontsize <- 10
  larger <- max(c(nObsCase,nObsControl))
  shiftDistPlot <- ggplot(ggpermCounts, aes(x=N_signif_tests)) +
    geom_density(size=0.25, fill="gray92") +
    geom_vline(mapping=aes(xintercept=nObsCase), size=0.25, colour="red", linetype=2) +
    geom_vline(mapping=aes(xintercept=nObsControl), size=0.25, colour="blue", linetype=2) +
    labs(x='Number of significant tests', y='Density', title=setName) +
    theme_bw() +
    theme(text=element_text(size=myfontsize), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.title=element_text(size=myfontsize, face="bold"), plot.title=element_text(hjust=0.5,size=myfontsize, face="bold"), axis.text=element_text(size=myfontsize), strip.text=element_text(size=myfontsize)) +
    annotate("text", x=nObsCase+(larger*0.02), y=Inf, vjust=4, hjust=0, label=pCase, colour="red") +
    annotate("text", x=nObsControl+(larger*0.02), y=Inf, vjust=2, hjust=0, label=pControl, colour="blue")
  print(shiftDistPlot)
  #
  # Add plot to figure document
  #
  # Add data to output object
  permCaseTab <- data.frame(cbind(permCaseTab,permCounts$case))
  permControlTab <- data.frame(cbind(permControlTab,permCounts$control))
  colnames(permControlTab)[ncol(permControlTab)] <- colnames(permCaseTab)[ncol(permCaseTab)] <- setName
  if (i==1) {
    obsTab <- data.frame("Category_set"=setName,"N_cats_case"=nObsCase,"P_case"=pCase,"N_cats_control"=nObsControl,"P_control"=pControl, stringsAsFactors=F)
  } else {
    obsTab <- rbind(obsTab,c("Category_set"=setName,"N_cats_case"=nObsCase,"P_case"=pCase,"N_cats_control"=nObsControl,"P_control"=pControl))
  }
}
dev.off()

# Write output data to file(s)
# File names
permCaseFile <- paste("nCatsCase_perm_p",maxP,"_",outtag,".txt",sep="")
permControlFile <- paste("nCatsControl_perm_p",maxP,"_",outtag,".txt",sep="")
obsFile <- paste("nCats_obs_p",maxP,"_",outtag,".txt",sep="")
#
# Write files
write.table(permCaseTab, file=permCaseFile, sep="\t", row.names=F, col.names=T, quote=F)
write.table(permControlTab, file=permControlFile, sep="\t", row.names=F, col.names=T, quote=F)
write.table(obsTab, file=obsFile, sep="\t", row.names=F, col.names=T, quote=F)
