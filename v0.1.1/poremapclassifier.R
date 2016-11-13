# =========================================================================== #
# poremapclassifier.R                                                         #
#                                                                             #
# Classify a test set of Oxford Nanopore read alignments by whether their     #
# frequency of insertions, length of consecutive correct bases and length of  #
# consecutive matching bases is more similar to alignments of "minion" or     #
# "rand" reads from a known control reference DNA for which we have an        #
# identical reference.                                                        #
#                                                                             #
# Usage:                                                                      #
#                                                                             #
#     Rscript poremapclassifier.R \                                           #
#         runid \                                                             #
#         training_alignstats_path \                                          #
#         testing_alignstats_path \                                           #
#         outdir \                                                            #
#         outprefix                                                           #
#                                                                             #
# Input:                                                                      #
#                                                                             #
# - The training and testing data are the alignstats.txt files output by the  #
#   porealignstats.py program when the isnonrandomaln status is not known for #
#   each alignment, and set to default values.                                #
#                                                                             #
# - The fields currently used to build the model are:                         #
#      datatype : Data type : minion, random, mixed, unknown                  #
#    freqIperbp : Number of insertion segments per aligned base               #
#   meanrunlenC : Mean length of same-as-ref segments in alignment            #
#    freqMperbp : Number of match segments per aligned base                   #
#          pctC : Pct bases with same-as-ref call in read                     #
#          pctI : Pct bases in insertions in read                             #
#                                                                             #
# Output:                                                                     #
#                                                                             #
#   outdir/outprefix_alignclass: classification of each target read   #
#   outdir/outprefix_model.txt : model used for the classificiation           #
#   outdir/outprefix_PLOT.png  : diff bw real & random control reads          #
#                                                                             #
# Camilla Ip                                                                  #
# camilla.ip@well.ox.ac.uk                                                    #
# March 2015                                                                  #
# =========================================================================== #

# =========================================================================== #
# Constants
# =========================================================================== #

probability_threshold <- 0.99

# =========================================================================== #
# Command-line arguments
# =========================================================================== #

args <- commandArgs(trailingOnly = TRUE)
print(sprintf("Length(args)=%s", length(args)))
if (length(args) < 5) {
    txt <- sprintf("Usage: Rscript poremapclassifier.R runid training_initstats_path testing_initstats_path outdir outprefix")
    print(txt)
    quit(status=1)
}
runid <- args[1]
training_path <- args[2]
testing_path <- args[3]
outdir <- args[4]
outprefix <- args[5]

# =========================================================================== #
# Output file paths                                                           #
# =========================================================================== #

training_relvars_png <- sprintf("%s/%s_training_relativevalues.png", outdir, outprefix)
training_pccont_png <- sprintf("%s/%s_training_pccont.png", outdir, outprefix)
training_pcascores_png <- sprintf("%s/%s_training_pcascores.png", outdir, outprefix)
testing_relvars_png <- sprintf("%s/%s_test_relativevalues", outdir, outprefix)
testing_pccont_png <- sprintf("%s/%s_test_pccont.png", outdir, outprefix)
result_outpath <- sprintf("%s/%s_alignclass.txt", outdir, outprefix)

# =========================================================================== #
# Read in all the training and test data and sort by datatype                 #
# =========================================================================== #

training_data <- read.table(training_path, header=TRUE, sep="\t")
training_data <- training_data[with(training_data, order(datatype)),]
testing_data <- read.table(testing_path, header=TRUE, sep="\t")
testing_data <- testing_data[with(testing_data, order(datatype)),]

# =========================================================================== #
# Created a scaled matrix of the training data                                #
# - training|testing_D is a data frame                                        #
# - training|testing_M is a matrix                                            #
# - training|testing_S is a scaled matrix of numeric vals with similar range  #
# - training|testing_relvars_png shows the relative values of the 2 datatypes #
#                                                                             #
# The ranges of the scaled values in the training and testing matrices, saved #
# to the log file, should be approximately the same.                          #
#                                                                             #
# A heat map of the scaled values for each variable should show a high degree #
# of correlation between the data type and the value. That is, each variable  #
# should have similar values for each datatype. Since the rows of the scaled  #
# matrix has already been sorted by datatype, the plot should show a block of #
# a different colour for each datatype of each variable. Preliminary tests    #
# suggested these variables will be able to predict the data type.            #
# =========================================================================== #

training_D <- data.frame(training_data$datatype, training_data$meanrunlenC,
    training_data$freqIperbp, training_data$freqMperbp, training_data$pctC,
    training_data$pctI)
training_M <- data.matrix(training_D)
training_S <- scale(training_M)
colnames(training_S) <- c("datatype", "meanrunlenC", "freqIperbp", "freqMperbp", "pctC", "pctI")

testing_D <- data.frame(testing_data$datatype, testing_data$meanrunlenC,
    testing_data$freqIperbp, testing_data$freqMperbp, testing_data$pctC,
    testing_data$pctI)
testing_M <- data.matrix(testing_D)
testing_S <- scale(testing_M)
colnames(testing_S) <- c("datatype", "meanrunlenC", "freqIperbp", "freqMperbp", "pctC", "pctI")

cat("# Range of values in scaled matrix:\n")
cat(sprintf("  Training : %f %f\n", range(training_S)[1], range(training_S)[2]))
cat(sprintf("  Testing  : %f %f\n", range(testing_S)[1], range(testing_S)[2]))
cat("\n")

png(filename=training_relvars_png, width=500, height=800)
plot.new()
par(mfrow=c(1,1))
image(training_S, col = heat.colors(20), axes=FALSE, xlab="",ylab="", srt=45, main=runid)
axis(side=2, at=seq(0.0, 1.0, by=0.2), labels=colnames(training_S), srt=45, tick=FALSE, xpd=TRUE)
xticklabels <- sort(unique(training_data$datatype))
x <- 0.0
y <- 1.0
n <- length(xticklabels)
axis(side=1, at=seq((y-x)/(n+1), y-(y-x)/(n+1), by=(y-x)/(n+1)), labels=xticklabels, srt=45, tick=FALSE, xpd=TRUE)
box()
garbage <- dev.off()

png(filename=testing_relvars_png, width=500, height=800)
plot.new()
par(mfrow=c(1,1))
image(testing_S, col = heat.colors(20), axes=FALSE, xlab="",ylab="", srt=45, main=runid)
axis(side=2, at=seq(0.0, 1.0, by=0.2), labels=colnames(testing_S), srt=45, tick=FALSE, xpd=TRUE)
xticklabels <- sort(unique(testing_data$datatype))
x <- 0.0
y <- 1.0
n <- length(xticklabels)
axis(side=1, at=seq((y-x)/(n+1), y-(y-x)/(n+1), by=(y-x)/(n+1)), labels=xticklabels, srt=45, tick=FALSE, xpd=TRUE)
box()
garbage <- dev.off()

# =========================================================================== #
# Compute the principal components for the independent variables in the       #
# scaled training matrix.                                                     #
#                                                                             #
# Save a summary of the principal components to the log file.                 #
# The "Proportion of Variance" tells us how much of the variance is captured  #
# by PC1, then PC2, and so forth.                                             #
#                                                                             #
# Produce a line plot showing the amount of variation in the training data is #
# captured by the i-ith principal component.                                  #
#
# Plot the scores (the values of PC1, PC2, PC3) of each datatype (where       #
# minion=blue, rand=red) to see how each variable discriminates between the   #
# two datatypes.
# =========================================================================== #

training_pca <- prcomp(training_S[,2-6], center=TRUE)

cat("# Principal components from the training data:\n")
print(summary(training_pca))

png(filename=training_pccont_png, width=500, height=500)
plot.new()
par(mfrow=c(1,1))
maintext <- sprintf('%s\nAmount of variation in training data captured by\neach of the principal components', runid)
plot(training_pca, type="l", main=maintext)
garbage <- dev.off()

training_pcascores <- data.frame(training_M[,1], training_pca$x[,1:3])
colors <- c("blue", "red")
png(filename=training_pcascores_png, width=1000, height=500)
#dev.new()
par(mfrow=c(1,2))
maintext1 <- sprintf("%s\nDiscrimination of alignments by PC1 and PC2\nwhere minion alignments are in blue, random in red", runid)
maintext2 <- sprintf("%s\nDiscrimination of alignments by PC2 an PC3\nwhere minion alignments are in blue, random in red", runid)
plot(training_pcascores[,2], training_pcascores[,3], pch=20, col=colors[training_pcascores[,1]], xlab="PC1", ylab="PC2", main=maintext1)
plot(training_pcascores[,3], training_pcascores[,4], pch=20, col=colors[training_pcascores[,1]], xlab="PC2", ylab="PC3", main=maintext2)
garbage <- dev.off()

# =========================================================================== #
# Print some properties of the scaled training data that will be used to      #
# build the logistic regression model.                                        #
# =========================================================================== #

training_R <- data.frame(training_S)
#colnames(training_R)

cat("\n# Summary of the scaled training data:\n")
summary(training_R)
cat("\n# Standard deviation of the scaled training data:\n")
sapply(training_R, sd)
cat("\n# Contingency table of the scaled training data:\n")
xtabs(~datatype, data=training_R)

# =========================================================================== #
# Separate the training data into two randomly-selected subsets:              #
# 2/3 for building the linear regression model and 1/3 for testing it         #
#                                                                             #
# - training_R a data frame of the scaled values                              #
# - sapply prints a summary and standard deviation of the values              #
# - xtab prints a contingency table of the category outcome and predictors    #
#   which we should check to ensure there are no zero cells.                  #
#                                                                             #
# Modified from http://www.gettinggeneticsdone.com/2011/02/                   #
#   split-data-frame-into-testing-and.html                                    #
# =========================================================================== #

splitdf <- function(dataframe, testpct, seed=NULL) {
    if (!is.null(seed)) set.seed(seed)
    index <- 1:nrow(dataframe)
    trainindex <- sample(index, trunc(length(index)*testpct/100.0))
    trainset <- dataframe[-trainindex, ]
    testset <- dataframe[trainindex, ]
    list(trainset=trainset,testset=testset)
}

training_A <- data.frame(training_S)
training_A$datatype <- as.integer(training_data$datatype == "minion")
training_splits <- splitdf(training_A, 33, seed=808)
training_modeldata <- training_splits$trainset
training_modeltest <- training_splits$testset

# =========================================================================== #
# Create the model from the modeldata and print the accuracy of predictions   #
# for the modeltest data.                                                     #
#                                                                             #
# Allocate the test data point to "minion" if the prediction of being "minion"#
# (i.e., of being datatype=1) has a probability >= probability_threshold      #
# which has been set to the quite stringent value of 0.99.                    #
# It means we have a low false positive rate (i.e., low chance of incorrectly #
# classifying a real read as random) but a high false negative rate (i.e., a  #
# high chance of incorrectly classifying a random read as real).              #
#                                                                             #
# Modified from                                                               #
# https://stat.ethz.ch/pipermail/r-help/2006-September/113037.html            #
# =========================================================================== #
 
training_modeldata.lg <- glm(datatype ~ freqIperbp+meanrunlenC+pctC+pctI+freqMperbp, family=binomial, data=training_modeldata)
training_modeltest.probabilityisminion <- predict(training_modeldata.lg, training_modeltest, type="response")
training_pctcorrect <- sum(as.integer(training_modeltest.probabilityisminion >= probability_threshold) == training_modeltest$datatype) / length(training_modeltest$datatype) * 100.0

cat("\nLogistic regression model constructed from 2/3 of the training data\n")
cat(sprintf("  Number of data points used to construct model: %d\n", length(training_modeldata[,1])))
summary(training_modeldata.lg)
cat("\nAccuracy of the logistic regression model:\n")
cat(sprintf("  Number of training alignments: %d\n", length(training_data[,1])))
cat(sprintf("  Number of training alignments used to build model: %d\n", length(training_modeldata[,1])))
cat(sprintf("  Number of training alignments used to test model: %d\n", length(training_modeltest[,1])))
cat(sprintf("  Prediction accuracy : %.2f%%\n", training_pctcorrect))

# =========================================================================== #
# Use the same model to predict the datatype of the testing set.              #
# =========================================================================== #

testing_A <- data.frame(testing_S)
testing_A$datatype <- as.integer(testing_data$datatype == "minion")
testing.probabilityisminion <- predict(training_modeldata.lg, testing_A, type="response")
testing_isminion <- as.integer(as.integer(testing.probabilityisminion >= probability_threshold) == testing_A$datatype)
testing_isrand <- as.integer(as.integer(testing.probabilityisminion < probability_threshold) == testing_A$datatype)
testing_pctminion <- sum(testing_isminion) / length(testing_A[,1]) * 100.0
testing_pctrand <- sum(testing_isrand) / length(testing_A[,1]) * 100.0

cat("\nPredictions for the testing set:\n")
cat(sprintf("  Number of testing alignments: %d\n", length(testing_data[,1])))
cat(sprintf("  Number of testing alignments predicted as datatype=minion: %.2f%%\n", testing_pctminion))
cat(sprintf("  Number of testing alignments predicted as datatype=random: %.2f%%\n", testing_pctrand))

# =========================================================================== #
# Create a new data.frame with the updated isnonrandomaln predictions and     #
# save it as an alignclass file, which is the first 11 columns of the original#
# testing_data plus an "isnonrandomaln" column.                               #
# =========================================================================== #

result_data <- data.frame(
    testing_data$runid, testing_data$readid, testing_data$readtype, testing_data$readclass,
    testing_data$datatype, testing_data$mapprog, testing_data$mapparams, testing_data$samflag,
    testing_data$refcontigid, testing_data$refcontigpos1, testing_data$mapq, testing_isminion)
colnames(result_data) <- c("runid", "readid", "readtype", "readclass", "datatype", "mapprog", "mapparams", "samflag", "refcontigid", "refcontigpos1", "mapq", "isnonrandomaln")
write.table(result_data, file=result_outpath, append=FALSE, quote=FALSE, sep="\t", row.names=FALSE)
quit(status=0)

# =========================================================================== #
