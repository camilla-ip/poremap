# =========================================================================== #
# poremapcoverage.R                                                           #
#                                                                             #
# Plot the read depth for random and non-random alignments across the         #
# refcontig.                                                                  #
#                                                                             #
# Usage:                                                                      #
#                                                                             #
#     Rscript poremapcoverage.R \                                             #
#         runid \                                                             #
#         refcontigid \                                                       #
#         refcontiglen \                                                      #
#         depthpath \                                                         #
#         outdir \                                                            #
#         outprefix                                                           #
#                                                                             #
# where                                                                       #
#     depthpath is a tab-sep file with 7 columns:                             #
#         pos, rawgood, rawrand, rawboth, wingood, winrand, winboth           #
#                                                                             #
# Output:                                                                     #
#                                                                             #
#   outdir/outprefix_coverage.png: Line plots for 3 windowed depths           #
#                                                                             #
# Camilla Ip                                                                  #
# camilla.ip@well.ox.ac.uk                                                    #
# May 2015                                                                    #
# =========================================================================== #

# =========================================================================== #
# Constants
# =========================================================================== #

# Palettes from http://blog.mollietaylor.com/2012/10/color-blindness-and-palette-choice.html
alpha <- 0.5
cbPaletteLt <- c(rgb(0.40,0.76,0.65, alpha), rgb(0.99,0.55,0.38, alpha), rgb(0.55,0.63,0.80, alpha), rgb(0.91,0.54,0.76, alpha),
                 rgb(0.65,0.85,0.33, alpha), rgb(1.00,0.85,0.18, alpha), rgb(0.90,0.77,0.58, alpha), rgb(0.70,0.70,0.70, alpha))
cbPaletteDk <- c(rgb(0.11,0.62,0.47, alpha), rgb(0.85,0.37,0.01, alpha), rgb(0.46,0.44,0.70, alpha), rgb(0.91,0.16,0.54, alpha),
                 rgb(0.40,0.65,0.12, alpha), rgb(0.90,0.67,0.01, alpha), rgb(0.65,0.46,0.11, alpha), rgb(0.26,0.26,0.26, alpha))

nonrandraw_color <- cbPaletteLt[1]
randraw_color <- cbPaletteLt[2]
bothraw_color <- cbPaletteLt[3]

nonrandwin_color <- cbPaletteDk[1]
randwin_color <- cbPaletteDk[2]
bothwin_color <- cbPaletteDk[3]

# =========================================================================== #
# Command-line arguments
# =========================================================================== #

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
    txt <- sprintf("Usage: Rscript poremapcoverage runid refcontigid refcontiglen depthpath outdir outprefix")
    print(txt)
    quit(status=1)
}
runid <- args[1]
refcontigid <- args[2]
refcontiglen <- as.integer(args[3])
depthpath <- args[4]
outdir <- args[5]
outprefix <- args[6]

if (refcontiglen > 1000000) {
    sampleinterval <- 1000
    unitinterval <- 1000000
    units <- 'Mb'
} else {
    sampleinterval <- 100
    unitinterval <- 100
    units <- 'Kb'
}
winflanksz <- 10

# =========================================================================== #
# Output file paths                                                           #
# =========================================================================== #

plot_path <- sprintf("%s/%s_coverage.png", outdir, outprefix)

# =========================================================================== #
# Create 6 line plots on one graph                                            #
# =========================================================================== #

data <- read.table(depthpath, header=TRUE, sep="\t", blank.lines.skip=TRUE)

maintext <- sprintf("Read depth: %s, %s", runid, refcontigid)
xlabeltext <- sprintf("Genome position (%s)", units)
legendlabels <- c('good', 'rand', 'both')
legendcolors <- c(cbPaletteDk[1], cbPaletteDk[2], cbPaletteDk[3])
ymax <- max(data$rawboth, data$winboth)

png(filename=plot_path,  width=1000, height=500)
par(mfrow=c(1,1))
plot(NULL, xlab=xlabeltext, ylab="Read depth", main=maintext, xlim=c(0, refcontiglen/unitinterval), ylim=c(0, ymax))
lines(data$pos[(data$pos %% sampleinterval)==0]/unitinterval, data$wingood[(data$pos %% sampleinterval)==0], col=legendcolors[1], type='l')
lines(data$pos[(data$pos %% sampleinterval)==0]/unitinterval, data$winrand[(data$pos %% sampleinterval)==0], col=legendcolors[2], type='l')
lines(data$pos[(data$pos %% sampleinterval)==0]/unitinterval, data$winboth[(data$pos %% sampleinterval)==0], col=legendcolors[3], type='l')
legend(refcontiglen/unitinterval*0.93, ymax, legendlabels, pch=15, pt.cex=1.5, col=legendcolors, cex=1.0, bty="n")
garbage <- dev.off()

# =========================================================================== #
