# Plot HMM output
# This is probably the most inefficient code I ever wrote,
#   but it works so I'm unmotivated to change it :(
# V0.1

install.packages("ggplot2", lib = "./", repos = "http://cran.wu.ac.at")
library(ggplot2,lib.loc = "./")
install.packages("vegan", lib = "./", repos = "http://cran.wu.ac.at")
library(vegan,lib.loc = "./")

PlotHmmOutput <- function(HMM.OUT, type = "eval", gene = "bchX") {
  require(ggplot2)
  HMM.OUT$V24 <- log(HMM.OUT$V7)
  if (type == "eval") {
    ggplot(HMM.OUT, aes(V24)) +
      geom_freqpoly(binwidth = 2) +
      geom_vline (xintercept = log(1e-50), colour = "red") +
      labs(x = "ln E-value", y = "Count", title = paste(gene, "HMM model E-value distribution"))
  } else {
    ggplot(HMM.OUT, aes(V8)) +
      geom_freqpoly(binwidth = 2) +
      geom_vline (xintercept = 150, colour = "red") +
      labs(x = "Score", y = "Count", title = paste(gene, "HMM model Score distribution"))
  }
}

pdf("nifH_bch_hmmEvaluation.pdf")
HMM.OUT <- read.table(commandArgs()[3])
PlotHmmOutput(HMM.OUT[grep("nifH", HMM.OUT[,1]),], 
              type = "eval", gene = "nifH")
PlotHmmOutput(HMM.OUT[grep("nifH", HMM.OUT[,1]),], 
              type = "score",
              gene = "nifH")
PlotHmmOutput(HMM.OUT[HMM.OUT[, 1] == "bchX", ], 
              type = "eval", gene = "bchX")
PlotHmmOutput(HMM.OUT[HMM.OUT[, 1] == "bchX", ], 
              type = "score",
              gene = "bchX")
PlotHmmOutput(HMM.OUT[HMM.OUT[, 1] == "chlL-bchL", ],
              type = "eval", gene = "chlL-bchL")
PlotHmmOutput(HMM.OUT[HMM.OUT[, 1] == "chlL-bchL", ], 
              type = "score", 
              gene = "chlL-bchL")

OTUnums <- sort(as.numeric(gsub("[^0-9]+", "", unique(HMM.OUT[, 4]))))
scoreTable <- matrix(0, nrow = length(OTUnums),ncol = 3)
rownames(scoreTable) <- OTUnums
colnames(scoreTable) <- c("nifH", "chlL", "bchX")
evalTable <- scoreTable
nifHSet <- grep("nifH", HMM.OUT[, 1])
chlLSet <- grep("chlL", HMM.OUT[, 1])
bchXSet <- grep("bchX", HMM.OUT[, 1])
for(i in seq(nrow(scoreTable))){
	thisOTUset <- which(as.numeric(gsub("[^0-9]+", "", HMM.OUT[, 4])) == rownames(scoreTable)[i])
	scoreTable[i, 1] <- HMM.OUT[intersect(thisOTUset, nifHSet)[1], 8]
	scoreTable[i, 2] <- HMM.OUT[intersect(thisOTUset, chlLSet)[1], 8]
	scoreTable[i, 3] <- HMM.OUT[intersect(thisOTUset, bchXSet)[1], 8]
	evalTable[i, 1] <- HMM.OUT[intersect(thisOTUset, nifHSet)[1], 7]
	evalTable[i, 2] <- HMM.OUT[intersect(thisOTUset, chlLSet)[1], 7]
	evalTable[i, 3] <- HMM.OUT[intersect(thisOTUset, bchXSet)[1], 7]

}
assignments <- rownames(which(apply(scoreTable, 1, rank) == 3, arr.ind = TRUE))

par(mfrow = c(1, 1))
plot(scoreTable[, 2], 
     scoreTable[, 1],
     pch = 16,
     col = as.numeric(as.factor(assignments)),
     xlab = colnames(scoreTable)[2],
     ylab = "nifH",
     main = "separation of ChlL")
abline(a = 0, 
       b = 1)
legend("topleft",
       levels(as.factor(assignments)),
       col = c(1:length(levels(as.factor(assignments)))), 
       pch = 16)
plot(scoreTable[, 3], 
     scoreTable[, 1],
     pch = 16,
     col = as.numeric(as.factor(assignments)),
     xlab = colnames(scoreTable)[3],
     ylab = "nifH", 
     main = "separation of BchX")
abline(a = 0, 
       b = 1)
legend("topleft",
       levels(as.factor(assignments)),
       col = c(1:length(levels(as.factor(assignments)))),
       pch = 16)

pcoA <- rda(scoreTable)
plot(pcoA$CA$u[, 1] * pcoA$CA$eig[1],
     pcoA$CA$u[, 2] * pcoA$CA$eig[2],
     col = as.numeric(as.factor(assignments)),
     pch = 16,
     xlab = paste("PcoA 1 (", round(pcoA$CA$eig[1] / sum(pcoA$CA$eig), 4) * 100, "%)"),
     ylab = paste("PcoA 2 (", round(pcoA$CA$eig[2] / sum(pcoA$CA$eig), 4) * 100, "%)"))
legend("topleft",
       levels(as.factor(assignments)),
       col = c(1:length(levels(as.factor(assignments)))),
       pch = 16)
dev.off()
