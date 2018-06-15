# Importing data --------
require(devtools)
source_url('https://raw.githubusercontent.com/FredHasselman/toolboxR/master/C-3PR.R')

## If you dowloaded the csv file to your harddrive use this code:
RPPdata<-read.csv('data-raw/rpp_data.csv',stringsAsFactors=F )
RPPdata<-df.Clean(RPPdata)
RPPdata<-RPPdata$df

sortPvals <- sort(RPPdata$T.pval.recalc.O)
sum(sortPvals[1:(length(sortPvals)-1)] == sortPvals[2:length(sortPvals)])
which(sortPvals[1:(length(sortPvals)-1)] == sortPvals[2:length(sortPvals)])

RPPdata <- dplyr::filter(RPPdata, !is.na(T.pval.USE.O), !is.na(T.pval.USE.R))

# We have 99 studies for which p-values and effect sizes could be calculated
nrow(RPPdata)
# We have 97 studies for which p-values of the original effect were significantly below .05
idOK <- complete.cases(RPPdata$T.r.O,RPPdata$T.r.R)
sum(idOK)

RPPdata$Power.Rn <- as.numeric(RPPdata$Power.R)

# Convert F tests to T tests
for(i in 1:nrow(RPPdata)) {
  if(RPPdata$T.Test.Statistic.O[i] == "F" & RPPdata$T.df1.O[i] == 1) {
    RPPdata$T.Test.Statistic.O[i] <- "t"
    RPPdata$T.Test.value.O[i] <- sqrt(RPPdata$T.Test.value.O[i])
  }

  if(RPPdata$T.Test.Statistic.R[i] == "F" & RPPdata$T.df1.O[i] == 1) {
    RPPdata$T.Test.Statistic.R[i] <- "t"
    RPPdata$T.Test.value.R[i] <- sqrt(RPPdata$T.Test.value.R[i]) * (1 - 2*(RPPdata$Direction.R[i] == "opposite"))
  }
}

# Correcting test statistics --------------
RPPdata$test_statistic <- RPPdata$T.Test.Statistic.R
RPPdata$test_statistic[RPPdata$test_statistic == "Chi2"] <- "chisq"
RPPdata$test_statistic[RPPdata$test_statistic == "beta"] <- "t"
RPPdata$test_statistic[RPPdata$test_statistic %in% c("r")] <- "correlation"

# Screening partial studies ------------
RPPdata <- subset(RPPdata, !is.na(T.Test.value.O))
RPPdata <- subset(RPPdata, !is.na(T.Test.value.R))

# Saving dataset
osc_data <- RPPdata
use_data(osc_data, overwrite = TRUE)


