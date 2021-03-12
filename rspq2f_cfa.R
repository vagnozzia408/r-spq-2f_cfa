# AUTHOR: Anna Marie Vagnozzi, M.S.
# CREATED: 5/7/2020
# LAST UPDATED: 5/28/2020
# DESCRIPTION: This R script was written to accompany the paper "Issues with the Revised Study
# Process Questionnaire (R-SPQ-2F) in Undergraduate Anatomy & Physiology Students" by
# Johnson, Gallagher, and Vagnozzi.

# Require Packages
require(MVN) # sROC 0.1-2
require(matrixStats)
require(lavaan) # version 0.6-6

# Set working directory. Must contain R-SPQ-2F item response data.
setwd("C:/Users/avagnoz/Dropbox/Clemson/ESED/GRA/SPQ_CFA/R_CFA_Script/")
rm(list=ls())

# Read in the R-SPQ-2F item-level responses.
# This CFA is conducted for complete responses only.
spq_data = read.table("cleaned_spq_data.txt", header = TRUE)
spq_data = subset(spq_data,select=-c(ID))

# Correlation and covariance matrices.
cor_matrix <- cor(spq_data)
stdevs <- as.matrix(colSds(as.matrix(spq_data)))
cov_matrix <- cov(spq_data)

# Check for univariate/multivariate normality.
mvn_results <- mvn(data = spq_data, mvnTest = "mardia", covariance=FALSE)

mvn_results$Descriptives
mvn_results$univariateNormality
mvn_results$multivariateNormality

# CFA MODEL #
spq_model <- 'deep =~ SPQ1 + SPQ2 + SPQ5 + SPQ6 + SPQ9 + SPQ10 + SPQ13 + SPQ14 + SPQ17 + SPQ18
              surface =~ SPQ3 + SPQ4 + SPQ7 + SPQ8 + SPQ11 + SPQ12 + SPQ15 + SPQ16 + SPQ19 + SPQ20'

fit <- list()

# MAXIMUM LIKELIHOOD METHOD #
fit$ml <- cfa(spq_model, data=spq_data, estimator="MLM", std.lv=TRUE)
#summary(fit$ml, fit.measures=TRUE)
# Produce the completely standardized solution.
standardizedSolution(fit$ml)

# WEIGHTED LEAST SQUARES METHOD #
fit$wls <- cfa(spq_model, data=spq_data, estimator="WLS", std.lv=TRUE)
#summary(fit$wls, fit.measures=TRUE)
# Produce the completely standardized solution.
standardizedSolution(fit$wls)

# EVALUATE MODEL FITS #
fitindices <- c("npar",
                "chisq", "df", "pvalue",
                "chisq.scaled", "df.scaled", "pvalue.scaled",
                "cfi", "cfi.scaled", 
                "tli", "tli.scaled", 
                "rmsea", "rmsea.scaled", 
                "srmr")

get_fits <- function(fit, fitidxs=fitindices, digits=4){
  x <- fitMeasures(fit)
  round(x[fitidxs],digits)
}

sapply(fit, function(X) get_fits(X))

# R^2 Values
r2_ml <- inspect(fit$ml, 'r2')
r2_wls <- inspect(fit$wls, 'r2')

# Produces standardized residuals for the model (see if |z|>2.0 indicating possible source of misfit)
resid(fit$ml,type='standardized') 
resid(fit$wls, type='standardized')

##### REMOVE OUTLIERS and Fit Again #####

# Note: The multivariateOutlierMethod option in the MVN package identifies a different number of multivariate
# outliers each time it is run. The author of the package was contacted on 5/11/2020. Updates are in progress.
mvnoutlrs <- mvn(data = spq_data, mvnTest = "mardia", covariance=FALSE, multivariateOutlierMethod = "quan", showOutliers=TRUE, showNewData=TRUE)

# Data set with no multivariate outliers.
spq_no_outlrs <- mvnoutlrs$newData

# Fit the model again with no multivariate outliers.

# MAXIMUM LIKELIHOOD METHOD #
fit$mlno <- cfa(spq_model, data=spq_no_outlrs, estimator="MLM", std.lv=TRUE)
#summary(fit$mlno, fit.measures=TRUE)
standardizedSolution(fit$mlno)

# WEIGHTED LEAST SQUARES METHOD #
fit$wlsno <- cfa(spq_model, data=spq_no_outlrs, estimator="WLS", std.lv=TRUE)
#summary(fit$wlsno, fit.measures=TRUE)
standardizedSolution(fit$wlsno)

# Evaluate again:
sapply(fit, function(X) get_fits(X))

# Remove the univariate outliers.
# Identify for each item using boxplot.stats(spq_no_outlrs$SPQ1)$out
spq_no_outlrs2 <- spq_no_outlrs[!(spq_no_outlrs$SPQ1==1),]
spq_no_outlrs2 <- spq_no_outlrs2[!(spq_no_outlrs2$SPQ2==1),]
spq_no_outlrs2 <- spq_no_outlrs2[!(spq_no_outlrs2$SPQ7==4),]
spq_no_outlrs2 <- spq_no_outlrs2[!(spq_no_outlrs2$SPQ7==5),]
spq_no_outlrs2 <- spq_no_outlrs2[!(spq_no_outlrs2$SPQ11==5),]
spq_no_outlrs2 <- spq_no_outlrs2[!(spq_no_outlrs2$SPQ12==5),]
spq_no_outlrs2 <- spq_no_outlrs2[!(spq_no_outlrs2$SPQ13==1),]
spq_no_outlrs2 <- spq_no_outlrs2[!(spq_no_outlrs2$SPQ15==4),]
spq_no_outlrs2 <- spq_no_outlrs2[!(spq_no_outlrs2$SPQ15==5),]
spq_no_outlrs2 <- spq_no_outlrs2[!(spq_no_outlrs2$SPQ17==5),]
spq_no_outlrs2 <- spq_no_outlrs2[!(spq_no_outlrs2$SPQ18==1),]

# Fit the model once more with both univariate and multivariate outliers removed.

# MAXIMUM LIKELIHOOD METHOD #
fit$mlno2 <- cfa(spq_model, data=spq_no_outlrs2, estimator="MLM", std.lv=TRUE)
#summary(fit$mlno2, fit.measures=TRUE)
standardizedSolution(fit$mlno2)

# WEIGHTED LEAST SQUARES METHOD #
fit$wlsno2 <- cfa(spq_model, data=spq_no_outlrs2, estimator="WLS", std.lv=TRUE)
#summary(fit$wlsno2, fit.measures=TRUE)
standardizedSolution(fit$wlsno2)

# Final Fit Indices
final_fits <- sapply(fit, function(X) get_fits(X))
final_fits

require(semTools)
reliability(fit$ml)
reliability(fit$wls)
