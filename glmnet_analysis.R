# How to install package from source. 
# install.packages("/cloud/home/r1816512/CERVICAL_CANCER/randomForestSRC_3.2.3.tar.gz", repos = NULL, type = "source")

# Install from CRAN mirror nearest to Uganda, Netherlands, Evoluso.com

install.packages("glmnet", repos = "https://mirrors.evoluso.com/CRAN")

install.packages("randomForestSRC", repos = "https://mirrors.evoluso.com/CRAN")

library(randomForestSRC)
packageVersion("randomForestSRC")
# [1] ‘3.2.3’

# fit a GLM with lasso or elasticnet regularization

# First detach survival package to avoid conflict in the Matrix library also found in glmnet
detach("package:survival", unload = TRUE)
library(Matrix)
library(glmnet)
# Loaded glmnet 4.1-8
packageVersion("glmnet")
#[1] ‘4.1.8’
 ?glmnet()

library(lattice)
library(latticeExtra)


TCGA_GTEx_CESC_OS_gene_read_counts_306_Tumor_22_Normal_7910 <- readRDS("TCGA_GTEx_CESC_OS_gene_read_counts_306_Tumor_22_Normal_7910.rds")
dim(TCGA_GTEx_CESC_OS_gene_read_counts_306_Tumor_22_Normal_7910)
# [1] 7910  328

# Specify settings for the glmnet analysis.
nfolds <- 20                    # Number of cross-validation folds.
alpha  <- 0.95                  # Elastic net mixing parameter.
lambda <- 10^(seq(0, -2, -0.05)) # Lambda sequence.
set.seed(1)

# Create a factor response variable
g2 <- factor(rep(c("Tumor", "Normal"), c(306, 22)))

# Create a matrix for predictor variables
x <- as.matrix(TCGA_GTEx_CESC_OS_gene_read_counts_306_Tumor_22_Normal_7910)

# Fit a binomial regression model using glmnet
fit.glmnet <- glmnet(x, g2, family = "binomial", lambda = lambda, alpha = alpha)

# Error in glmnet(x, g2, family = "binomial", lambda = lambda, alpha = alpha) : 
# number of observations in y (328) not equal to the number of rows of x (7910)


fit2r = glmnet(x,g2, family = "binomial", relax=TRUE)


# This is the cross-validation step.
r <- system.time(out.cv.glmnet <-
                   cv.glmnet(X,y,family = "binomial",type.measure = "class",
                             alpha = alpha,nfolds = nfolds,lambda = lambda))
lambda <- out.cv.glmnet$lambda
cat(sprintf("Cross-validation took %0.2f seconds.\n",r["elapsed"]))
# Cross-validation took 0.59 seconds.
