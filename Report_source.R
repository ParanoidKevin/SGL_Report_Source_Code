#install.packages("MASS")
#install.packages("SGL")
#install.packages("glmnet")
library(glmnet)
library(MASS)
library(SGL)
set.seed(123)

##### define a function that calculates the performance of SGL in a single MC run
SGL_LASSO_performance <- function(num_obs, num_par, num_grp, g, true_beta, SNR = 2){
  # generate x
  x <- mvrnorm(1, rep(0,num_obs), Sigma = diag(num_obs))
  for (i in 1:(num_par-1)){
    x <- cbind(x, mvrnorm(1, rep(0,num_obs), Sigma = diag(num_obs)))
  }
  # generate residuals
  e <- SNR*mvrnorm(1, rep(0,num_obs), Sigma = diag(num_obs))
  # compute y for our model
  y <- x %*% true_beta + e
  
  # fit SGL model
  data1   <- list(x=x, y=y)
  index1  <- ceiling(1:num_par/(num_par/num_grp))
  fit1    <- SGL(data1, index1)
  betahat <- fit1$beta
  # fit LASSO model
  fit2    <- glmnet(x,y,family="mgaussian")
  betahat2 <- as.matrix(fit2$beta)
  
  # choose the penalty parameters so that the number of nonzero coefficients is g*5 if number exactly matches
  # otherwise, select the minimum extra number of nonzero coefficients
  nonzero <- g*5
  
  ####### precision for SGL
  if (nonzero %in% colSums(betahat != 0)) {
    index2     <- which(colSums(betahat != 0) == nonzero)
    precisions <- vector(mode = "numeric")
    for (i in 1:length(index2)){
      beta_est    <- betahat[, index2[i]]
      #best_lambda <- fit1$lambdas[index2[i]] 
      # calculate the proportion of correct nonzero coefficient indentifications
      precisions <- c(precisions, 1-length(setdiff(which(beta_est != 0), which(true_beta != 0)))/length(which(true_beta != 0)))
      precision  <- max(precisions)
    }
  } else{ # if extra nonzero coefficients are chosen
    precisions <- vector(mode = "numeric")
    num_nonzeros_pos <- colSums(betahat != 0) - g*5
    num_nonzeros     <- min(num_nonzeros_pos[num_nonzeros_pos > 0]) + g*5
    index3           <- which(num_nonzeros_pos == min(num_nonzeros_pos[num_nonzeros_pos > 0]))
    for (i in 1:length(index3)){
      beta_est         <- betahat[, index3[i]]
      #best_lambda      <- fit1$lambdas[index3] 
      # calculate the proportion of correct nonzero coefficient indentifications
      num_match        <- length(which(beta_est != 0)) - length(setdiff(which(beta_est != 0), which(true_beta != 0)))
      precisions       <- c(precisions, num_match/num_nonzeros)
      precision        <- max(precisions)
    }
  }
  
  ######### pricision for LASSO
  if (nonzero %in% colSums(betahat2 != 0)) {
    indexL      <- which(colSums(betahat2 != 0) == nonzero)
    precisionsL <- vector(mode = "numeric")
    for (i in 1:length(indexL)){
      beta_estL    <- betahat2[, indexL[i]]
      #best_lambda <- fit2$lambdas[index2[i]] 
      # calculate the proportion of correct nonzero coefficient indentifications
      precisionsL     <- c(precisionsL, 1-length(setdiff(which(beta_estL != 0), which(true_beta != 0)))/length(which(true_beta != 0)))
      precisionLASSO  <- max(precisionsL)
    }
  } else{ # if extra nonzero coefficients are chosen
    precisionsL <- vector(mode = "numeric")
    num_nonzeros_posL <- colSums(betahat2 != 0) - g*5
    num_nonzerosL     <- min(num_nonzeros_posL[num_nonzeros_posL > 0]) + g*5
    indexL            <- which(num_nonzeros_posL == min(num_nonzeros_posL[num_nonzeros_posL > 0]))
    for (i in 1:length(indexL)){
      beta_estL         <- betahat2[, indexL[i]]
      # calculate the proportion of correct nonzero coefficient indentifications
      num_matchL     <- length(which(beta_estL != 0)) - length(setdiff(which(beta_estL != 0), which(true_beta != 0)))
      precisionsL    <- c(precisionsL, num_matchL/num_nonzerosL)
      precisionLASSO <- max(precisionsL)
    }
  }
  
  result <- list(precision, precisionLASSO)
  return(result)
}


####### define a function to define the true beta according to g
true_betas <- function(g,p,m){
  if (g == 1){
    true_beta <- c(c(1,2,3,4,5), rep(0, p-5))
  } else if (g == 2){
    true_beta <- c(c(1,2,3,4,5), rep(0, p/m-5), 
                   c(1,2,3,4,5), rep(0, p-5-p/m))
  } else if (g == 3){
    true_beta <- c(c(1,2,3,4,5), rep(0, p/m-5),
                   c(1,2,3,4,5), rep(0, p/m-5),
                   c(1,2,3,4,5), rep(0, p-5-2*p/m))
  }
}

performance_SGL   <- matrix(nrow = 4, ncol = 3)
performance_LASSO <- matrix(nrow = 4, ncol = 3)
########## g = 1,2,3
for (g in 1:3){
  
  ############################## scenario 1
  n <- 20
  p <- 200
  m <- 5
  true_beta <- true_betas(g,p,m)
  
  # run 100 MC runs of the procedure and take the average performance for SGL
  total <- 0
  for (i in 1:100){
    total = total + SGL_LASSO_performance(n, p, m, g, true_beta, SNR = 2)[1][[1]]
  }
  avg_precisionSGL1 <- total/100
  
  # run 100 MC runs of the procedure and take the average performance for LASSO
  total <- 0
  for (i in 1:100){
    total = total + SGL_LASSO_performance(n, p, m, g, true_beta, SNR = 2)[2][[1]]
  }
  avg_precisionLASSO1 <- total/100
  
  print(paste("----------------", round((4*(g-1)+1)*100/12,2), "% Completed", "-----------------"))
  
  ############################## scenario 2
  n <- 30
  p <- 500
  m <- 50
  true_beta <- true_betas(g,p,m)
  
  # run 100 MC runs of the procedure and take the average performance for SGL
  total <- 0
  for (i in 1:100){
    total = total + SGL_LASSO_performance(n, p, m, g, true_beta, SNR = 2)[1][[1]]
  }
  avg_precisionSGL2 <- total/100
  
  # run 100 MC runs of the procedure and take the average performance for LASSO
  total <- 0
  for (i in 1:100){
    total = total + SGL_LASSO_performance(n, p, m, g, true_beta, SNR = 2)[2][[1]]
  }
  avg_precisionLASSO2 <- total/100
  
  print(paste("----------------", round((4*(g-1)+2)*100/12,2), "% Completed", "-----------------"))
  
  ############################# scenario 3
  n <- 50
  p <- 1000
  m <- 20
  true_beta <- true_betas(g,p,m)
  
  # run 100 MC runs of the procedure and take the average performance for SGL
  total <- 0
  for (i in 1:100){
    total = total + SGL_LASSO_performance(n, p, m, g, true_beta, SNR = 2)[1][[1]]
  }
  avg_precisionSGL3 <- total/100
  
  # run 100 MC runs of the procedure and take the average performance for LASSO
  total <- 0
  for (i in 1:100){
    total = total + SGL_LASSO_performance(n, p, m, g, true_beta, SNR = 2)[2][[1]]
  }
  avg_precisionLASSO3 <- total/100
  
  print(paste("----------------", round((4*(g-1)+3)*100/12,2), "% Completed", "-----------------"))
  
  ############################ scenario 4
  n <- 100
  p <- 2000
  m <- 50
  true_beta <- true_betas(g,p,m)
  
  # run 100 MC runs of the procedure and take the average performance for SGL
  total <- 0
  for (i in 1:100){
    total = total + SGL_LASSO_performance(n, p, m, g, true_beta, SNR = 2)[1][[1]]
  }
  avg_precisionSGL4 <- total/100
  
  # run 100 MC runs of the procedure and take the average performance for LASSO
  total <- 0
  for (i in 1:100){
    total = total + SGL_LASSO_performance(n, p, m, g, true_beta, SNR = 2)[2][[1]]
  }
  avg_precisionLASSO4 <- total/100
  
  print(paste("----------------", round(g*100/3,2), "% Completed", "-----------------"))
  
  performance_SGL[,g]   <- c(avg_precisionSGL1, avg_precisionSGL2, avg_precisionSGL3, avg_precisionSGL4)
  performance_LASSO[,g] <- c(avg_precisionLASSO1, avg_precisionLASSO2, avg_precisionLASSO3, avg_precisionLASSO4)
}

performance_SGL
performance_LASSO
