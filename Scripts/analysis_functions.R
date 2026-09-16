### Analysis functions for three types of models to analyze cross sectional questionnaire
### data 
### (1) TS (total score) estimator: regress total score on treatment + covariates
### (2) IL-ME (item-level main effects): Item-level GEE to estimate treatment effect on 
###     total score while adjusting for item
### (3) IL-IS (item-level item-specific): Item-level GEE to estimate item-specific treatment effects,
###.    which we will combine together to estimate the treatment effect on the total score

library(geepack)
library(dplyr)
library(GenOrd)
library(tidyr)
library(MASS)
library(stringr)


### Helper functions
logit <- function(x) {
  return(log(x/(1-x)))
}

expit <- function(x) {
  return(exp(x)/(1+exp(x)))
}

## Derivative of expit(x^T\beta) with respect to beta, evaluated at beta hat
## where x is a vector, i.e. a row of the design matrix, and beta is also a vector
expit_derivative <- function(x, beta) {
  x_vec <- as.matrix(x)
  return(x_vec*exp(c(t(x_vec)%*%beta))/((1+exp(c(t(x_vec)%*%beta)))^2))
}

## Calculate the variance of a weighted sum of random variables given their
## covariance matrix (varmat) and weights (a)
sum_variance <- function(varmat, a) {
  sum <- 0
  n <- nrow(varmat)
  ## Sum the individual variances
  for (i in 1:n) {
    sum <- sum + (a[i]^2)*varmat[i,i]
  }
  ## Add the cross covariances
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      sum <- sum + 2*(a[i]*a[j]*varmat[i,j])
    }
  }
  return(sum)
}

## Return person-level residuals for the total score given a geeglm object, estimator
## (item specific or standard), the data, and the max value of the item response
## if a binomial family model was used,
## K = max value of response for support {0,1,...,K}
get_residuals <- function(model, estimator, data, K=NA, type="id") {
  newdata <- data
  if (model$family$family=="binomial" & model$family$link=="logit") {
    if (estimator=="item_itemspecific" | estimator=="item_maineffects") {
      if (is.na(K)) {
        K <- (data$response + data$failures)[1]
      }
      newdata$fitted_response <- K*model$fitted.values[,1] # The fitted values are the fitted p = expit(xTbeta)
      newdata$residual <- newdata$response - newdata$fitted_response
    }
    if (estimator=="standard") {
      if (is.na(K)) {
        K <- (data$score + data$failures)[1]
      }
      newdata$fitted_score <- K*model$fitted.values[,1]
      newdata$residual <- newdata$score - newdata$fitted_score
    }
  } else { # gaussian family
    newdata$residual <- residuals(model, "response")[,1]
  }
  
  id_residuals <- newdata %>% group_by(id) %>% summarize(residual = sum(residual)) %>% ungroup()
  if (type=="item") {
    return(newdata)
  } else {
    return(id_residuals)
  }
}

### Calculate treatment effect estimates for the three model types
### data: Dataframe with item responses or score values and any covariates. If estimator
###.      is "item_itemspecific" or "item_maineffects", should be long format with a row
###.      for each item, and the variable for the item response should be called "response".
###.      If estimator is "standard", there should be a row per participant and the variable
###.      for the total score should be called "score".
### 
### estimator: "item_itemspecific"=IL-IS, "item_maineffects"=IL-ME, "standard"=TS
### working_cor: Working correlation for the GEE. If "userdefined" is specified (i.e. to use
###              the block unstructured or block independence correlations), zcor must
###              also be specified to define the desired structure. Specify "none"
###              if using the TS estimator.
### zcor: Design matrix for the working correlation. Will be used with the geepack geeglm
###       function. Leave as NA if independence, exchangeable, or "none" correlations are specified.
### working_cor_name: Can specify a descriptive name when using a user-defined working correlation
###                   (e.g. "Block Unstructured").
### baseline: Specify "Y" if baseline adjustment is used. The type of adjustment
###           ("baseline_adjustment") can then be specified; this simply causes the
###           output data frame to display the type of baseline adjustment used.
### baseline_adjustment: Specify "is" if estimator="item_itemspecific" or "item_maineffects".
###.                     Specify "me" if estimator="standard".
### covariates: Vector of covariates to be used in the model, including the baseline response or score.
### covariate_adjustment: Vector of the same length as "covariates" describing how they
###                       are to be adjusted. "is"=item_specific; "me"=main effects.
###                       If estimator="standard", use "me" for all adjustments.
### round_digits: Digits to round output to if format="Y".
### format: "Y" to round digits and display in a more table-friendly manner. The estimates
###         will be stored as strings in this case. "N" to leave numbers as-is.
### dist_family: mean-variance relationship family and link function for the GEE. This
###              function is set up to handle gaussian(link="identity") or binomial(link="logit").
### interaction_names: Names of the treatment*item terms in the GEE model to be combined into
###                    the total treatment effect. Usually left as NA, in which case it will
###                    sum together the treatment effects on all the items in the dataset.
### var_ratio: If TRUE, calculates the ratio of baseline response variance to follow up
###.           residual variance and returns this along with the estimates.
get_estimates <- function(data, estimator=c("item_itemspecific", "item_maineffects", "standard"),
                          working_cor=c("independence", "exchangeable", "userdefined", "none"), 
                          zcor=NA, working_cor_name=NA, baseline="Y", baseline_adjustment=NA, 
                          covariates="", covariate_adjustment="", round_digits=3,
                          format="N", dist_family=gaussian(link="identity"),
                          interaction_names=NA, var_ratio=F) {
  n_items <- length(unique(data$item))
  n_id <- length(unique(data$id))
  if (is.na(working_cor_name) & working_cor=="userdefined") {
    working_cor_name <- "User-Defined"
  }
  if (working_cor=="independence") {
    working_cor_name <- "Independence"
  }
  if (working_cor=="exchangeable") {
    working_cor_name <- "Exchangeable"
  }
  if (working_cor=="none") {
    working_cor_name <- "N/A"
  }
  
  covariates_formula <- rep("", length(covariates))
  covariates_formula[which(covariate_adjustment=="is")] <- "*factor(item)"
  covariates_formula[which(covariate_adjustment=="me")] <- ""
  
  if (baseline=="Y") {
    if (!is.na(baseline_adjustment)) {
      if (baseline_adjustment=="is") {
        baseline_name <- "Item-Specific"
      }
      if (baseline_adjustment=="me") {
        baseline_name <- "Main Effects"
      }
    } else {
      baseline_name <- "Main Effects"
    }
  }
  
  if (estimator=="item_itemspecific") {
    formula <- paste(
      gsub(" \\+ $", "", paste("response ~ trt*factor(item)",
                               paste(covariates, covariates_formula, collapse=" + ", sep=""), sep=" + ")), 
      " - 1", sep="")
    
    if (dist_family$family=="binomial") {
      formula <- gsub(pattern="^response", replacement="cbind(response, failures)", formula)
    }
    ## Formula for reduced model without baseline response as a covariate
    formula_reduced <- gsub(pattern=" \\+ baseline_response\\*factor\\(item\\)", replacement="", formula)
    
    if (working_cor=="userdefined") {
      model <- geeglm(as.formula(formula),
                      id=id, data=data, corstr=working_cor, zcor=zcor, family=dist_family)
      model_reduced <- geeglm(as.formula(formula_reduced),
                              id=id, data=data, corstr=working_cor, zcor=zcor, family=dist_family)
    } else {
      model <- geeglm(as.formula(formula),
                      id=id, data=data, corstr=working_cor, family=dist_family)
      model_reduced <- geeglm(as.formula(formula_reduced),
                              id=id, data=data, corstr=working_cor, family=dist_family)
    }
    coefs <- summary(model)$coefficients
    
    if (identical(NA, interaction_names)) {
      interaction_names <- rownames(coefs)[!is.na(str_extract(rownames(coefs), pattern="trt:factor\\(item\\)[0-9]+"))]
    }
    if (dist_family$family=="binomial" & dist_family$link=="logit") {
      K <- (data$response + data$failures)[1] # max value for each response (Bin(K, p))
      coefs <- coef(model)
      vbeta <- vcov(model)
      
      ## For each person, get linear predictor with & without treatment and find expected
      ## treatment effect on the total score
      design_matrix <- model.matrix(model)
      
      id_estimates <- lapply(unique(data$id), function(id_value) {
        data_rows <- which(data$id==id_value) 
        id_design_matrix <- design_matrix[data_rows,]
        id_data <- data %>% filter(id==id_value)
        
        ## Design matrices for this person under trt and ctrl scenarios
        id_design_matrix_trt <- id_design_matrix
        id_design_matrix_trt[,"trt"] <- 1
        id_design_matrix_trt[,interaction_names] <- diag(n_items)[,-1]
        id_design_matrix_ctrl <- id_design_matrix
        id_design_matrix_ctrl[,"trt"] <- 0
        id_design_matrix_ctrl[,interaction_names] <- 0
        
        id_data$linear_predictor_trt <- c(id_design_matrix_trt%*%coefs)
        id_data$linear_predictor_ctrl <- c(id_design_matrix_ctrl%*%coefs)
        id_data <- id_data %>% 
          mutate(
            trt_effect = expit(linear_predictor_trt) - expit(linear_predictor_ctrl)
          )
        id_trt_effect <- sum(id_data$trt_effect)
        
        ## For each row of the design matrix (each item), subtract the derivative of expit(xijTbeta) 
        ## under control from the corresponding derivative under treatment. In other words, find
        ## the derivative of the estimated linear treatment effect for this individual evaluated at beta hat.
        ## Then sum across items to get the derivative for this individual.
        ## Later we will sum across individuals, and then finally apply the delta method
        id_expit_derivative <- rowSums(sapply(1:nrow(id_design_matrix), function(i) {
          expit_derivative(x=id_design_matrix_trt[i,], beta=coefs) - expit_derivative(id_design_matrix_ctrl[i,], coefs)
        }))
        return(list(trt_effect = id_trt_effect, expit_derivative = id_expit_derivative))
      })
      estimate <- K*mean(sapply(id_estimates, "[[", "trt_effect"))
      
      ## Variance w/ delta method
      expit_derivative_sum <- rowSums(sapply(id_estimates, "[[", "expit_derivative"))
      estimate_se <- (K/n_id)*sqrt(c(t(expit_derivative_sum)%*%vbeta%*%expit_derivative_sum))
    } else {
      robust_cov_use <- vcov(model)[c("trt", interaction_names),c("trt", interaction_names)]
      if (length(interaction_names)==n_items-1) {
        estimate <- n_items*coefs["trt","Estimate"] + sum(coefs[interaction_names,"Estimate"])
        a <- c(n_items, rep(1, n_items-1))
      } else {
        estimate <- length(interaction_names)*coefs["trt","Estimate"] + sum(coefs[interaction_names,"Estimate"])
        a <- c(length(interaction_names), rep(1, length(interaction_names)))
      }
      estimate_se <- sqrt(sum_variance(robust_cov_use, a))
    }
  }
  
  if (estimator=="item_maineffects") {
    formula <- paste(
      gsub(" \\+ $", "", paste("response ~ I(trt/n_items) + factor(item)",
                               paste(covariates, covariates_formula, collapse=" + ", sep=""), 
                               sep=" + ")), 
      " - 1", sep="")
    if (dist_family$family=="binomial") {
      formula <- gsub(pattern="^response", replacement="cbind(response, failures)", formula)
    }
    formula_reduced <- gsub(pattern=" \\+ baseline_response\\*factor\\(item\\)", replacement="", formula)
    
    if (working_cor=="userdefined") {
      model <- geeglm(as.formula(formula), id=id, data=data, 
                      corstr=working_cor, zcor=zcor, family=dist_family)
      model_reduced <- geeglm(as.formula(formula_reduced), id=id, data=data, 
                              corstr=working_cor, zcor=zcor, family=dist_family)
    } else {
      model <- geeglm(as.formula(formula), id=id, data=data, 
                      corstr=working_cor, family=dist_family)
      model_reduced <- geeglm(as.formula(formula_reduced), id=id, data=data, 
                              corstr=working_cor, family=dist_family)
    }
    
    if (dist_family$family=="binomial" & dist_family$link=="logit") {
      K <- (data$response + data$failures)[1] # max value for each response (Bin(K, p))
      coefs <- coef(model)
      vbeta <- vcov(model)
      
      ## For each person, get linear predictor with & without treatment and find expected
      ## treatment effect on the total score
      design_matrix <- model.matrix(model)
      
      id_estimates <- lapply(unique(data$id), function(id_value) {
        data_rows <- which(data$id==id_value) 
        id_design_matrix <- design_matrix[data_rows,]
        id_data <- data %>% filter(id==id_value)
        
        ## Design matrices for this person under trt and ctrl scenarios
        id_design_matrix_trt <- id_design_matrix
        id_design_matrix_trt[,"I(trt/n_items)"] <- 1/n_items
        id_design_matrix_ctrl <- id_design_matrix
        id_design_matrix_ctrl[,"I(trt/n_items)"] <- 0
        
        id_data$linear_predictor_trt <- c(id_design_matrix_trt%*%coefs)
        id_data$linear_predictor_ctrl <- c(id_design_matrix_ctrl%*%coefs)
        id_data <- id_data %>% 
          mutate(
            trt_effect = expit(linear_predictor_trt) - expit(linear_predictor_ctrl)
          )
        id_trt_effect <- sum(id_data$trt_effect)
        
        ## For each row of the design matrix (each item), subtract the derivative of expit(xijTbeta) 
        ## under control from the corresponding derivative under treatment. In other words, find
        ## the derivative of the estimated linear treatment effect for this individual evaluated at beta hat.
        ## Then sum across items to get the derivative for this individual.
        ## Later we will sum across individuals, and then finally apply the delta method
        id_expit_derivative <- rowSums(sapply(1:nrow(id_design_matrix), function(i) {
          expit_derivative(x=id_design_matrix_trt[i,], beta=coefs) - expit_derivative(id_design_matrix_ctrl[i,], coefs)
        }))
        return(list(trt_effect = id_trt_effect, expit_derivative = id_expit_derivative))
      })
      estimate <- K*mean(sapply(id_estimates, "[[", "trt_effect"))
      
      ## Variance w/ delta method
      expit_derivative_sum <- rowSums(sapply(id_estimates, "[[", "expit_derivative"))
      estimate_se <- (K/n_id)*sqrt(c(t(expit_derivative_sum)%*%vbeta%*%expit_derivative_sum))
    } else {
      estimate <- summary(model)$coefficients["I(trt/n_items)","Estimate"]
      estimate_se <- summary(model)$coefficients["I(trt/n_items)","Std.err"]
    }
  }
  
  if (estimator=="standard") {
    formula <- gsub(" \\+ $", "", paste("score ~ trt",
                                        paste(covariates, covariates_formula, collapse=" + ", sep=""), 
                                        sep=" + "))
    if (dist_family$family=="binomial") {
      formula <- gsub(pattern="^score", replacement="cbind(score, failures)", formula)
    }
    formula_reduced <- gsub(" \\+ baseline_score", "", formula)
    
    model <- geeglm(as.formula(formula), id=id, data=data, family=dist_family)
    model_reduced <- geeglm(as.formula(formula_reduced), id=id, data=data, family=dist_family)
    
    if (dist_family$family=="binomial" & dist_family$link=="logit") {
      K <- (data$score + data$failures)[1] # max value for each response (Bin(K, p))
      coefs <- coef(model)
      vbeta <- vcov(model)
      
      ## For each person, get linear predictor with & without treatment and find expected
      ## treatment effect on the total score
      design_matrix <- model.matrix(model)
      
      id_estimates <- lapply(unique(data$id), function(id_value) {
        data_rows <- which(data$id==id_value) 
        id_design_matrix <- matrix(nrow=length(data_rows), data=design_matrix[data_rows,])
        colnames(id_design_matrix) <- colnames(design_matrix)
        #id_design_matrix <- design_matrix[data_rows,]
        id_data <- data %>% filter(id==id_value)
        
        ## Design matrices for this person under trt and ctrl scenarios
        id_design_matrix_trt <- id_design_matrix
        id_design_matrix_trt[,"trt"] <- 1
        id_design_matrix_ctrl <- id_design_matrix
        id_design_matrix_ctrl[,"trt"] <- 0
        
        linear_predictor_trt <- c(id_design_matrix_trt%*%coefs)
        linear_predictor_ctrl <- c(id_design_matrix_ctrl%*%coefs)
        id_trt_effect = expit(linear_predictor_trt) - expit(linear_predictor_ctrl)
        ## For each row of the design matrix (each item), subtract the derivative of expit(xijTbeta) 
        ## under control from the corresponding derivative under treatment. In other words, find
        ## the derivative of the estimated linear treatment effect for this individual evaluated at beta hat.
        ## Then sum across items to get the derivative for this individual.
        ## Later we will sum across individuals, and then finally apply the delta method
        id_expit_derivative <- rowSums(sapply(1:nrow(id_design_matrix), function(i) {
          expit_derivative(x=id_design_matrix_trt[i,], beta=coefs) - expit_derivative(id_design_matrix_ctrl[i,], coefs)
        }))
        return(list(trt_effect = id_trt_effect, expit_derivative = id_expit_derivative))
      })
      estimate <- K*mean(sapply(id_estimates, "[[", "trt_effect"))
      
      ## Variance w/ delta method
      expit_derivative_sum <- rowSums(sapply(id_estimates, "[[", "expit_derivative"))
      estimate_se <- (K/n_id)*sqrt(c(t(expit_derivative_sum)%*%vbeta%*%expit_derivative_sum))
    } else {
      estimate <- summary(model)$coefficients["trt","Estimate"]
      estimate_se <- summary(model)$coefficients["trt","Std.err"]
    }
  }

  estimate_ci <- c(estimate - qnorm(0.975)*estimate_se, estimate + qnorm(0.975)*estimate_se)
  p <- 2*pnorm(abs(estimate/estimate_se), lower.tail=F) 
  
  if (p < 0.001) {
    p_present <- "< 0.001"
  } else {
    p_present <- format(round(p, round_digits), nsmall=round_digits)
  }
  
  ## Get QIC and QICu
  qic <- format(round(unname(QIC(model)["QIC"]), 2), nsmall=2)
  qicu <- format(round(unname(QIC(model)["QICu"]), 2), nsmall=2)

  id_residuals <- get_residuals(model, estimator, data)
  residual_sd <- sd(id_residuals$residual)
  
  if (var_ratio) {
    item_residuals <- get_residuals(model, estimator, data, type="item")
    ### item level residual variance of baseline and post responses
    residual_var_data <- item_residuals %>% group_by(item) %>% summarize(residual_var=var(residual),
                                                                         baseline_var=var(baseline_response)) %>%
      ungroup() %>% data.frame() %>%
      mutate(resid_var_ratio = baseline_var/residual_var)
  }

  if (format=="Y") {
    if (baseline=="Y") {
      result <- data.frame(
        estimator=estimator,
        baseline_name=baseline_name,
        estimate=format(estimate, digits=round_digits, nsmall=round_digits),
        se=format(estimate_se, digits=round_digits, nsmall=round_digits),
        lowerci=format(estimate_ci[1], digits=round_digits, nsmall=round_digits),
        upperci=format(estimate_ci[2], digits=round_digits, nsmall=round_digits),
        p=p_present,
        working_cor=working_cor_name,
        qic=qic,
        qicu=qicu,
        residual_sd=format(residual_sd, digits=round_digits, nsmall=round_digits))
    } else {
      result <- data.frame(
        estimator=estimator,
        estimate=format(estimate, digits=round_digits, nsmall=round_digits),
        se=format(estimate_se, digits=round_digits, nsmall=round_digits),
        lowerci=format(estimate_ci[1], digits=round_digits, nsmall=round_digits),
        upperci=format(estimate_ci[2], digits=round_digits, nsmall=round_digits),
        p=p_present,
        working_cor=working_cor_name,
        qic=qic,
        qicu=qicu,
        residual_sd=format(residual_sd, digits=round_digits, nsmall=round_digits))
    }
  } else { # Don't round the digits
    if (baseline=="Y") {
      result <- data.frame(
        estimator=estimator,
        baseline_name=baseline_name,
        estimate=estimate,
        se=estimate_se,
        lowerci=estimate_ci[1],
        upperci=estimate_ci[2],
        p=p,
        working_cor=working_cor_name,
        qic=qic,
        qicu=qicu,
        residual_sd=residual_sd
      )
    } else {
      result <- data.frame(
        estimator=estimator,
        estimate=estimate,
        se=estimate_se,
        lowerci=estimate_ci[1],
        upperci=estimate_ci[2],
        p=p,
        working_cor=working_cor_name,
        qic=qic,
        qicu=qicu,
        residual_sd=residual_sd)
    }
  }
  if (var_ratio) {
    return(list(result=result, residual_var_data=residual_var_data))
  } else {
    return(result)
  }
}


### Function to produce the correlation parameter design matrix for two structures:
### corstr: Type of working correlation to produce the design matrix for.
###         "block_ind" = Block Independent; "block_unstr" = Block Unstructured
### subscale_list: List with length equal to the number of subscales. Each element
###                should be a vector with the item numbers in that subscale.
get_zcor <- function(data, corstr=c("block_ind", "block_unstr"), subscale_list) {
  zcor <- genZcor(clusz = table(data$id), waves = data$item, corstrv="unstr")
  n_subscale <- length(subscale_list)
  
  if (corstr=="block_ind") {
    n_params <- n_subscale
  }
  if (corstr=="block_unstr") {
    n_params <- n_subscale + choose(n_subscale, 2)
  }
  
  col_names_use <- vector(mode="list", length=n_params)
  zcor_result <- matrix(nrow=nrow(zcor), ncol=n_params, data=NA)
  
  if (corstr=="block_ind") {
    for (subscale in 1:n_params) {
      curr_subscale <- subscale_list[[subscale]]
      for (item in 1:(length(curr_subscale)-1) ) {
        
        col_names_use[[subscale]] <- c(col_names_use[[subscale]],
                                       paste0("alpha.", curr_subscale[item], ":", curr_subscale[(item+1):length(curr_subscale)])
        )
      }
    }
  }
  
  if (corstr=="block_unstr") {
    param <- 1
    for (subscale1 in 1:n_subscale) {
      for (subscale2 in subscale1:n_subscale) {
        curr_subscale1 <- subscale_list[[subscale1]]
        curr_subscale2 <- subscale_list[[subscale2]]
        
        if (subscale1==subscale2) {
          for (item in 1:(length(curr_subscale1)-1) ) {
            
            col_names_use[[param]] <- c(col_names_use[[param]],
                                        paste0("alpha.", curr_subscale1[item], ":",
                                               curr_subscale1[(item+1):length(curr_subscale1)])
            )
          }
        } else {
          for (item1 in 1:length(curr_subscale1)) {
            for (item2 in 1:length(curr_subscale2)) {
              if (curr_subscale1[item1] < curr_subscale2[item2]) {
                col_names_use[[param]] <- c(col_names_use[[param]],
                                            paste0("alpha.", curr_subscale1[item1], ":",
                                                   curr_subscale2[item2]))
              } else {
                col_names_use[[param]] <- c(col_names_use[[param]],
                                            paste0("alpha.", curr_subscale2[item2], ":",
                                                   curr_subscale1[item1]))
              }
            }
          }
        }
        param <- param + 1
      }
    }
  }
  for (col in 1:n_params) {
    zcor_result[,col] <- apply(zcor[,which(colnames(zcor) %in% col_names_use[[col]])], 1, sum)
  }
  
  return(zcor_result)
}


#==============================================================================#
#==============================================================================#
#===================FUNCTIONS FOR SIMULATING DATA =============================#
#==============================================================================#
#==============================================================================#

### Calculates expected value and level-specific response probabilities for an ordinal 
### variable given the set of cumulative probabilities for the possible response levels.
### If the support is not specified, it is assumed to be {0,1,...,K} where K is the
### length of the cumulative probabilities vector.

### cumulative_probs: Vector of cumulative response probabilities. Length should be
###                   one less than the number of response levels, as the last cumulative
###                   probability will always be 1.
### support: Can specify the support (e.g. c(1:5)) if it is not {0,...,K} for some K.
expected_response <- function(cumulative_probs, support=NA) {
  if (identical(support, NA)) {
    support <- c(0:length(cumulative_probs))
  }
  probs <- c(0, cumulative_probs, 1)
  response_probs <- sapply(2:length(probs), function(i) {
    probs[i] - probs[i-1]
  })
  expected_response <- 0
  for (i in 1:length(support)) {
    expected_response <- expected_response + support[i]*response_probs[i]
  }
  return(list(
    response_probs=response_probs,
    expected_response=expected_response
  ))
}


## Takes in dataset of baseline responses and regression coefficients, and one ID.
## Gets cumulative probabilities list for each person in the dataset. This has to be
## done separately for each individual because the probabilities depend on their baseline
## responses.
## alpha0 should be a vector of the length of the support of the ordinal response (e.g. 0-6),
## while alpha1 and alpha2 should be a list of length n_items containing the treatment effect
## and baseline effect coefficients (respectively) associated with each item.
get_cprob <- function(data, id_number, alpha0, alpha1, alpha2) {
  curr_data <- data %>% filter(id==id_number)
  cprob_list <- lapply(1:nrow(curr_data), function(i) {
    expit(alpha0 + as.numeric(curr_data[i,]["trt"])*alpha1[[as.numeric(curr_data[i,]["item"])]] +
            as.numeric(curr_data[i,]["baseline_response"])*alpha2[[as.numeric(curr_data[i,]["item"])]])
  })
  return(cprob_list)
}


## Function to simulate a single correlated ordinal dataset using a cumulative logit
## model (CLM) given:
## n: sample size
## n_items: number of items in the questionnaire
## alpha0: response specific intercepts for CLM
## alpha1: treatment effect coefficients for CLM (list of length n_items)
## alpha2: baseline effect coefficients for CLM
## support_list: list of length n_items with the support for each item, e.g. c(0:6)
## within_corr: correlation between items within a subscale
## without_corr: correlation between items from different subscales
## subscale_sizes: vector with the number of items per subscale. Should sum to n_items.
simulate_dataset <- function(n, n_items, alpha0, alpha2, support_list, treatment_effects,
                             within_corr, without_corr, subscale_sizes) {
  ## Correlation matrix between items
  ## Set items within a subscale to have a correlation of within_corr and
  ## items from different subscales to have a correlation of without_corr
  Sigma <- matrix(without_corr, nrow=n_items, ncol=n_items)
  for (i in 1:length(subscale_sizes)) {
    if (i==1) {
      Sigma[1:subscale_sizes[i],1:subscale_sizes[i]] <- matrix(data=within_corr,
                                                               nrow=subscale_sizes[i],
                                                               ncol=subscale_sizes[i])
    } else {
      Sigma[(cumsum(subscale_sizes)[i-1]+1):cumsum(subscale_sizes)[i],
            (cumsum(subscale_sizes)[i-1]+1):cumsum(subscale_sizes)[i]] <- matrix(data=within_corr,
                                                                                 nrow=subscale_sizes[i],
                                                                                 ncol=subscale_sizes[i])
    }
  }
  diag(Sigma) <- 1
  
  ## Generate the baseline responses
  prob_list_baseline <- lapply(1:n_items, function(i) return(expit(alpha0)))
  sim_baseline_sample <- ordsample(n, marginal=prob_list_baseline, support=support_list,
                                   Sigma=Sigma)
  sim_baseline <- sim_baseline_sample %>%
    data.frame() %>%
    rename_with(~gsub("^X", "item", .x)) %>%
    mutate(
      id=row_number(),
      baseline_score = rowSums(across(starts_with("item")))
    )
  sim_baseline$trt <- rbinom(n, 1, 0.5) # generate 1:1 treatment assignments
  sim_baseline_long <- sim_baseline %>%
    pivot_longer(
      cols=starts_with("item"),
      names_to="item",
      values_to="baseline_response",
      names_pattern="item([0-9]+)"
    ) %>%
    mutate(
      item=as.numeric(item)
    ) %>%
    arrange(id, item)
  
  mean_baseline_responses <- (sim_baseline_long %>%
                                group_by(item) %>%
                                summarize(mean=mean(baseline_response)) %>%
                                data.frame())$mean
  
  # Given the item specific intercepts (alpha0) and baseline response coefficients,
  # find the treatment effect coefficients corresponding to the desired treatment
  # effect when the baseline response is equal to the observed mean baseline response
  # in the data.
  
  ## Value of alpha1 that makes difference between the expected response with
  ## treatment and the expected response without treatment equal to the target
  ## treatment effect
  treatment_effects_used <- treatment_effects
  alpha1_info <- vector(mode="list", length=length(treatment_effects))
  for (i in 1:length(alpha1_info)) {
    # Expected response without treatment
    expected_response_ctrl <- expected_response(expit(alpha0 + lapply(1:length(alpha2), function(j)
      alpha2[[j]]*mean_baseline_responses[j])[[i]]), support_list[[i]])$expected_response
    
    ## If the expected response without treatment is too high, such that the expected response
    ## with treatment would end up larger than the max value of the support (i.e. > 6), 
    ## uniroot will be unable to find the root and we will get an error. So I will
    ## program it such that in this scenario, it will find the coefficients corresponding
    ## to an expected response of 5.99 (close to the max) with treatment.
    if (expected_response_ctrl + treatment_effects[i] > (max(support_list[[i]])-0.01) ) {
      trt_effect_use <- max(support_list[[i]]) - 0.01 - expected_response_ctrl
      alpha1_info[[i]]$trt_effect_adjusted <- T
      treatment_effects_used[i] <- trt_effect_use
      ## Similarly, if the expected response without treatment is too low and the treatment
      ## effect is negative, we need to adjust it so that the expected response with 
      ## treatment won't be below 0 (or the minimum value of the support).
    } else if (expected_response_ctrl + treatment_effects[i] < (min(support_list[[i]])+0.01) ) {
      trt_effect_use <- min(support_list[[i]]) + 0.01 - expected_response_ctrl
      alpha1_info[[i]]$trt_effect_adjusted <- T
      treatment_effects_used[i] <- trt_effect_use
    } else {
      trt_effect_use <- treatment_effects[i]
      alpha1_info[[i]]$trt_effect_adjusted <- F
    }
    
    alpha1_info[[i]]$alpha1 <- rep(uniroot(function(alpha0, alpha1, alpha2) {
      # Expected response with treatment
      expected_response(expit(alpha0 + alpha1 + lapply(1:length(alpha2), function(j)
        alpha2[[j]]*mean_baseline_responses[j])[[i]]), support_list[[i]])$expected_response -
        
        (expected_response_ctrl + trt_effect_use)
    }, c(-10, 10), alpha0=alpha0, alpha2=alpha2)$root,
    length(alpha0))
  }
  alpha1 <- lapply(alpha1_info, "[[", "alpha1")
  adjustment_indicator <- sapply(alpha1_info, "[[", "trt_effect_adjusted")
  
  ## Use the baseline responses to generate the cumulative prob vector for each
  ## individual in the dataset
  all_cprobs <- lapply(1:length(unique(sim_baseline$id)), function(i) {
    get_cprob(sim_baseline_long, i, alpha0, alpha1, alpha2)
  })
  
  sim_sample_list <- lapply(1:length(all_cprobs), function(i) {
    samples <- ordsample(1, marginal=all_cprobs[[i]], support=support_list,
                         Sigma=Sigma, cormat="continuous")
    return(list(samples=samples))
  })
  
  sim_sample <- do.call(rbind, lapply(sim_sample_list, "[[", "samples"))
  sim_data_wide <- sim_sample %>%
    data.frame() %>%
    rename_with(~gsub("^X", "item", .x)) %>%
    mutate(
      id=row_number(),
      score = rowSums(across(starts_with("item")))
    ) %>%
    left_join(sim_baseline %>% dplyr::select(id, trt, baseline_score), by="id")
  
  sim_data_long <- sim_data_wide %>%
    pivot_longer(
      cols=starts_with("item"),
      names_to="item",
      values_to="response",
      names_pattern="item([0-9]+)"
    ) %>%
    mutate(item=as.numeric(item)) %>%
    left_join(sim_baseline_long %>% dplyr::select(id, item, baseline_response), by=c("id", "item")) %>%
    dplyr::select(id, trt, item, response, baseline_response, score, baseline_score) %>%
    arrange(id, item)
  
  return(list(
    sim_data_wide=sim_data_wide,
    sim_data_long=sim_data_long,
    treatment_effects_used=treatment_effects_used,
    treatment_effects_adjusted=adjustment_indicator
  ))
}


## Simulate a multivariate normal dataset of item responses, given:
## n: Sample size
## n_items: Number of items
## beta0: Intercept. Should be a scalar; this function uses the same intercept for each item.
## baseline_effects: Baseline response coefficients. Length equal to n_items.
## treatment_effects: Treatment effect coefficients. Length equal to n_items.
## sigma2_baseline: Marginal variance of baseline item responses
## sigma2_post: Marginal variance of follow up item responses
## within_corr: Correlation between items in the same subscale
## without_corr: Correlation between items from different subscales
## subscale_sizes: Vector with length equal to the number of subscales, specifying the
##                 number of items within each subscale. The numbers should add up to n_items.
simulate_dataset_mvn <- function(n, n_items, beta0, baseline_effects, treatment_effects, 
                                 sigma2_baseline=1, sigma2_post=1,
                                 within_corr, without_corr, subscale_sizes) {
  ## Correlation matrix between items
  ## Set items within a subscale to have a correlation of within_corr and
  ## items from different subscales to have a correlation of without_corr
  Sigma <- matrix(without_corr, nrow=n_items, ncol=n_items)
  for (i in 1:length(subscale_sizes)) {
    if (i==1) {
      Sigma[1:subscale_sizes[i],1:subscale_sizes[i]] <- matrix(data=within_corr,
                                                               nrow=subscale_sizes[i],
                                                               ncol=subscale_sizes[i])
    } else {
      Sigma[(cumsum(subscale_sizes)[i-1]+1):cumsum(subscale_sizes)[i],
            (cumsum(subscale_sizes)[i-1]+1):cumsum(subscale_sizes)[i]] <- matrix(data=within_corr,
                                                                                 nrow=subscale_sizes[i],
                                                                                 ncol=subscale_sizes[i])
    }
  }
  diag(Sigma) <- 1
  Sigma_baseline <- sigma2_baseline*Sigma
  Sigma_post <- sigma2_post*Sigma
  
  baseline_mu <- rep(beta0, n_items)
  
  ## Generate the baseline responses
  sim_baseline_sample <- mvrnorm(n=n, mu=baseline_mu, Sigma=Sigma_baseline)
  
  sim_baseline <- sim_baseline_sample %>%
    data.frame() %>%
    rename_with(~gsub("^X", "item", .x)) %>%
    mutate(
      id=row_number(),
      baseline_score = rowSums(across(starts_with("item")))
    )
  sim_baseline$trt <- rbinom(n, 1, 0.5) # generate 1:1 treatment assignments
  sim_baseline_long <- sim_baseline %>%
    pivot_longer(
      cols=starts_with("item"),
      names_to="item",
      values_to="baseline_response",
      names_pattern="item([0-9]+)"
    ) %>%
    mutate(
      item=as.numeric(item)
    ) %>%
    arrange(id, item)
  
  effects_map <- data.frame(
    item=1:length(treatment_effects),
    beta0=beta0,
    beta1=treatment_effects,
    beta2=baseline_effects
  )
  
  means_data <- sim_baseline_long %>% 
    left_join(effects_map, by="item") %>%
    mutate(
      mu = beta0 + beta1*trt + beta2*baseline_response
    )
  
  sim_data_long <- do.call(rbind, lapply(unique(sim_baseline$id), function(i) {
    curr_data <- means_data %>% filter(id==i)
    curr_data$response <- mvrnorm(n=1, mu=curr_data$mu, Sigma=Sigma_post)
    return(curr_data)
  })) %>%
    dplyr::select(id, baseline_score, trt, item, baseline_response, response)
  
  sim_data_wide <- sim_data_long %>%
    group_by(id, trt) %>%
    summarize(score=sum(response), baseline_score=sum(baseline_response)) %>%
    ungroup() %>%
    data.frame()
  
  return(list(
    sim_data_wide=sim_data_wide,
    sim_data_long=sim_data_long
  ))
}
