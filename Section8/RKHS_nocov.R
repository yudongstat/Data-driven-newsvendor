rm(list = ls())
library(parallel)
simulation <- function(cut_off){
  my_simulation_procedure <- function(cut_off){
    library(splines2)
    library(pracma)
    library(Rsolnp)
    chosen_index <- 2
    load("bakery_data.rdata")
    T <- 6
    h_cost <- rep(1, T)
    b_cost <- rep(1, T)
    n <- nrow(my_data)/T
    
    full_D_mat <- matrix(my_data[[products_all[indexing_order[chosen_index]]]],
                         n, T, byrow = TRUE)
    
    full_X_mat <- matrix(NA, n, T)
    
    for(t in 2:T){
      full_X_mat[, t] <- full_D_mat[, t-1] 
    }
    
    start_time <- proc.time()
    sigma_candidate <- seq(50, 100, 10)#(seq(10, 50, length.out = 5))^2
    lambda_candidate <- 5^(-seq(1, 10, by = 2))
    
    
    gaussian_kernel<-function(x1,x2,sigma=1){
      return(exp(-(x1-x2)^2/sigma^2))
    }
    
    
    ind_selected <- rep(NA, 2)
    CV_time <- NA
    run_time <- NA
    
    if(cut_off == 0){
      test_decision_our <- rep(100, T)
      test_D_mat <- full_D_mat[cut_off + 1, ]
    }
    
    if(cut_off == 1){
      test_decision_our <- full_D_mat[1, ]
      test_D_mat <- full_D_mat[cut_off + 1, ]
    }
    
    if(cut_off > 1){
      D_mat <- full_D_mat[1:cut_off, ]
      X_mat <- full_X_mat[1:cut_off, ]
      if(cut_off < 5){
        sigma <- sample(x = sigma_candidate, size = 1)
        mylambda <- sample(x = lambda_candidate, size = 1)
      }else{
        ########################################################################################################################################################################################################################################
        ########################################################################################################################################################################################################################################
        CV_array <- array(0, dim = c(length(sigma_candidate), length(lambda_candidate), 5))
        my_split<-1:cut_off
        if(cut_off%%5 == 0){
          my_split<-matrix(my_split,5, cut_off/5, byrow = TRUE)
        }else{
          my_split <- c(my_split, rep(1e4, (5*ceiling(cut_off/5) - cut_off)))
          my_split <- matrix(my_split, 5, ceiling(cut_off/5))
        }
        
        
        for(fold in 1:5){
          test_index <- c(my_split[fold,])
          test_index <- test_index[which(test_index < 5e3)]
          test_D_mat<- D_mat[test_index, ]
          test_X_mat<- X_mat[test_index, ]
          
          train_index <- c(my_split[-fold,])
          train_index <- train_index[which(train_index < 5e3)]
          train_D_mat <- D_mat[train_index, ]
          train_X_mat <- X_mat[train_index, ]
          
          train_n <- nrow(train_D_mat)
          test_n <- nrow(test_D_mat)
          if(is.null(test_n)){
            test_n <- 1
          }
          
          for(ppp in 1:length(sigma_candidate)){
            for(qqq in 1:length(lambda_candidate)){
              k_mat <- array(NA, dim = c(train_n, train_n, T))
              for(t in 2:T){
                feature_mat <- train_X_mat[, t]
                for(i in 1:train_n){
                  for(j in 1:train_n){
                    k_mat[i, j, t] <- gaussian_kernel(x1 = feature_mat[i],
                                                      x2 = feature_mat[j],
                                                      sigma = sigma_candidate[ppp])
                  }
                }
              }
              
              keropt_cost <- function(a_vec, mylambda = NA){
                period1_constants <- a_vec[1]
                starting_level <- period1_constants
                total_cost <- sum(h_cost[1] * pmax(0, starting_level - train_D_mat[, 1]) + b_cost[1] * pmax(0, train_D_mat[, 1] - starting_level)) 
                ending_level <- pmax(0, starting_level - train_D_mat[, 1])
                
                a_mat <- matrix(a_vec[-1], T - 1, train_n)
                for(t in 2:T){
                  order_up_to <- c(k_mat[, , t] %*% a_mat[t-1, ])
                  if(min(order_up_to) < 0){
                    total_cost <- total_cost + 1e4
                  }
                  starting_level <- pmax(ending_level, order_up_to)
                  total_cost <- total_cost + sum(h_cost[t] * pmax(0, starting_level - train_D_mat[, t]) + b_cost[t] * pmax(0, train_D_mat[, t] - starting_level))
                  ending_level <- pmax(0, starting_level - train_D_mat[, t])
                  total_cost <- total_cost + c(train_n * mylambda * t(a_mat[t-1, ] ) %*% k_mat[, , t] %*% a_mat[t-1, ])
                }
                
                return(total_cost)
              }
              
              opt_res<-solnp(pars = rep(0.1, 1 + train_n*(T - 1)),
                             fun = keropt_cost,
                             mylambda = lambda_candidate[qqq],
                             control = list(trace = 0))
              a_vec <- opt_res$pars
              period1_constants <- a_vec[1]
              a_mat <- matrix(a_vec[-1], T - 1, train_n)
              
              if(test_n > 1){
                k_mat_test <- array(NA, dim = c(test_n, train_n, T))
                for(t in 2:T){
                  feature_mat_train <- train_X_mat[, t]
                  feature_mat_test <- test_X_mat[, t]
                  for(i in 1:test_n){
                    for(j in 1:train_n){
                      k_mat_test[i, j, t] <- gaussian_kernel(x1 = feature_mat_test[i],
                                                             x2 = feature_mat_train[j],
                                                             sigma = sigma_candidate[ppp])
                    }
                  }
                }
                
                starting_level <- period1_constants
                total_cost <- sum(h_cost[1] * pmax(0, starting_level - test_D_mat[, 1]) + b_cost[1] * pmax(0, test_D_mat[, 1] - starting_level)) 
                ending_level <- pmax(0, starting_level - test_D_mat[, 1])
                
                for(t in 2:T){
                  order_up_to <- c(k_mat_test[, , t] %*% a_mat[t-1, ])
                  starting_level <- pmax(ending_level, order_up_to)
                  total_cost <- total_cost + sum(h_cost[t] * pmax(0, starting_level - test_D_mat[, t]) + b_cost[t] * pmax(0, test_D_mat[, t] - starting_level))
                  ending_level <- pmax(0, starting_level - test_D_mat[, t])
                }
              }else{
                k_mat_test <- matrix(NA, train_n, T)
                for(t in 2:T){
                  feature_mat_train <- train_X_mat[, t]
                  feature_mat_test <- test_X_mat[t]
                  for(j in 1:train_n){
                    k_mat_test[j, t] <- gaussian_kernel(x1 = feature_mat_test,
                                                        x2 = feature_mat_train[j],
                                                        sigma = sigma_candidate[ppp])
                  }
                }
                
                starting_level <- period1_constants
                total_cost <- h_cost[1] * pmax(0, starting_level - test_D_mat[1]) + b_cost[1] * pmax(0, test_D_mat[1] - starting_level) 
                ending_level <- pmax(0, starting_level - test_D_mat[1])
                
                for(t in 2:T){
                  order_up_to <- sum(k_mat_test[, t] * a_mat[t-1, ])
                  starting_level <- pmax(ending_level, order_up_to)
                  total_cost <- total_cost + h_cost[t] * pmax(0, starting_level - test_D_mat[t]) + b_cost[t] * pmax(0, test_D_mat[t] - starting_level)
                  ending_level <- pmax(0, starting_level - test_D_mat[t])
                }
              }
              
              
              CV_array[ppp, qqq, fold] <- total_cost
              print(c(fold, ppp, qqq, opt_res$convergence))
            }
          }
        }
        CV_mat <- apply(CV_array, MARGIN = c(1, 2), sum)
        ind_selected <- which(CV_mat == min(CV_mat), arr.ind = TRUE)[1,]
        sigma <- sigma_candidate[ind_selected[1]] 
        mylambda <- lambda_candidate[ind_selected[2]]
        
        CV_time <- (proc.time() - start_time)["elapsed"]
      }
      
      start_time <- proc.time()
      train_D_mat <- D_mat[1:cut_off, ]
      train_X_mat <- X_mat[1:cut_off, ]
      
      test_D_mat <- full_D_mat[cut_off + 1, ]
      test_X_mat <- full_X_mat[cut_off + 1, ]
      
      train_n <- nrow(train_D_mat)
      test_n <- 1
      
      k_mat <- array(NA, dim = c(train_n, train_n, T))
      for(t in 2:T){
        feature_mat <- train_X_mat[, t]
        for(i in 1:train_n){
          for(j in 1:train_n){
            k_mat[i, j, t] <- gaussian_kernel(x1 = feature_mat[i],
                                              x2 = feature_mat[j],
                                              sigma = sigma)
          }
        }
      }
      
      keropt_cost <- function(a_vec, mylambda = NA){
        period1_constants <- a_vec[1]
        starting_level <- period1_constants
        total_cost <- sum(h_cost[1] * pmax(0, starting_level - train_D_mat[, 1]) + b_cost[1] * pmax(0, train_D_mat[, 1] - starting_level)) 
        ending_level <- pmax(0, starting_level - train_D_mat[, 1])
        
        a_mat <- matrix(a_vec[-1], T - 1, train_n)
        for(t in 2:T){
          order_up_to <- c(k_mat[, , t] %*% a_mat[t-1, ])
          if(min(order_up_to) < 0){
            total_cost <- total_cost + 1e4
          }
          starting_level <- pmax(ending_level, order_up_to)
          total_cost <- total_cost + sum(h_cost[t] * pmax(0, starting_level - train_D_mat[, t]) + b_cost[t] * pmax(0, train_D_mat[, t] - starting_level))
          ending_level <- pmax(0, starting_level - train_D_mat[, t])
          total_cost <- total_cost + c(train_n * mylambda * t(a_mat[t-1, ] ) %*% k_mat[, , t] %*% a_mat[t-1, ])
        }
        
        return(total_cost)
      }
      
      opt_res<-solnp(pars = rep(0.1, 1 + train_n*(T - 1)),
                     fun = keropt_cost,
                     mylambda = mylambda,
                     control = list(trace = 0))
      a_vec <- opt_res$pars
      period1_constants <- a_vec[1]
      a_mat <- matrix(a_vec[-1], T - 1, train_n)
      
      k_mat_test <- matrix(NA, train_n, T)
      for(t in 2:T){
        feature_mat_train <- train_X_mat[, t]
        feature_mat_test <- test_X_mat[t]
        for(j in 1:train_n){
          k_mat_test[j, t] <- gaussian_kernel(x1 = feature_mat_test,
                                              x2 = feature_mat_train[j],
                                              sigma = sigma)
        }
      }
      
      test_decision_our <- period1_constants
      
      for(t in 2:T){
        order_up_to <- sum(k_mat_test[, t] * a_mat[t-1, ])
        test_decision_our <- c(test_decision_our, order_up_to)
      }
      
      run_time <- (proc.time() - start_time)["elapsed"]
    }
    
    test_cost_our <- 0
    temp_level_our <- test_decision_our[1]
    for(qqqq in 1:T){
      test_cost_our<-test_cost_our+(h_cost[qqqq]*pmax(0,temp_level_our-test_D_mat[qqqq])+b_cost[qqqq]*pmax(0,test_D_mat[qqqq]-temp_level_our))
      
      if(qqqq<T){
        temp_level_our<-pmax(temp_level_our-test_D_mat[qqqq],test_decision_our[qqqq+1])
      }
    }
    
    
    return(c(cut_off, test_cost_our, CV_time, run_time, ind_selected))
  }
  res <- try(my_simulation_procedure(cut_off),silent = TRUE)
  return(res)
}

clnum <- 28
cl <- makeCluster(getOption("cl.cores", clnum))
res <- parLapply(cl, 0:83,  simulation)
stopCluster(cl)

save(res, file = "RKHS_nocov.rdata")
