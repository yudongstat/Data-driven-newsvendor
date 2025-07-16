rm(list = ls())
library(parallel)
simulation <- function(cut_off){
  library(R.utils)
  my_simulation_procedure<-function(cut_off){
    library(splines2)
    library(pracma)
    library(Rsolnp)
    chosen_index <- 2
    multiplier1 <- 2.5
    multiplier2 <- 0.2
    load("bakery_data.rdata")
    m_n_candidate <- 2:10
    
    
    
    T <- 6
    h_cost <- rep(1, T)
    b_cost <- rep(1, T)
    n <- nrow(my_data)/T
    
    full_D_mat <- matrix(my_data[[products_all[indexing_order[chosen_index]]]],
                         n, T, byrow = TRUE)
    
    
    
    m_n <- NA
    CV_time_our <- NA
    run_time_our <- NA
    
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
      
      lower_knots <- rep(NA, T)
      upper_knots <- rep(NA, T)
      
      for(t in 2:T){
        upper_knots[t] <- max(D_mat[, t-1]) 
        lower_knots[t] <- min(D_mat[, t-1]) 
        upper_knots[t] <- upper_knots[t] * multiplier1
        if(lower_knots[t] > 0){
          lower_knots[t] <- lower_knots[t] * multiplier2
        }else{
          lower_knots[t] <- lower_knots[t] * multiplier1
        }
      }
      if(cut_off < 5){
        m_n <- 1
      }else{
        start_time <- proc.time()
        ########################################################################################################################################################################################################################################
        ########################################################################################################################################################################################################################################
        CV_mat <- matrix(0, 5, length(m_n_candidate))
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
          
          train_index <- c(my_split[-fold,])
          train_index <- train_index[which(train_index < 5e3)]
          train_D_mat <- D_mat[train_index, ]
          
          for(ppp in 1:length(m_n_candidate)){
            m_n <- m_n_candidate[ppp]
            
            cut_quantile <- seq(from = 0,to = 1,length.out = m_n+2)[-c(1, m_n+2)]
            
            knots_6 <- quantile(x = train_D_mat[, 5], probs = cut_quantile)
            
            price_6 <- function(alpha){
              basis_mat_6 <- bSpline(x = train_D_mat[,5],
                                     knots = knots_6, degree = 3, 
                                     Boundary.knots = c(lower_knots[6], upper_knots[6]),
                                     intercept = TRUE)
              temp_level <- c(basis_mat_6 %*% alpha)
              return(sum(h_cost[6] * pmax(0, temp_level - train_D_mat[, 6]) + b_cost[6] * pmax(0, train_D_mat[, 6] - temp_level)))
            }
            
            price_6_empirical <- function(temp_level){
              return(sum(h_cost[6] * pmax(0, temp_level - train_D_mat[, 6]) + b_cost[6] * pmax(0, train_D_mat[, 6] - temp_level)))
            }
            
            
            res_6 <- solnp(pars = rep(1, m_n+4),
                           fun = price_6,
                           LB = rep(0, m_n+4),
                           UB = rep(Inf,m_n+4),
                           control = list(trace = 0))
            
            basis_mat_6 <- bSpline(x = train_D_mat[, 5],
                                   knots = knots_6, degree = 3,
                                   Boundary.knots = c(lower_knots[6], upper_knots[6]),
                                   intercept = TRUE)
            best_decision_6 <- c(basis_mat_6 %*% res_6$par)
            
            knots_5 <- quantile(x = train_D_mat[, 4], probs = cut_quantile)
            
            price_5 <- function(alpha){
              basis_mat_5 <- bSpline(x = train_D_mat[,4],
                                     knots = knots_5, degree = 3, 
                                     Boundary.knots = c(lower_knots[5], upper_knots[5]),
                                     intercept = TRUE)
              temp_level <- c(basis_mat_5 %*% alpha)
              return(sum(h_cost[5] * pmax(0, temp_level - train_D_mat[, 5]) + b_cost[5] * pmax(0, train_D_mat[, 5] - temp_level)) + price_6_empirical(pmax(temp_level - train_D_mat[,5], best_decision_6)))
            }
            
            price_5_empirical <- function(temp_level){
              return(sum(h_cost[5] * pmax(0, temp_level - train_D_mat[, 5]) + b_cost[5] * pmax(0, train_D_mat[, 5] - temp_level)) + price_6_empirical(pmax(temp_level - train_D_mat[,5], best_decision_6)))
            }
            
            
            res_5 <- solnp(pars = rep(1, m_n+4),
                           fun = price_5,
                           LB = rep(0, m_n+4),
                           UB = rep(Inf,m_n+4),
                           control = list(trace = 0))
            
            basis_mat_5 <- bSpline(x = train_D_mat[, 4],
                                   knots = knots_5, degree = 3,
                                   Boundary.knots = c(lower_knots[5], upper_knots[5]),
                                   intercept = TRUE)
            best_decision_5 <- c(basis_mat_5 %*% res_5$par)
            
            
            knots_4 <- quantile(x = train_D_mat[, 3], probs = cut_quantile)
            
            
            price_4 <- function(alpha){
              basis_mat_4 <- bSpline(x = train_D_mat[,3],
                                     knots = knots_4, degree = 3, 
                                     Boundary.knots = c(lower_knots[4], upper_knots[4]),
                                     intercept = TRUE)
              temp_level <- c(basis_mat_4 %*% alpha)
              return(sum(h_cost[4] * pmax(0, temp_level - train_D_mat[, 4]) + b_cost[4] * pmax(0, train_D_mat[, 4] - temp_level)) + price_5_empirical(pmax(temp_level - train_D_mat[,4], best_decision_5)))
            }
            
            price_4_empirical <- function(temp_level){
              return(sum(h_cost[4] * pmax(0, temp_level - train_D_mat[, 4]) + b_cost[4] * pmax(0, train_D_mat[, 4] - temp_level)) + price_5_empirical(pmax(temp_level - train_D_mat[,4], best_decision_5)))
            }
            
            
            res_4 <- solnp(pars = rep(1, m_n+4),
                           fun = price_4,
                           LB = rep(0, m_n+4),
                           UB = rep(Inf,m_n+4),
                           control = list(trace = 0))
            
            basis_mat_4 <- bSpline(x = train_D_mat[, 3],
                                   knots = knots_4, degree = 3,
                                   Boundary.knots = c(lower_knots[4], upper_knots[4]),
                                   intercept = TRUE)
            best_decision_4 <- c(basis_mat_4 %*% res_4$par)
            
            knots_3 <- quantile(x = train_D_mat[, 2], probs = cut_quantile)
            
            
            price_3 <- function(alpha){
              basis_mat_3 <- bSpline(x = train_D_mat[,2],
                                     knots = knots_3, degree = 3, 
                                     Boundary.knots = c(lower_knots[3], upper_knots[3]),
                                     intercept = TRUE)
              temp_level <- c(basis_mat_3 %*% alpha)
              return(sum(h_cost[3] * pmax(0, temp_level - train_D_mat[, 3]) + b_cost[3] * pmax(0, train_D_mat[, 3] - temp_level)) + price_4_empirical(pmax(temp_level - train_D_mat[,3], best_decision_4)))
            }
            
            price_3_empirical <- function(temp_level){
              return(sum(h_cost[3] * pmax(0, temp_level - train_D_mat[, 3]) + b_cost[3] * pmax(0, train_D_mat[, 3] - temp_level)) + price_4_empirical(pmax(temp_level - train_D_mat[,3], best_decision_4)))
            }
            
            
            res_3 <- solnp(pars = rep(1, m_n+4),
                           fun = price_3,
                           LB = rep(0, m_n+4),
                           UB = rep(Inf,m_n+4),
                           control = list(trace = 0))
            
            basis_mat_3 <- bSpline(x = train_D_mat[, 2],
                                   knots = knots_3, degree = 3,
                                   Boundary.knots = c(lower_knots[3], upper_knots[3]),
                                   intercept = TRUE)
            best_decision_3 <- c(basis_mat_3 %*% res_3$par)
            
            
            knots_2 <- quantile(x = train_D_mat[, 1], probs = cut_quantile)
            
            
            price_2 <- function(alpha){
              basis_mat_2 <- bSpline(x = train_D_mat[,1],
                                     knots = knots_2, degree = 3, 
                                     Boundary.knots = c(lower_knots[2], upper_knots[2]),
                                     intercept = TRUE)
              temp_level <- c(basis_mat_2 %*% alpha)
              return(sum(h_cost[2] * pmax(0, temp_level - train_D_mat[, 2]) + b_cost[2] * pmax(0, train_D_mat[, 2] - temp_level)) + price_3_empirical(pmax(temp_level - train_D_mat[,2], best_decision_3)))
            }
            
            price_2_empirical <- function(temp_level){
              return(sum(h_cost[2] * pmax(0, temp_level - train_D_mat[, 2]) + b_cost[2] * pmax(0, train_D_mat[, 2] - temp_level)) + price_3_empirical(pmax(temp_level - train_D_mat[,2], best_decision_3)))
            }
            
            
            res_2 <- solnp(pars = rep(1, m_n+4),
                           fun = price_2,
                           LB = rep(0, m_n+4),
                           UB = rep(Inf,m_n+4),
                           control = list(trace = 0))
            
            basis_mat_2 <- bSpline(x = train_D_mat[, 1],
                                   knots = knots_2, degree = 3,
                                   Boundary.knots = c(lower_knots[2], upper_knots[2]),
                                   intercept = TRUE)
            best_decision_2 <- c(basis_mat_2 %*% res_2$par)
            
            price_1 <- function(alpha){
              return(sum(h_cost[1] * pmax(0, alpha - train_D_mat[, 1]) + b_cost[1] * pmax(0, train_D_mat[, 1] - alpha)) + price_2_empirical(pmax(alpha - train_D_mat[,1], best_decision_2)))
            }
            
            
            res_1 <- optimize(f = price_1, interval = c(0, 1000))
            
            best_decision_1 <- res_1$minimum
            
            test_decision_1 <- best_decision_1
            
            if(is.null(nrow(test_D_mat))){
              temp_x <- test_D_mat[1]
              temp_base <- bSpline(x = temp_x, knots = knots_2, degree = 3,
                                   Boundary.knots = c(lower_knots[2], upper_knots[2]),
                                   intercept = TRUE)
              test_decision_2 <- c(temp_base%*%res_2$par)
              
              temp_x <- test_D_mat[2]
              temp_base <- bSpline(x = temp_x, knots = knots_3, degree = 3,
                                   Boundary.knots = c(lower_knots[3], upper_knots[3]),
                                   intercept = TRUE)
              test_decision_3 <- c(temp_base%*%res_3$par)
              
              
              temp_x <- test_D_mat[3]
              temp_base <- bSpline(x = temp_x, knots = knots_4, degree = 3,
                                   Boundary.knots = c(lower_knots[4], upper_knots[4]),
                                   intercept = TRUE)
              test_decision_4 <- c(temp_base%*%res_4$par)
              
              
              temp_x <- test_D_mat[4]
              temp_base <- bSpline(x = temp_x, knots = knots_5, degree = 3,
                                   Boundary.knots = c(lower_knots[5], upper_knots[5]),
                                   intercept = TRUE)
              test_decision_5 <- c(temp_base%*%res_5$par)
              
              
              temp_x <- test_D_mat[5]
              temp_base <- bSpline(x = temp_x, knots = knots_6, degree = 3,
                                   Boundary.knots = c(lower_knots[6], upper_knots[6]),
                                   intercept = TRUE)
              test_decision_6 <- c(temp_base%*%res_6$par)
              
              test_decision <- c(test_decision_2, test_decision_3, test_decision_4, test_decision_5, test_decision_6)
              
              temp_res<-0
              temp_level <- test_decision_1
              for(qqqq in 1:T){
                temp_res <- temp_res + sum(h_cost[qqqq] * pmax(0, temp_level - test_D_mat[qqqq]) + b_cost[qqqq] * pmax(0, test_D_mat[qqqq] - temp_level))
                if(qqqq<T){
                  temp_level <- pmax(temp_level - test_D_mat[qqqq], test_decision[qqqq]) 
                }
              }
              CV_mat[fold,ppp] <- temp_res
              
            }else{
              temp_x <- test_D_mat[, 1]
              temp_base <- bSpline(x = temp_x,
                                   knots = knots_2,
                                   degree = 3,
                                   Boundary.knots = c(lower_knots[2], upper_knots[2]),
                                   intercept = TRUE)
              test_decision_2 <- c(temp_base%*%res_2$par)
              
              temp_x <- test_D_mat[, 2]
              temp_base <- bSpline(x = temp_x,
                                   knots = knots_3,
                                   degree = 3,
                                   Boundary.knots = c(lower_knots[3], upper_knots[3]),
                                   intercept = TRUE)
              test_decision_3 <- c(temp_base%*%res_3$par)
              
              
              temp_x <- test_D_mat[, 3]
              temp_base <- bSpline(x = temp_x,
                                   knots = knots_4,
                                   degree = 3,
                                   Boundary.knots = c(lower_knots[4], upper_knots[4]),
                                   intercept = TRUE)
              test_decision_4 <- c(temp_base%*%res_4$par)
              
              temp_x <- test_D_mat[, 4]
              temp_base <- bSpline(x = temp_x,
                                   knots = knots_5,
                                   degree = 3,
                                   Boundary.knots = c(lower_knots[5], upper_knots[5]),
                                   intercept = TRUE)
              test_decision_5 <- c(temp_base%*%res_5$par)
              
              temp_x <- test_D_mat[, 5]
              temp_base <- bSpline(x = temp_x,
                                   knots = knots_6,
                                   degree = 3,
                                   Boundary.knots = c(lower_knots[6], upper_knots[6]),
                                   intercept = TRUE)
              test_decision_6 <- c(temp_base%*%res_6$par)
              
              
              test_decision <- rbind(test_decision_2, test_decision_3, test_decision_4, test_decision_5, test_decision_6)
              
              temp_res<-0
              temp_level <- test_decision_1
              for(qqqq in 1:T){
                temp_res <- temp_res + sum(h_cost[qqqq] * pmax(0, temp_level - test_D_mat[,qqqq]) + b_cost[qqqq] * pmax(0, test_D_mat[,qqqq] - temp_level))
                if(qqqq<T){
                  temp_level<-pmax(temp_level-test_D_mat[,qqqq],test_decision[qqqq,]) 
                }
              }
              CV_mat[fold,ppp]<-temp_res
            }
            
            
            print(c(fold,ppp))
          }
        }
        CV_seq <- apply(CV_mat, MARGIN = 2, sum)
        
        m_n <- m_n_candidate[(which.min(CV_seq))[1]]
        
        CV_time_our <- (proc.time() - start_time)["elapsed"]
      }
      start_time <- proc.time()
      train_D_mat <- D_mat[1:cut_off, ]
      
      
      cut_quantile <- seq(from = 0,to = 1,length.out = m_n+2)[-c(1, m_n+2)]
      
      knots_6 <- quantile(x = train_D_mat[, 5], probs = cut_quantile)
      
      price_6 <- function(alpha){
        basis_mat_6 <- bSpline(x = train_D_mat[,5],
                               knots = knots_6, degree = 3, 
                               Boundary.knots = c(lower_knots[6], upper_knots[6]),
                               intercept = TRUE)
        temp_level <- c(basis_mat_6 %*% alpha)
        return(sum(h_cost[6] * pmax(0, temp_level - train_D_mat[, 6]) + b_cost[6] * pmax(0, train_D_mat[, 6] - temp_level)))
      }
      
      price_6_empirical <- function(temp_level){
        return(sum(h_cost[6] * pmax(0, temp_level - train_D_mat[, 6]) + b_cost[6] * pmax(0, train_D_mat[, 6] - temp_level)))
      }
      
      
      res_6 <- solnp(pars = rep(1, m_n+4),
                     fun = price_6,
                     LB = rep(0, m_n+4),
                     UB = rep(Inf,m_n+4),
                     control = list(trace = 0))
      
      basis_mat_6 <- bSpline(x = train_D_mat[, 5],
                             knots = knots_6, degree = 3,
                             Boundary.knots = c(lower_knots[6], upper_knots[6]),
                             intercept = TRUE)
      best_decision_6 <- c(basis_mat_6 %*% res_6$par)
      
      knots_5 <- quantile(x = train_D_mat[, 4], probs = cut_quantile)
      
      price_5 <- function(alpha){
        basis_mat_5 <- bSpline(x = train_D_mat[,4],
                               knots = knots_5, degree = 3, 
                               Boundary.knots = c(lower_knots[5], upper_knots[5]),
                               intercept = TRUE)
        temp_level <- c(basis_mat_5 %*% alpha)
        return(sum(h_cost[5] * pmax(0, temp_level - train_D_mat[, 5]) + b_cost[5] * pmax(0, train_D_mat[, 5] - temp_level)) + price_6_empirical(pmax(temp_level - train_D_mat[,5], best_decision_6)))
      }
      
      price_5_empirical <- function(temp_level){
        return(sum(h_cost[5] * pmax(0, temp_level - train_D_mat[, 5]) + b_cost[5] * pmax(0, train_D_mat[, 5] - temp_level)) + price_6_empirical(pmax(temp_level - train_D_mat[,5], best_decision_6)))
      }
      
      
      res_5 <- solnp(pars = rep(1, m_n+4),
                     fun = price_5,
                     LB = rep(0, m_n+4),
                     UB = rep(Inf,m_n+4),
                     control = list(trace = 0))
      
      basis_mat_5 <- bSpline(x = train_D_mat[, 4],
                             knots = knots_5, degree = 3,
                             Boundary.knots = c(lower_knots[5], upper_knots[5]),
                             intercept = TRUE)
      best_decision_5 <- c(basis_mat_5 %*% res_5$par)
      
      
      knots_4 <- quantile(x = train_D_mat[, 3], probs = cut_quantile)
      
      
      price_4 <- function(alpha){
        basis_mat_4 <- bSpline(x = train_D_mat[,3],
                               knots = knots_4, degree = 3, 
                               Boundary.knots = c(lower_knots[4], upper_knots[4]),
                               intercept = TRUE)
        temp_level <- c(basis_mat_4 %*% alpha)
        return(sum(h_cost[4] * pmax(0, temp_level - train_D_mat[, 4]) + b_cost[4] * pmax(0, train_D_mat[, 4] - temp_level)) + price_5_empirical(pmax(temp_level - train_D_mat[,4], best_decision_5)))
      }
      
      price_4_empirical <- function(temp_level){
        return(sum(h_cost[4] * pmax(0, temp_level - train_D_mat[, 4]) + b_cost[4] * pmax(0, train_D_mat[, 4] - temp_level)) + price_5_empirical(pmax(temp_level - train_D_mat[,4], best_decision_5)))
      }
      
      
      res_4 <- solnp(pars = rep(1, m_n+4),
                     fun = price_4,
                     LB = rep(0, m_n+4),
                     UB = rep(Inf,m_n+4),
                     control = list(trace = 0))
      
      basis_mat_4 <- bSpline(x = train_D_mat[, 3],
                             knots = knots_4, degree = 3,
                             Boundary.knots = c(lower_knots[4], upper_knots[4]),
                             intercept = TRUE)
      best_decision_4 <- c(basis_mat_4 %*% res_4$par)
      
      knots_3 <- quantile(x = train_D_mat[, 2], probs = cut_quantile)
      
      
      price_3 <- function(alpha){
        basis_mat_3 <- bSpline(x = train_D_mat[,2],
                               knots = knots_3, degree = 3, 
                               Boundary.knots = c(lower_knots[3], upper_knots[3]),
                               intercept = TRUE)
        temp_level <- c(basis_mat_3 %*% alpha)
        return(sum(h_cost[3] * pmax(0, temp_level - train_D_mat[, 3]) + b_cost[3] * pmax(0, train_D_mat[, 3] - temp_level)) + price_4_empirical(pmax(temp_level - train_D_mat[,3], best_decision_4)))
      }
      
      price_3_empirical <- function(temp_level){
        return(sum(h_cost[3] * pmax(0, temp_level - train_D_mat[, 3]) + b_cost[3] * pmax(0, train_D_mat[, 3] - temp_level)) + price_4_empirical(pmax(temp_level - train_D_mat[,3], best_decision_4)))
      }
      
      
      res_3 <- solnp(pars = rep(1, m_n+4),
                     fun = price_3,
                     LB = rep(0, m_n+4),
                     UB = rep(Inf,m_n+4),
                     control = list(trace = 0))
      
      basis_mat_3 <- bSpline(x = train_D_mat[, 2],
                             knots = knots_3, degree = 3,
                             Boundary.knots = c(lower_knots[3], upper_knots[3]),
                             intercept = TRUE)
      best_decision_3 <- c(basis_mat_3 %*% res_3$par)
      
      
      knots_2 <- quantile(x = train_D_mat[, 1], probs = cut_quantile)
      
      
      price_2 <- function(alpha){
        basis_mat_2 <- bSpline(x = train_D_mat[,1],
                               knots = knots_2, degree = 3, 
                               Boundary.knots = c(lower_knots[2], upper_knots[2]),
                               intercept = TRUE)
        temp_level <- c(basis_mat_2 %*% alpha)
        return(sum(h_cost[2] * pmax(0, temp_level - train_D_mat[, 2]) + b_cost[2] * pmax(0, train_D_mat[, 2] - temp_level)) + price_3_empirical(pmax(temp_level - train_D_mat[,2], best_decision_3)))
      }
      
      price_2_empirical <- function(temp_level){
        return(sum(h_cost[2] * pmax(0, temp_level - train_D_mat[, 2]) + b_cost[2] * pmax(0, train_D_mat[, 2] - temp_level)) + price_3_empirical(pmax(temp_level - train_D_mat[,2], best_decision_3)))
      }
      
      
      res_2 <- solnp(pars = rep(1, m_n+4),
                     fun = price_2,
                     LB = rep(0, m_n+4),
                     UB = rep(Inf,m_n+4),
                     control = list(trace = 0))
      
      basis_mat_2 <- bSpline(x = train_D_mat[, 1],
                             knots = knots_2, degree = 3,
                             Boundary.knots = c(lower_knots[2], upper_knots[2]),
                             intercept = TRUE)
      best_decision_2 <- c(basis_mat_2 %*% res_2$par)
      
      price_1 <- function(alpha){
        return(sum(h_cost[1] * pmax(0, alpha - train_D_mat[, 1]) + b_cost[1] * pmax(0, train_D_mat[, 1] - alpha)) + price_2_empirical(pmax(alpha - train_D_mat[,1], best_decision_2)))
      }
      
      
      res_1 <- optimize(f = price_1, interval = c(0, 1000))
      
      best_decision_1 <- res_1$minimum
      
      test_D_mat <- full_D_mat[cut_off + 1, ]
      
      test_decision_1 <- best_decision_1
      
      temp_x <- test_D_mat[1]
      temp_base <- bSpline(x = temp_x, knots = knots_2, degree = 3,
                           Boundary.knots = c(lower_knots[2], upper_knots[2]),
                           intercept = TRUE)
      test_decision_2 <- c(temp_base%*%res_2$par)
      
      temp_x <- test_D_mat[2]
      temp_base <- bSpline(x = temp_x, knots = knots_3, degree = 3,
                           Boundary.knots = c(lower_knots[3], upper_knots[3]),
                           intercept = TRUE)
      test_decision_3 <- c(temp_base%*%res_3$par)
      
      
      temp_x <- test_D_mat[3]
      temp_base <- bSpline(x = temp_x, knots = knots_4, degree = 3,
                           Boundary.knots = c(lower_knots[4], upper_knots[4]),
                           intercept = TRUE)
      test_decision_4 <- c(temp_base%*%res_4$par)
      
      
      temp_x <- test_D_mat[4]
      temp_base <- bSpline(x = temp_x, knots = knots_5, degree = 3,
                           Boundary.knots = c(lower_knots[5], upper_knots[5]),
                           intercept = TRUE)
      test_decision_5 <- c(temp_base%*%res_5$par)
      
      
      temp_x <- test_D_mat[5]
      temp_base <- bSpline(x = temp_x, knots = knots_6, degree = 3,
                           Boundary.knots = c(lower_knots[6], upper_knots[6]),
                           intercept = TRUE)
      test_decision_6 <- c(temp_base%*%res_6$par)
      
      test_decision_our <- c(test_decision_1,
                             test_decision_2,
                             test_decision_3,
                             test_decision_4,
                             test_decision_5,
                             test_decision_6)
      run_time_our <- (proc.time() - start_time)["elapsed"]
    } 
    
    
    
    test_cost_our <- 0
    temp_level_our <- test_decision_our[1]
    for(qqqq in 1:T){
      test_cost_our<-test_cost_our+(h_cost[qqqq]*pmax(0,temp_level_our-test_D_mat[qqqq])+b_cost[qqqq]*pmax(0,test_D_mat[qqqq]-temp_level_our))
      
      if(qqqq<T){
        temp_level_our<-pmax(temp_level_our-test_D_mat[qqqq],test_decision_our[qqqq+1])
      }
    }
    cost_list_our <- test_cost_our
    
    
    
    return(c(cut_off, test_cost_our, CV_time_our, run_time_our, m_n))
  }
  res <- withTimeout({try(my_simulation_procedure(cut_off), silent = TRUE)}, 
                     timeout = 3600, onTimeout = "error")
  return(res)
}


clnum <- 28
cl <- makeCluster(getOption("cl.cores", clnum))
res <- parLapply(cl, 0:83,  simulation)
stopCluster(cl)

save(res, file = "spline_nocov.rdata")