rm(list = ls())
library(parallel)
simulation <- function(cut_off){
  library(R.utils)
  my_simulation_procedure<-function(cut_off){
    library(splines2)
    library(pracma)
    library(Rsolnp)
    chosen_index <- 1
    yw_time <- 180
    yw_time1 <- 600
    
    para_range <- 3
     
    multiplier1 <- 1
    multiplier2 <- 1 
    
    m_n_candidate <- 2:15
    
    
    
    
    
    load("bakery_data.rdata") 
    
    T <- 6
    h_cost <- rep(1, T)
    b_cost <- rep(1, T)
    n <- nrow(my_data)/T
    
    
    
    full_D_mat <- matrix(my_data[[products_all[indexing_order[chosen_index]]]],
                         n, T, byrow = TRUE)
    
    other_D_mat <- matrix(apply(my_data[products_all[-indexing_order[chosen_index]]], MARGIN = 1, sum),
                          n, T, byrow = TRUE)
    # indexing_order <- indexing_order[1:11]
    # indexing_order <- indexing_order[-chosen_index]
    # other_D_mat <- matrix(apply(my_data[products_all[indexing_order]], MARGIN = 1, sum),
    #                       n, T, byrow = TRUE)
    
    cov_dim <- c(0, rep(2, 5))
    
    full_X_array <- array(NA, dim = c(n, max(cov_dim), T))
    
    for(t in 2:T){
      full_X_array[, 1:2, t] <- cbind(full_D_mat[, t-1], other_D_mat[, t-1]) 
    }
    
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
      X_array <- full_X_array[1:cut_off, , ]
      
      lower_knots <- rep(NA, T)
      upper_knots <- rep(NA, T)
      
      for(t in 2:T){
        upper_knots[t] <- max(X_array[, 1, t]) + para_range * max(abs(X_array[, 2, t])) 
        lower_knots[t] <- min(X_array[, 1, t]) - para_range * max(abs(X_array[, 2, t])) 
        upper_knots[t] <- upper_knots[t] * multiplier1
        if(lower_knots[t] > 0){
          lower_knots[t] <- lower_knots[t] * multiplier2
        }else{
          lower_knots[t] <- lower_knots[t] * multiplier1
        }
      }
      
      if(cut_off < 5){
        m_n <- 2
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
          test_X_array <- X_array[test_index, , ]
          
          train_index <- c(my_split[-fold,])
          train_index <- train_index[which(train_index < 5e3)]
          train_D_mat <- D_mat[train_index, ]
          train_X_array <- X_array[train_index, , ]
          
          for(ppp in 1:length(m_n_candidate)){
            m_n <- m_n_candidate[ppp]
            
            knots_6 <- seq(from = lower_knots[6],
                           to = upper_knots[6],
                           length.out = m_n+2)[-c(1, m_n+2)]
            
            price_6 <- function(alpha){
              beta_yw <- c(1, alpha[1:(cov_dim[6] - 1)])
              
              basis_mat <- bSpline(x = c(train_X_array[, 1:cov_dim[6], 6] %*% beta_yw),
                                   knots = knots_6, degree = 3, 
                                   Boundary.knots = c(lower_knots[6], upper_knots[6]),
                                   intercept = TRUE)
              
              temp_level <- c(basis_mat %*% alpha[-c(1:(cov_dim[6] - 1))])
              return(sum(h_cost[6] * pmax(0, temp_level - train_D_mat[, 6]) + b_cost[6] * pmax(0, train_D_mat[, 6] - temp_level)))
            }
            
            price_6_empirical <- function(temp_level){
              return(sum(h_cost[6] * pmax(0, temp_level - train_D_mat[, 6]) + b_cost[6] * pmax(0, train_D_mat[, 6] - temp_level)))
            }
            
            
            res_6 <- withTimeout({solnp(pars = c(rep(0, cov_dim[6]-1), rep(1, m_n+4)),
                                        fun = price_6,
                                        LB = c(rep(-para_range, cov_dim[6]-1), rep(0, m_n+4)),
                                        UB = c(rep(para_range, cov_dim[6]-1), rep(Inf,m_n+4)),
                                        control = list(trace = 0))}, 
                                 timeout = yw_time, onTimeout = "warning")
            
            if(is.null(res_6)){
              res_6 <- optim(par = c(rep(0, cov_dim[6]-1), rep(1, m_n+4)), 
                             fn = price_6,
                             method = "L-BFGS-B",
                             lower = c(rep(-para_range, cov_dim[6]-1), rep(0, m_n+4)),
                             upper = c(rep(para_range, cov_dim[6]-1), rep(Inf,m_n+4)))
              res_6$pars <- res_6$par
            }
            
            beta_yw <- c(1, res_6$pars[1:(cov_dim[6] - 1)])
            basis_mat <- bSpline(x = c(train_X_array[, 1:cov_dim[6], 6] %*% beta_yw),
                                 knots = knots_6, degree = 3,
                                 Boundary.knots = c(lower_knots[6], upper_knots[6]),
                                 intercept = TRUE)
            best_decision_6 <- c(basis_mat %*% res_6$par[-c(1:(cov_dim[6] - 1))])
            
            
            knots_5 <- seq(from = lower_knots[5],
                           to = upper_knots[5],
                           length.out = m_n+2)[-c(1, m_n+2)]
            
            price_5 <- function(alpha){
              beta_yw <- c(1, alpha[1:(cov_dim[5] - 1)])
              
              basis_mat <- bSpline(x = c(train_X_array[, 1:cov_dim[5], 5] %*% beta_yw),
                                   knots = knots_5, degree = 3, 
                                   Boundary.knots = c(lower_knots[5], upper_knots[5]),
                                   intercept = TRUE)
              
              temp_level <- c(basis_mat %*% alpha[-c(1:(cov_dim[5] - 1))])
              return(sum(h_cost[5] * pmax(0, temp_level - train_D_mat[, 5]) + b_cost[5] * pmax(0, train_D_mat[, 5] - temp_level)) + price_6_empirical(pmax(temp_level - train_D_mat[, 5], best_decision_6)))
            }
            
            price_5_empirical <- function(temp_level){
              return(sum(h_cost[5] * pmax(0, temp_level - train_D_mat[, 5]) + b_cost[5] * pmax(0, train_D_mat[, 5] - temp_level)) + price_6_empirical(pmax(temp_level - train_D_mat[, 5], best_decision_6)))
            }
            
            
            res_5 <- withTimeout({solnp(pars = c(rep(0, cov_dim[5]-1), rep(1, m_n+4)),
                                        fun = price_5,
                                        LB = c(rep(-para_range, cov_dim[5]-1), rep(0, m_n+4)),
                                        UB = c(rep(para_range, cov_dim[5]-1), rep(Inf,m_n+4)),
                                        control = list(trace = 0))}, 
                                 timeout = yw_time, onTimeout = "warning")
            
            if(is.null(res_5)){
              res_5 <- optim(par = c(rep(0, cov_dim[5]-1), rep(1, m_n+4)), 
                             fn = price_5,
                             method = "L-BFGS-B",
                             lower = c(rep(-para_range, cov_dim[5]-1), rep(0, m_n+4)),
                             upper = c(rep(para_range, cov_dim[5]-1), rep(Inf,m_n+4)))
              res_5$pars <- res_5$par
            }
            
            beta_yw <- c(1, res_5$pars[1:(cov_dim[5] - 1)])
            basis_mat <- bSpline(x = c(train_X_array[, 1:cov_dim[5], 5] %*% beta_yw),
                                 knots = knots_5, degree = 3,
                                 Boundary.knots = c(lower_knots[5], upper_knots[5]),
                                 intercept = TRUE)
            best_decision_5 <- c(basis_mat %*% res_5$par[-c(1:(cov_dim[5] - 1))])
            
            
            
            knots_4 <- seq(from = lower_knots[4],
                           to = upper_knots[4],
                           length.out = m_n+2)[-c(1, m_n+2)]
            
            price_4 <- function(alpha){
              beta_yw <- c(1, alpha[1:(cov_dim[4] - 1)])
              
              basis_mat <- bSpline(x = c(train_X_array[, 1:cov_dim[4], 4] %*% beta_yw),
                                   knots = knots_4, degree = 3, 
                                   Boundary.knots = c(lower_knots[4], upper_knots[4]),
                                   intercept = TRUE)
              
              temp_level <- c(basis_mat %*% alpha[-c(1:(cov_dim[4] - 1))])
              return(sum(h_cost[4] * pmax(0, temp_level - train_D_mat[, 4]) + b_cost[4] * pmax(0, train_D_mat[, 4] - temp_level)) + price_5_empirical(pmax(temp_level - train_D_mat[, 4], best_decision_5)))
            }
            
            price_4_empirical <- function(temp_level){
              return(sum(h_cost[4] * pmax(0, temp_level - train_D_mat[, 4]) + b_cost[4] * pmax(0, train_D_mat[, 4] - temp_level)) + price_5_empirical(pmax(temp_level - train_D_mat[, 4], best_decision_5)))
            }
            
            
            res_4 <- withTimeout({solnp(pars = c(rep(0, cov_dim[4]-1), rep(1, m_n+4)),
                                        fun = price_4,
                                        LB = c(rep(-para_range, cov_dim[4]-1), rep(0, m_n+4)),
                                        UB = c(rep(para_range, cov_dim[4]-1), rep(Inf,m_n+4)),
                                        control = list(trace = 0))}, 
                                 timeout = yw_time, onTimeout = "warning")
            
            if(is.null(res_4)){
              res_4 <- optim(par = c(rep(0, cov_dim[4]-1), rep(1, m_n+4)), 
                             fn = price_4,
                             method = "L-BFGS-B",
                             lower = c(rep(-para_range, cov_dim[4]-1), rep(0, m_n+4)),
                             upper = c(rep(para_range, cov_dim[4]-1), rep(Inf,m_n+4)))
              res_4$pars <- res_4$par
            }
            
            beta_yw <- c(1, res_4$pars[1:(cov_dim[4] - 1)])
            basis_mat <- bSpline(x = c(train_X_array[, 1:cov_dim[4], 4] %*% beta_yw),
                                 knots = knots_4, degree = 3,
                                 Boundary.knots = c(lower_knots[4], upper_knots[4]),
                                 intercept = TRUE)
            best_decision_4 <- c(basis_mat %*% res_4$par[-c(1:(cov_dim[4] - 1))])
            
            
            knots_3 <- seq(from = lower_knots[3],
                           to = upper_knots[3],
                           length.out = m_n+2)[-c(1, m_n+2)]
            
            price_3 <- function(alpha){
              beta_yw <- c(1, alpha[1:(cov_dim[3] - 1)])
              
              basis_mat <- bSpline(x = c(train_X_array[, 1:cov_dim[3], 3] %*% beta_yw),
                                   knots = knots_3, degree = 3, 
                                   Boundary.knots = c(lower_knots[3], upper_knots[3]),
                                   intercept = TRUE)
              
              temp_level <- c(basis_mat %*% alpha[-c(1:(cov_dim[3] - 1))])
              return(sum(h_cost[3] * pmax(0, temp_level - train_D_mat[, 3]) + b_cost[3] * pmax(0, train_D_mat[, 3] - temp_level)) + price_4_empirical(pmax(temp_level - train_D_mat[, 3], best_decision_4)))
            }
            
            price_3_empirical <- function(temp_level){
              return(sum(h_cost[3] * pmax(0, temp_level - train_D_mat[, 3]) + b_cost[3] * pmax(0, train_D_mat[, 3] - temp_level)) + price_4_empirical(pmax(temp_level - train_D_mat[, 3], best_decision_4)))
            }
            
            
            res_3 <- withTimeout({solnp(pars = c(rep(0, cov_dim[3]-1), rep(1, m_n+4)),
                                        fun = price_3,
                                        LB = c(rep(-para_range, cov_dim[3]-1), rep(0, m_n+4)),
                                        UB = c(rep(para_range, cov_dim[3]-1), rep(Inf,m_n+4)),
                                        control = list(trace = 0))}, 
                                 timeout = yw_time, onTimeout = "warning")
            
            if(is.null(res_3)){
              res_3 <- optim(par = c(rep(0, cov_dim[3]-1), rep(1, m_n+4)), 
                             fn = price_3,
                             method = "L-BFGS-B",
                             lower = c(rep(-para_range, cov_dim[3]-1), rep(0, m_n+4)),
                             upper = c(rep(para_range, cov_dim[3]-1), rep(Inf,m_n+4)))
              res_3$pars <- res_3$par
            }
            
            
            beta_yw <- c(1, res_3$pars[1:(cov_dim[3] - 1)])
            basis_mat <- bSpline(x = c(train_X_array[, 1:cov_dim[3], 3] %*% beta_yw),
                                 knots = knots_3, degree = 3,
                                 Boundary.knots = c(lower_knots[3], upper_knots[3]),
                                 intercept = TRUE)
            best_decision_3 <- c(basis_mat %*% res_3$par[-c(1:(cov_dim[3] - 1))])
            
            knots_2 <- seq(from = lower_knots[2],
                           to = upper_knots[2],
                           length.out = m_n+2)[-c(1, m_n+2)]
            
            price_2 <- function(alpha){
              beta_yw <- c(1, alpha[1:(cov_dim[2] - 1)])
              
              basis_mat <- bSpline(x = c(train_X_array[, 1:cov_dim[2], 2] %*% beta_yw),
                                   knots = knots_2, degree = 3, 
                                   Boundary.knots = c(lower_knots[2], upper_knots[2]),
                                   intercept = TRUE)
              
              temp_level <- c(basis_mat %*% alpha[-c(1:(cov_dim[2] - 1))])
              return(sum(h_cost[2] * pmax(0, temp_level - train_D_mat[, 2]) + b_cost[2] * pmax(0, train_D_mat[, 2] - temp_level)) + price_3_empirical(pmax(temp_level - train_D_mat[, 2], best_decision_3)))
            }
            
            price_2_empirical <- function(temp_level){
              return(sum(h_cost[2] * pmax(0, temp_level - train_D_mat[, 2]) + b_cost[2] * pmax(0, train_D_mat[, 2] - temp_level)) + price_3_empirical(pmax(temp_level - train_D_mat[, 2], best_decision_3)))
            }
            
            
            res_2 <- withTimeout({solnp(pars = c(rep(0, cov_dim[2]-1), rep(1, m_n+4)),
                                        fun = price_2,
                                        LB = c(rep(-para_range, cov_dim[2]-1), rep(0, m_n+4)),
                                        UB = c(rep(para_range, cov_dim[2]-1), rep(Inf,m_n+4)),
                                        control = list(trace = 0))}, 
                                 timeout = yw_time, onTimeout = "warning")
            
            if(is.null(res_2)){
              res_2 <- optim(par = c(rep(0, cov_dim[2]-1), rep(1, m_n+4)), 
                             fn = price_2,
                             method = "L-BFGS-B",
                             lower = c(rep(-para_range, cov_dim[2]-1), rep(0, m_n+4)),
                             upper = c(rep(para_range, cov_dim[2]-1), rep(Inf,m_n+4)))
              res_2$pars <- res_2$par
            }
            
            beta_yw <- c(1, res_2$pars[1:(cov_dim[2] - 1)])
            basis_mat <- bSpline(x = c(train_X_array[, 1:cov_dim[2], 2] %*% beta_yw),
                                 knots = knots_2, degree = 3,
                                 Boundary.knots = c(lower_knots[2], upper_knots[2]),
                                 intercept = TRUE)
            best_decision_2 <- c(basis_mat %*% res_2$par[-c(1:(cov_dim[2] - 1))])
            
            
            price_1 <- function(alpha){
              return(sum(h_cost[1] * pmax(0, alpha - train_D_mat[, 1]) + b_cost[1] * pmax(0, train_D_mat[, 1] - alpha)) + price_2_empirical(pmax(alpha - train_D_mat[,1], best_decision_2)))
            }
            
            
            res_1 <- optimize(f = price_1, interval = c(0, 1000))
            
            best_decision_1 <- res_1$minimum
            
            test_decision_1 <- best_decision_1
            
            if(is.null(nrow(test_D_mat))){
              temp_x <- sum(test_X_array[1:cov_dim[2], 2] * c(1, res_2$pars[1:(cov_dim[2] - 1)])) 
              temp_base <- bSpline(x = temp_x, knots = knots_2, degree = 3,
                                   Boundary.knots = c(lower_knots[2], upper_knots[2]),
                                   intercept = TRUE)
              test_decision_2 <- c(temp_base%*%res_2$par[-c(1:(cov_dim[2] - 1))])
              
              temp_x <- sum(test_X_array[1:cov_dim[3], 3] * c(1, res_3$pars[1:(cov_dim[3] - 1)])) 
              temp_base <- bSpline(x = temp_x, knots = knots_3, degree = 3,
                                   Boundary.knots = c(lower_knots[3], upper_knots[3]),
                                   intercept = TRUE)
              test_decision_3 <- c(temp_base%*%res_3$par[-c(1:(cov_dim[3] - 1))])
              
              temp_x <- sum(test_X_array[1:cov_dim[4], 4] * c(1, res_4$pars[1:(cov_dim[4] - 1)]))           
              temp_base <- bSpline(x = temp_x, knots = knots_4, degree = 3,
                                   Boundary.knots = c(lower_knots[4], upper_knots[4]),
                                   intercept = TRUE)
              test_decision_4 <- c(temp_base%*%res_4$par[-c(1:(cov_dim[4] - 1))])
              
              
              temp_x <- sum(test_X_array[1:cov_dim[5], 5] * c(1, res_5$pars[1:(cov_dim[5] - 1)]))
              temp_base <- bSpline(x = temp_x, knots = knots_5, degree = 3,
                                   Boundary.knots = c(lower_knots[5], upper_knots[5]),
                                   intercept = TRUE)
              test_decision_5 <- c(temp_base%*%res_5$par[-c(1:(cov_dim[5] - 1))])
              
              
              temp_x <- sum(test_X_array[1:cov_dim[6], 6] * c(1, res_6$pars[1:(cov_dim[6] - 1)]))
              temp_base <- bSpline(x = temp_x, knots = knots_6, degree = 3,
                                   Boundary.knots = c(lower_knots[6], upper_knots[6]),
                                   intercept = TRUE)
              test_decision_6 <- c(temp_base%*%res_6$par[-c(1:(cov_dim[6] - 1))])
              
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
              temp_x <- c(test_X_array[, 1:cov_dim[2], 2] %*% c(1, res_2$pars[1:(cov_dim[2] - 1)]))  
              temp_base <- bSpline(x = temp_x,
                                   knots = knots_2,
                                   degree = 3,
                                   Boundary.knots = c(lower_knots[2], upper_knots[2]),
                                   intercept = TRUE)
              test_decision_2 <- c(temp_base%*%res_2$par[-c(1:(cov_dim[2] - 1))])
              
              temp_x <- c(test_X_array[, 1:cov_dim[3], 3] %*% c(1, res_3$pars[1:(cov_dim[3] - 1)]))
              temp_base <- bSpline(x = temp_x,
                                   knots = knots_3,
                                   degree = 3,
                                   Boundary.knots = c(lower_knots[3], upper_knots[3]),
                                   intercept = TRUE)
              test_decision_3 <- c(temp_base%*%res_3$par[-c(1:(cov_dim[3] - 1))])
              
              
              temp_x <- c(test_X_array[, 1:cov_dim[4], 4] %*% c(1, res_4$pars[1:(cov_dim[4] - 1)]))
              temp_base <- bSpline(x = temp_x,
                                   knots = knots_4,
                                   degree = 3,
                                   Boundary.knots = c(lower_knots[4], upper_knots[4]),
                                   intercept = TRUE)
              test_decision_4 <- c(temp_base%*%res_4$par[-c(1:(cov_dim[4] - 1))])
              
              temp_x <- c(test_X_array[, 1:cov_dim[5], 5] %*% c(1, res_5$pars[1:(cov_dim[5] - 1)]))
              temp_base <- bSpline(x = temp_x,
                                   knots = knots_5,
                                   degree = 3,
                                   Boundary.knots = c(lower_knots[5], upper_knots[5]),
                                   intercept = TRUE)
              test_decision_5 <- c(temp_base%*%res_5$par[-c(1:(cov_dim[5] - 1))])
              
              temp_x <- c(test_X_array[, 1:cov_dim[6], 6] %*% c(1, res_6$pars[1:(cov_dim[6] - 1)]))
              temp_base <- bSpline(x = temp_x,
                                   knots = knots_6,
                                   degree = 3,
                                   Boundary.knots = c(lower_knots[6], upper_knots[6]),
                                   intercept = TRUE)
              test_decision_6 <- c(temp_base%*%res_6$par[-c(1:(cov_dim[6] - 1))])
              
              
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
      train_X_array <- X_array[1:cut_off, , ]
      
      knots_6 <- seq(from = lower_knots[6],
                     to = upper_knots[6],
                     length.out = m_n+2)[-c(1, m_n+2)]
      
      price_6 <- function(alpha){
        beta_yw <- c(1, alpha[1:(cov_dim[6] - 1)])
        
        basis_mat <- bSpline(x = c(train_X_array[, 1:cov_dim[6], 6] %*% beta_yw),
                             knots = knots_6, degree = 3, 
                             Boundary.knots = c(lower_knots[6], upper_knots[6]),
                             intercept = TRUE)
        
        temp_level <- c(basis_mat %*% alpha[-c(1:(cov_dim[6] - 1))])
        return(sum(h_cost[6] * pmax(0, temp_level - train_D_mat[, 6]) + b_cost[6] * pmax(0, train_D_mat[, 6] - temp_level)))
      }
      
      price_6_empirical <- function(temp_level){
        return(sum(h_cost[6] * pmax(0, temp_level - train_D_mat[, 6]) + b_cost[6] * pmax(0, train_D_mat[, 6] - temp_level)))
      }
      
      
      res_6 <- withTimeout({solnp(pars = c(rep(0, cov_dim[6]-1), rep(1, m_n+4)),
                                  fun = price_6,
                                  LB = c(rep(-para_range, cov_dim[6]-1), rep(0, m_n+4)),
                                  UB = c(rep(para_range, cov_dim[6]-1), rep(Inf,m_n+4)),
                                  control = list(trace = 0))}, 
                           timeout = yw_time1, onTimeout = "warning")
      
      if(is.null(res_6)){
        res_6 <- optim(par = c(rep(0, cov_dim[6]-1), rep(1, m_n+4)), 
                       fn = price_6,
                       method = "L-BFGS-B",
                       lower = c(rep(-para_range, cov_dim[6]-1), rep(0, m_n+4)),
                       upper = c(rep(para_range, cov_dim[6]-1), rep(Inf,m_n+4)))
        res_6$pars <- res_6$par
      }
      
      beta_yw <- c(1, res_6$pars[1:(cov_dim[6] - 1)])
      basis_mat <- bSpline(x = c(train_X_array[, 1:cov_dim[6], 6] %*% beta_yw),
                           knots = knots_6, degree = 3,
                           Boundary.knots = c(lower_knots[6], upper_knots[6]),
                           intercept = TRUE)
      best_decision_6 <- c(basis_mat %*% res_6$par[-c(1:(cov_dim[6] - 1))])
      
      
      knots_5 <- seq(from = lower_knots[5],
                     to = upper_knots[5],
                     length.out = m_n+2)[-c(1, m_n+2)]
      
      price_5 <- function(alpha){
        beta_yw <- c(1, alpha[1:(cov_dim[5] - 1)])
        
        basis_mat <- bSpline(x = c(train_X_array[, 1:cov_dim[5], 5] %*% beta_yw),
                             knots = knots_5, degree = 3, 
                             Boundary.knots = c(lower_knots[5], upper_knots[5]),
                             intercept = TRUE)
        
        temp_level <- c(basis_mat %*% alpha[-c(1:(cov_dim[5] - 1))])
        return(sum(h_cost[5] * pmax(0, temp_level - train_D_mat[, 5]) + b_cost[5] * pmax(0, train_D_mat[, 5] - temp_level)) + price_6_empirical(pmax(temp_level - train_D_mat[, 5], best_decision_6)))
      }
      
      price_5_empirical <- function(temp_level){
        return(sum(h_cost[5] * pmax(0, temp_level - train_D_mat[, 5]) + b_cost[5] * pmax(0, train_D_mat[, 5] - temp_level)) + price_6_empirical(pmax(temp_level - train_D_mat[, 5], best_decision_6)))
      }
      
      
      res_5 <- withTimeout({solnp(pars = c(rep(0, cov_dim[5]-1), rep(1, m_n+4)),
                                  fun = price_5,
                                  LB = c(rep(-para_range, cov_dim[5]-1), rep(0, m_n+4)),
                                  UB = c(rep(para_range, cov_dim[5]-1), rep(Inf,m_n+4)),
                                  control = list(trace = 0))}, 
                           timeout = yw_time1, onTimeout = "warning")
      
      if(is.null(res_5)){
        res_5 <- optim(par = c(rep(0, cov_dim[5]-1), rep(1, m_n+4)), 
                       fn = price_5,
                       method = "L-BFGS-B",
                       lower = c(rep(-para_range, cov_dim[5]-1), rep(0, m_n+4)),
                       upper = c(rep(para_range, cov_dim[5]-1), rep(Inf,m_n+4)))
        res_5$pars <- res_5$par
      }
      
      beta_yw <- c(1, res_5$pars[1:(cov_dim[5] - 1)])
      basis_mat <- bSpline(x = c(train_X_array[, 1:cov_dim[5], 5] %*% beta_yw),
                           knots = knots_5, degree = 3,
                           Boundary.knots = c(lower_knots[5], upper_knots[5]),
                           intercept = TRUE)
      best_decision_5 <- c(basis_mat %*% res_5$par[-c(1:(cov_dim[5] - 1))])
      
      
      
      knots_4 <- seq(from = lower_knots[4],
                     to = upper_knots[4],
                     length.out = m_n+2)[-c(1, m_n+2)]
      
      price_4 <- function(alpha){
        beta_yw <- c(1, alpha[1:(cov_dim[4] - 1)])
        
        basis_mat <- bSpline(x = c(train_X_array[, 1:cov_dim[4], 4] %*% beta_yw),
                             knots = knots_4, degree = 3, 
                             Boundary.knots = c(lower_knots[4], upper_knots[4]),
                             intercept = TRUE)
        
        temp_level <- c(basis_mat %*% alpha[-c(1:(cov_dim[4] - 1))])
        return(sum(h_cost[4] * pmax(0, temp_level - train_D_mat[, 4]) + b_cost[4] * pmax(0, train_D_mat[, 4] - temp_level)) + price_5_empirical(pmax(temp_level - train_D_mat[, 4], best_decision_5)))
      }
      
      price_4_empirical <- function(temp_level){
        return(sum(h_cost[4] * pmax(0, temp_level - train_D_mat[, 4]) + b_cost[4] * pmax(0, train_D_mat[, 4] - temp_level)) + price_5_empirical(pmax(temp_level - train_D_mat[, 4], best_decision_5)))
      }
      
      
      res_4 <- withTimeout({solnp(pars = c(rep(0, cov_dim[4]-1), rep(1, m_n+4)),
                                  fun = price_4,
                                  LB = c(rep(-para_range, cov_dim[4]-1), rep(0, m_n+4)),
                                  UB = c(rep(para_range, cov_dim[4]-1), rep(Inf,m_n+4)),
                                  control = list(trace = 0))}, 
                           timeout = yw_time1, onTimeout = "warning")
      
      if(is.null(res_4)){
        res_4 <- optim(par = c(rep(0, cov_dim[4]-1), rep(1, m_n+4)), 
                       fn = price_4,
                       method = "L-BFGS-B",
                       lower = c(rep(-para_range, cov_dim[4]-1), rep(0, m_n+4)),
                       upper = c(rep(para_range, cov_dim[4]-1), rep(Inf,m_n+4)))
        res_4$pars <- res_4$par
      }
      
      beta_yw <- c(1, res_4$pars[1:(cov_dim[4] - 1)])
      basis_mat <- bSpline(x = c(train_X_array[, 1:cov_dim[4], 4] %*% beta_yw),
                           knots = knots_4, degree = 3,
                           Boundary.knots = c(lower_knots[4], upper_knots[4]),
                           intercept = TRUE)
      best_decision_4 <- c(basis_mat %*% res_4$par[-c(1:(cov_dim[4] - 1))])
      
      
      knots_3 <- seq(from = lower_knots[3],
                     to = upper_knots[3],
                     length.out = m_n+2)[-c(1, m_n+2)]
      
      price_3 <- function(alpha){
        beta_yw <- c(1, alpha[1:(cov_dim[3] - 1)])
        
        basis_mat <- bSpline(x = c(train_X_array[, 1:cov_dim[3], 3] %*% beta_yw),
                             knots = knots_3, degree = 3, 
                             Boundary.knots = c(lower_knots[3], upper_knots[3]),
                             intercept = TRUE)
        
        temp_level <- c(basis_mat %*% alpha[-c(1:(cov_dim[3] - 1))])
        return(sum(h_cost[3] * pmax(0, temp_level - train_D_mat[, 3]) + b_cost[3] * pmax(0, train_D_mat[, 3] - temp_level)) + price_4_empirical(pmax(temp_level - train_D_mat[, 3], best_decision_4)))
      }
      
      price_3_empirical <- function(temp_level){
        return(sum(h_cost[3] * pmax(0, temp_level - train_D_mat[, 3]) + b_cost[3] * pmax(0, train_D_mat[, 3] - temp_level)) + price_4_empirical(pmax(temp_level - train_D_mat[, 3], best_decision_4)))
      }
      
      
      res_3 <- withTimeout({solnp(pars = c(rep(0, cov_dim[3]-1), rep(1, m_n+4)),
                                  fun = price_3,
                                  LB = c(rep(-para_range, cov_dim[3]-1), rep(0, m_n+4)),
                                  UB = c(rep(para_range, cov_dim[3]-1), rep(Inf,m_n+4)),
                                  control = list(trace = 0))}, 
                           timeout = yw_time1, onTimeout = "warning")
      
      if(is.null(res_3)){
        res_3 <- optim(par = c(rep(0, cov_dim[3]-1), rep(1, m_n+4)), 
                       fn = price_3,
                       method = "L-BFGS-B",
                       lower = c(rep(-para_range, cov_dim[3]-1), rep(0, m_n+4)),
                       upper = c(rep(para_range, cov_dim[3]-1), rep(Inf,m_n+4)))
        res_3$pars <- res_3$par
      }
      
      
      beta_yw <- c(1, res_3$pars[1:(cov_dim[3] - 1)])
      basis_mat <- bSpline(x = c(train_X_array[, 1:cov_dim[3], 3] %*% beta_yw),
                           knots = knots_3, degree = 3,
                           Boundary.knots = c(lower_knots[3], upper_knots[3]),
                           intercept = TRUE)
      best_decision_3 <- c(basis_mat %*% res_3$par[-c(1:(cov_dim[3] - 1))])
      
      knots_2 <- seq(from = lower_knots[2],
                     to = upper_knots[2],
                     length.out = m_n+2)[-c(1, m_n+2)]
      
      price_2 <- function(alpha){
        beta_yw <- c(1, alpha[1:(cov_dim[2] - 1)])
        
        basis_mat <- bSpline(x = c(train_X_array[, 1:cov_dim[2], 2] %*% beta_yw),
                             knots = knots_2, degree = 3, 
                             Boundary.knots = c(lower_knots[2], upper_knots[2]),
                             intercept = TRUE)
        
        temp_level <- c(basis_mat %*% alpha[-c(1:(cov_dim[2] - 1))])
        return(sum(h_cost[2] * pmax(0, temp_level - train_D_mat[, 2]) + b_cost[2] * pmax(0, train_D_mat[, 2] - temp_level)) + price_3_empirical(pmax(temp_level - train_D_mat[, 2], best_decision_3)))
      }
      
      price_2_empirical <- function(temp_level){
        return(sum(h_cost[2] * pmax(0, temp_level - train_D_mat[, 2]) + b_cost[2] * pmax(0, train_D_mat[, 2] - temp_level)) + price_3_empirical(pmax(temp_level - train_D_mat[, 2], best_decision_3)))
      }
      
      
      res_2 <- withTimeout({solnp(pars = c(rep(0, cov_dim[2]-1), rep(1, m_n+4)),
                                  fun = price_2,
                                  LB = c(rep(-para_range, cov_dim[2]-1), rep(0, m_n+4)),
                                  UB = c(rep(para_range, cov_dim[2]-1), rep(Inf,m_n+4)),
                                  control = list(trace = 0))}, 
                           timeout = yw_time1, onTimeout = "warning")
      
      if(is.null(res_2)){
        res_2 <- optim(par = c(rep(0, cov_dim[2]-1), rep(1, m_n+4)), 
                       fn = price_2,
                       method = "L-BFGS-B",
                       lower = c(rep(-para_range, cov_dim[2]-1), rep(0, m_n+4)),
                       upper = c(rep(para_range, cov_dim[2]-1), rep(Inf,m_n+4)))
        res_2$pars <- res_2$par
      }
      
      beta_yw <- c(1, res_2$pars[1:(cov_dim[2] - 1)])
      basis_mat <- bSpline(x = c(train_X_array[, 1:cov_dim[2], 2] %*% beta_yw),
                           knots = knots_2, degree = 3,
                           Boundary.knots = c(lower_knots[2], upper_knots[2]),
                           intercept = TRUE)
      best_decision_2 <- c(basis_mat %*% res_2$par[-c(1:(cov_dim[2] - 1))])
      
      
      price_1 <- function(alpha){
        return(sum(h_cost[1] * pmax(0, alpha - train_D_mat[, 1]) + b_cost[1] * pmax(0, train_D_mat[, 1] - alpha)) + price_2_empirical(pmax(alpha - train_D_mat[,1], best_decision_2)))
      }
      
      
      res_1 <- optimize(f = price_1, interval = c(0, 1000))
      
      best_decision_1 <- res_1$minimum
      
      test_D_mat <- full_D_mat[cut_off + 1, ]
      test_X_array <- full_X_array[cut_off + 1, , ]
      
      test_decision_1 <- best_decision_1
      
      temp_x <- sum(test_X_array[1:cov_dim[2], 2] * c(1, res_2$pars[1:(cov_dim[2] - 1)])) 
      temp_base <- bSpline(x = temp_x, knots = knots_2, degree = 3,
                           Boundary.knots = c(lower_knots[2], upper_knots[2]),
                           intercept = TRUE)
      test_decision_2 <- c(temp_base%*%res_2$par[-c(1:(cov_dim[2] - 1))])
      
      temp_x <- sum(test_X_array[1:cov_dim[3], 3] * c(1, res_3$pars[1:(cov_dim[3] - 1)])) 
      temp_base <- bSpline(x = temp_x, knots = knots_3, degree = 3,
                           Boundary.knots = c(lower_knots[3], upper_knots[3]),
                           intercept = TRUE)
      test_decision_3 <- c(temp_base%*%res_3$par[-c(1:(cov_dim[3] - 1))])
      
      temp_x <- sum(test_X_array[1:cov_dim[4], 4] * c(1, res_4$pars[1:(cov_dim[4] - 1)]))           
      temp_base <- bSpline(x = temp_x, knots = knots_4, degree = 3,
                           Boundary.knots = c(lower_knots[4], upper_knots[4]),
                           intercept = TRUE)
      test_decision_4 <- c(temp_base%*%res_4$par[-c(1:(cov_dim[4] - 1))])
      
      
      temp_x <- sum(test_X_array[1:cov_dim[5], 5] * c(1, res_5$pars[1:(cov_dim[5] - 1)]))
      temp_base <- bSpline(x = temp_x, knots = knots_5, degree = 3,
                           Boundary.knots = c(lower_knots[5], upper_knots[5]),
                           intercept = TRUE)
      test_decision_5 <- c(temp_base%*%res_5$par[-c(1:(cov_dim[5] - 1))])
      
      
      temp_x <- sum(test_X_array[1:cov_dim[6], 6] * c(1, res_6$pars[1:(cov_dim[6] - 1)]))
      temp_base <- bSpline(x = temp_x, knots = knots_6, degree = 3,
                           Boundary.knots = c(lower_knots[6], upper_knots[6]),
                           intercept = TRUE)
      test_decision_6 <- c(temp_base%*%res_6$par[-c(1:(cov_dim[6] - 1))])
      
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
  res <- try(my_simulation_procedure(cut_off),silent = TRUE)
  return(res)
}


clnum <- 28
cl <- makeCluster(getOption("cl.cores", clnum))
res <- parLapply(cl, 0:83,  simulation)
stopCluster(cl)

save(res, file = "spline_cov.rdata")