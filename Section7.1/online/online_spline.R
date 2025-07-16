rm(list=ls())
library(parallel)
iter <- 1
simulation<-function(iter){
  my_simulation_procedure<-function(iter){
    library(splines2)
    library(pracma)
    library(Rsolnp)
    load("ye_cov_n200.rdata")
    full_D_mat <- res[[iter]][[26]]
    full_X_mat <- res[[iter]][[27]]
    
    load("cov_2_true1.rdata")
    
    opt_knots_1 <- res[[1]][[1]]
    opt_knots_2 <- res[[1]][[2]]
    
    
    opt_pars1 <- res[[2]][[1]]
    opt_pars2 <- res[[2]][[2]]
    
    b_cost<-c(5,1,1)
    h_cost<-c(1,5,1)
    
    tuning_length <- 10
    
    window_length <- 50
    my_online_res <- matrix(NA, window_length, 3)
    opt_costs <- rep(NA, window_length)
    
    for(cut_off in 0:(window_length - 1)){
      if(cut_off > 0){
        D_mat <- full_D_mat[1:cut_off, ]
        X_mat <- full_X_mat[1:cut_off, ]
      }
      
      X_test <- full_X_mat[cut_off + 1, ]
      D_test <- full_D_mat[cut_off + 1, ]
      
      ptm <- proc.time()
      # we update the tuning parameters every 10 selling seasons
      if((cut_off + 1) %% tuning_length == 1){
        m_n_candidate <- 1:10
        
        if(cut_off == 0){
          # at the first selling season, we randomly choose the tuning parameters
          m_n <- sample(x = m_n_candidate, size = 1)
          
        }else{
          # update the tuning parameters based on cross validation
          # proposed spline approach
          CV_mat<-matrix(0,5,length(m_n_candidate))
          my_split<-1:cut_off
          if(cut_off%%5 == 0){
            my_split<-matrix(my_split,5, cut_off/5,byrow = TRUE)
          }else{
            my_split <- c(my_split, rep(1e4, (5*ceiling(cut_off/5) - cut_off)))
            my_split <- matrix(my_split, 5, ceiling(cut_off/5))
          }
          for(fold in 1:5){
            test_index <- c(my_split[fold,])
            test_index <- test_index[which(test_index < 5e3)]
            
            train_index <- c(my_split[-fold,])
            train_index <- train_index[which(train_index < 5e3)]
            
            test_X<-D_mat[test_index,1]
            test_Y<-D_mat[test_index,2]
            test_Z1<-X_mat[test_index,1]
            test_Z2<-X_mat[test_index,2]
            
            train_X<-D_mat[train_index,1]
            train_Y<-D_mat[train_index,2]
            train_Z1<-X_mat[train_index,1]
            train_Z2<-X_mat[train_index,2]
            
            for(ppp in 1:length(m_n_candidate)){
              m_n<-m_n_candidate[ppp]
              knots_2<-seq(from=0.9,to=13.5,length.out=m_n+2)[-c(1,m_n+2)]
              basis_mat_2<-bSpline(x=train_X,knots = knots_2,degree = 3,Boundary.knots = c(0.9,13.5),intercept = TRUE )
              
              price_2<-function(alpha){
                basis_mat_2<-bSpline(x=train_X+alpha[1]*train_Z2,knots = knots_2,degree = 3,Boundary.knots = c(0.9,13.5),intercept = TRUE )
                temp_level<-c(basis_mat_2%*%alpha[-1])
                return(sum(h_cost[2]*pmax(0,temp_level-train_Y)+b_cost[2]*pmax(0,train_Y-temp_level)))
              }
              
              price_2_empirical<-function(temp_level_2){
                return(sum(h_cost[2]*pmax(0,temp_level_2-train_Y)+b_cost[2]*pmax(0,train_Y-temp_level_2)))
              }
              
              res_2 <- optim(par = c(0.1,rep(0.1,m_n+4)), fn = price_2,
                             method = "L-BFGS-B",
                             lower = c(0,rep(0,m_n+4)),
                             upper = c(2,rep(Inf,m_n+4)))
              
              # res_2<-solnp(pars = c(0.1,rep(0.1,m_n+4)),fun = price_2,LB=c(0,rep(0,m_n+4)),UB=c(2,rep(Inf,m_n+4)),control = list(trace=0))
              basis_mat_2<-bSpline(x=train_X+res_2$par[1]*train_Z2,knots = knots_2,degree = 3,Boundary.knots = c(0.9,13.5),intercept = TRUE )
              best_decision_2<-c(basis_mat_2%*%res_2$par[-1])
              
              knots_1<-seq(from=1,to=5,length.out=m_n+2)[-c(1,m_n+2)]
              basis_mat_1<-bSpline(x=train_Z1,knots = knots_1,degree = 3,Boundary.knots = c(1,5),intercept = TRUE )
              price_1<-function(alpha){
                temp_level<-c(basis_mat_1%*%alpha)
                return(sum(h_cost[1]*pmax(0,temp_level-train_X)+b_cost[1]*pmax(0,train_X-temp_level))+price_2_empirical(pmax(temp_level-train_X,best_decision_2)))
              }
              
              res_1 <- optim(par = rep(0.1,m_n+4), fn = price_1,
                             method = "L-BFGS-B",
                             lower = rep(0,m_n+4))
              # res_1<-solnp(pars = rep(0.1,m_n+4),fun = price_1,LB=rep(0,m_n+4),control = list(trace=0))
              
              basis_test_1<-bSpline(x=test_Z1,knots = knots_1,degree = 3,Boundary.knots = c(1,5),intercept = TRUE )
              test_decision_1<-c(basis_test_1%*%res_1$par)
              
              temp_x<-c(test_Z2*res_2$par[1]+test_X)
              temp_base<-bSpline(x=temp_x,knots = knots_2,degree = 3,Boundary.knots = c(0.9,13.5),intercept = TRUE )
              test_decision_2<-c(temp_base%*%res_2$par[-1])
              
              temp_res<-0
              for(ooo in 1:length(test_X)){
                temp_res<-temp_res+h_cost[1]*max(0,test_decision_1[ooo]-test_X[ooo])+b_cost[1]*max(0,test_X[ooo]-test_decision_1[ooo])
                temp_res<-temp_res+h_cost[2]*max(0,max(test_decision_1[ooo]-test_X[ooo],test_decision_2[ooo])-test_Y[ooo])+b_cost[2]*max(0,test_Y[ooo]-max(test_decision_1[ooo]-test_X[ooo],test_decision_2[ooo]))
              }
              CV_mat[fold,ppp]<-temp_res
              print(c(fold,ppp))
            }
          }
          CV_seq<-apply(CV_mat,MARGIN = 2,sum)
          
          m_n<-m_n_candidate[(which.min(CV_seq))[1]]
        }
        
      }
      CVtime <- as.numeric((proc.time()-ptm)[3])
      
      ptm<-proc.time()
      
      if(cut_off == 0){
        test_decision_1 <- 6
        test_decision_2 <- 1
      }
      
      if(cut_off == 1){
        test_decision_1 <- D_mat[1]
        test_decision_2 <- D_mat[2]
      }
      
      if(cut_off > 1){
        train_X<-D_mat[,1]
        train_Y<-D_mat[,2]
        train_Z1<-X_mat[,1]
        train_Z2<-X_mat[,2]
        
        knots_2<-seq(from=0.9,to=13.5,length.out=m_n+2)[-c(1,m_n+2)]
        basis_mat_2<-bSpline(x=train_X,knots = knots_2,degree = 3,Boundary.knots = c(0.9,13.5),intercept = TRUE )
        
        price_2<-function(alpha){
          basis_mat_2<-bSpline(x=train_X+alpha[1]*train_Z2,knots = knots_2,degree = 3,Boundary.knots = c(0.9,13.5),intercept = TRUE )
          temp_level<-c(basis_mat_2%*%alpha[-1])
          return(sum(h_cost[2]*pmax(0,temp_level-train_Y)+b_cost[2]*pmax(0,train_Y-temp_level)))
        }
        
        price_2_empirical<-function(temp_level_2){
          return(sum(h_cost[2]*pmax(0,temp_level_2-train_Y)+b_cost[2]*pmax(0,train_Y-temp_level_2)))
        }
        
        # res_2<-solnp(pars = c(0.1,rep(0.1,m_n+4)),fun = price_2,LB=c(0,rep(0,m_n+4)),UB=c(2,rep(Inf,m_n+4)),control = list(trace=0))
        res_2 <- optim(par = c(0.1,rep(0.1,m_n+4)), fn = price_2,
                       method = "L-BFGS-B",
                       lower = c(0,rep(0,m_n+4)),
                       upper = c(2,rep(Inf,m_n+4)))
        basis_mat_2<-bSpline(x=train_X+res_2$par[1]*train_Z2,knots = knots_2,degree = 3,Boundary.knots = c(0.9,13.5),intercept = TRUE )
        best_decision_2<-c(basis_mat_2%*%res_2$par[-1])
        
        knots_1<-seq(from=1,to=5,length.out=m_n+2)[-c(1,m_n+2)]
        basis_mat_1<-bSpline(x=train_Z1,knots = knots_1,degree = 3,Boundary.knots = c(1,5),intercept = TRUE )
        price_1<-function(alpha){
          temp_level<-c(basis_mat_1%*%alpha)
          return(sum(h_cost[1]*pmax(0,temp_level-train_X)+b_cost[1]*pmax(0,train_X-temp_level))+price_2_empirical(pmax(temp_level-train_X,best_decision_2)))
        }
        
        # res_1<-solnp(pars = rep(0.1,m_n+4),fun = price_1,LB=rep(0,m_n+4),control = list(trace=0))
        res_1 <- optim(par = rep(0.1,m_n+4), fn = price_1,
                       method = "L-BFGS-B",
                       lower = rep(0,m_n+4))
        
        basis_test_1<-bSpline(x=X_test[1],knots = knots_1,degree = 3,Boundary.knots = c(1,5),intercept = TRUE )
        test_decision_1<-sum(basis_test_1 * res_1$par)
        
        temp_x<-c(X_test[2]*res_2$par[1]+D_test[1])
        temp_base<-bSpline(x=temp_x,knots = knots_2,degree = 3,Boundary.knots = c(0.9,13.5),intercept = TRUE )
        test_decision_2<-sum(temp_base*res_2$par[-1])
      }
      
      current_level <- test_decision_1
      temp_res <- h_cost[1]*max(0,current_level - D_test[1]) + b_cost[1]*max(0,D_test[1]-current_level)
      current_level <- max(current_level - D_test[1], test_decision_2)
      temp_res <- temp_res + h_cost[2]*max(0, current_level - D_test[2])+b_cost[2]*max(0,D_test[2] - current_level)
      cost <- temp_res
      time <- as.numeric((proc.time()-ptm)[3])
      my_online_res[cut_off + 1, ] <- c(cost, time, CVtime)
      
      
      
      basis_test_1<-bSpline(x=X_test[1],knots = opt_knots_1,degree = 3,Boundary.knots = c(1,5),intercept = TRUE )
      test_decision_1<-c(basis_test_1%*%opt_pars1)
      
      temp_x<-c(X_test[2]*opt_pars2[1]+D_test[1])
      temp_base<-bSpline(x=temp_x,knots = opt_knots_2,degree = 3,Boundary.knots = c(0.9,13.5),intercept = TRUE )
      test_decision_2<-c(temp_base%*%opt_pars2[-1])
      
      current_level <- test_decision_1
      temp_res <- h_cost[1]*max(0,current_level - D_test[1]) + b_cost[1]*max(0,D_test[1]-current_level)
      current_level <- max(current_level - D_test[1], test_decision_2)
      temp_res <- temp_res + h_cost[2]*max(0, current_level - D_test[2])+b_cost[2]*max(0,D_test[2] - current_level)
      opt_costs[cut_off + 1] <- temp_res
      
      print(c(cut_off + 1, cost, time, CVtime))
    }
    
    return(list(my_online_res, opt_costs) )
  }
  res<-try(my_simulation_procedure(iter),silent = TRUE)
  return(res)
}


clnum <- 2#detectCores()/2
cl <- makeCluster(getOption("cl.cores", clnum))
res<-parLapply(cl, 1:100,  simulation)
stopCluster(cl)

save(res,file = "online_spline.rdata")
