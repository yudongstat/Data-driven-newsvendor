rm(list=ls())
library(parallel)
iter <- 1
simulation<-function(iter){
  my_simulation_procedure<-function(iter){
    library(splines2)
    library(pracma)
    library(Rsolnp)
    load("Markov_10pcov_200.rdata")
    n<-200
    T<-10
    my_initial <- 0.1
    X_mat <- res[[iter]][[13]]
    D_mat <- res[[iter]][[14]]
    
    load("cov_10_true_Markov.rdata")
    
    b_cost<-rep(1,T)
    b_cost[1]<-5
    h_cost<-rep(1,T)
    h_cost[2]<-5
    full_D_mat <- D_mat
    full_X_mat <- X_mat
    lower_knots <- c(1,0.5,-0.9,0,0.8,1.5,2.5,3.5,4.5,5.5)
    upper_knots <- c(5,15,6,7.5,9,10,11.5,13,14.5,16)
    
    my_final_res <- as.list(rep(NA, n))
    
    for(cut_off in 0 : (n - 1)){
      
      
      knots_1 <- res[[1]][[1]]
      knots_2 <- res[[1]][[2]]
      knots_3 <- res[[1]][[3]]
      knots_4 <- res[[1]][[4]]
      knots_5 <- res[[1]][[5]]
      knots_6 <- res[[1]][[6]]
      knots_7 <- res[[1]][[7]]
      knots_8 <- res[[1]][[8]]
      knots_9 <- res[[1]][[9]]
      knots_10 <- res[[1]][[10]]
      
      pars1 <- res[[2]][[1]]
      pars2 <- res[[2]][[2]]
      pars3 <- res[[2]][[3]]
      pars4 <- res[[2]][[4]]
      pars5 <- res[[2]][[5]]
      pars6 <- res[[2]][[6]]
      pars7 <- res[[2]][[7]]
      pars8 <- res[[2]][[8]]
      pars9 <- res[[2]][[9]]
      pars10 <- res[[2]][[10]]
      
      test_D_mat <- full_D_mat[cut_off + 1, ]
      test_X_mat <- full_X_mat[cut_off + 1, ]
      
      
      basis_test_1<-bSpline(x=test_X_mat[1],knots = knots_1,degree = 3,Boundary.knots = c(lower_knots[1],upper_knots[1]),intercept = TRUE )
      test_decision_1<-c(basis_test_1%*%pars1)
      
      temp_x<-c(test_X_mat[2]*pars2[1]+test_D_mat[1])
      temp_base<-bSpline(x=temp_x,knots = knots_2,degree = 3,Boundary.knots = c(lower_knots[2],upper_knots[2]),intercept = TRUE )
      test_decision_2<-c(temp_base%*%pars2[-1])
      
      temp_x<-c(test_X_mat[3]*pars3[1]+test_D_mat[2])
      temp_base<-bSpline(x=temp_x,knots = knots_3,degree = 3,Boundary.knots = c(lower_knots[3],upper_knots[3]),intercept = TRUE )
      test_decision_3<-c(temp_base%*%pars3[-1])
      
      temp_x<-c(test_X_mat[4]*pars4[1]+test_D_mat[3])
      temp_base<-bSpline(x=temp_x,knots = knots_4,degree = 3,Boundary.knots = c(lower_knots[4],upper_knots[4]),intercept = TRUE )
      test_decision_4<-c(temp_base%*%pars4[-1])
      
      temp_x<-c(test_X_mat[5]*pars5[1]+test_D_mat[4])
      temp_base<-bSpline(x=temp_x,knots = knots_5,degree = 3,Boundary.knots = c(lower_knots[5],upper_knots[5]),intercept = TRUE )
      test_decision_5<-c(temp_base%*%pars5[-1])
      
      temp_x<-c(test_X_mat[6]*pars6[1]+test_D_mat[5])
      temp_base<-bSpline(x=temp_x,knots = knots_6,degree = 3,Boundary.knots = c(lower_knots[6],upper_knots[6]),intercept = TRUE )
      test_decision_6<-c(temp_base%*%pars6[-1])
      
      temp_x<-c(test_X_mat[7]*pars7[1]+test_D_mat[6])
      temp_base<-bSpline(x=temp_x,knots = knots_7,degree = 3,Boundary.knots = c(lower_knots[7],upper_knots[7]),intercept = TRUE )
      test_decision_7<-c(temp_base%*%pars7[-1])
      
      temp_x<-c(test_X_mat[8]*pars8[1]+test_D_mat[7])
      temp_base<-bSpline(x=temp_x,knots = knots_8,degree = 3,Boundary.knots = c(lower_knots[8],upper_knots[8]),intercept = TRUE )
      test_decision_8<-c(temp_base%*%pars8[-1])
      
      temp_x<-c(test_X_mat[9]*pars9[1]+test_D_mat[8])
      temp_base<-bSpline(x=temp_x,knots = knots_9,degree = 3,Boundary.knots = c(lower_knots[9],upper_knots[9]),intercept = TRUE )
      test_decision_9<-c(temp_base%*%pars9[-1])
      
      temp_x<-c(test_X_mat[10]*pars10[1]+test_D_mat[9])
      temp_base<-bSpline(x=temp_x,knots = knots_10,degree = 3,Boundary.knots = c(lower_knots[10],upper_knots[10]),intercept = TRUE )
      test_decision_10<-c(temp_base%*%pars10[-1])
      
      test_decision_oracle <- c(test_decision_1,test_decision_2,test_decision_3,test_decision_4,test_decision_5,test_decision_6,test_decision_7,test_decision_8,test_decision_9,test_decision_10)
      
      start_time <- proc.time()
      
      if(cut_off == 0){
        test_decision_our <- rep(6, T)
      }
      
      if(cut_off == 1){
        test_decision_our <- full_D_mat[1, ]
      }
      
      
      
      if(cut_off > 1){
        m_n_candidate <- 1:10
        D_mat <- full_D_mat[1:cut_off, ]
        X_mat <- full_X_mat[1:cut_off, ]
        if(cut_off < 5){
          m_n <- sample(m_n_candidate, 1)
        }else{
          ########################################################################################################################################################################################################################################
          ########################################################################################################################################################################################################################################
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
            test_D_mat<- D_mat[test_index, ]
            test_X_mat <- X_mat[test_index, ]
            
            train_index <- c(my_split[-fold,])
            train_index <- train_index[which(train_index < 5e3)]
            train_D_mat <- D_mat[train_index, ]
            train_X_mat <- X_mat[train_index, ]
            
            for(ppp in 1:length(m_n_candidate)){
              m_n<-m_n_candidate[ppp]
              
              
              knots_10<-seq(from=lower_knots[10],to=upper_knots[10],length.out=m_n+2)[-c(1,m_n+2)]
              
              price_10<-function(alpha){
                basis_mat_10<-bSpline(x=train_D_mat[,9]+alpha[1]*train_X_mat[,10],knots = knots_10,degree = 3,Boundary.knots = c(lower_knots[10],upper_knots[10]),intercept = TRUE )
                temp_level<-c(basis_mat_10%*%alpha[-1])
                return(sum(h_cost[10]*pmax(0,temp_level-train_D_mat[,10])+b_cost[10]*pmax(0,train_D_mat[,10]-temp_level)))
              }
              
              price_10_empirical<-function(temp_level){
                return(sum(h_cost[10]*pmax(0,temp_level-train_D_mat[,10])+b_cost[10]*pmax(0,train_D_mat[,10]-temp_level)))
              }
              
              price_10_var_est<-function(temp_level){
                return(h_cost[10]*pmax(0,temp_level-train_D_mat[,10])+b_cost[10]*pmax(0,train_D_mat[,10]-temp_level))
              }
              
              res_10<-solnp(pars = c(my_initial,rep(1,m_n+4)),fun = price_10,LB=c(0,rep(0,m_n+4)),UB=c(2,rep(Inf,m_n+4)),control = list(trace=0,outer.iter=10,inner.iter=20))
              
              basis_mat_10<-bSpline(x=train_D_mat[,9]+res_10$par[1]*train_X_mat[,10],knots = knots_10,degree = 3,Boundary.knots = c(lower_knots[10],upper_knots[10]),intercept = TRUE )
              best_decision_10<-c(basis_mat_10%*%res_10$par[-1])
              beta_10<-res_10$par[1]
              
              
              knots_9<-seq(from=lower_knots[9],to=upper_knots[9],length.out=m_n+2)[-c(1,m_n+2)]
              
              price_9<-function(alpha){
                basis_mat_9<-bSpline(x=train_D_mat[,8]+alpha[1]*train_X_mat[,9],knots = knots_9,degree = 3,Boundary.knots = c(lower_knots[9],upper_knots[9]),intercept = TRUE )
                temp_level<-c(basis_mat_9%*%alpha[-1])
                return(sum(h_cost[9]*pmax(0,temp_level-train_D_mat[,9])+b_cost[9]*pmax(0,train_D_mat[,9]-temp_level))+price_10_empirical(pmax(temp_level-train_D_mat[,9],best_decision_10) ))
              }
              
              price_9_empirical<-function(temp_level){
                return(sum(h_cost[9]*pmax(0,temp_level-train_D_mat[,9])+b_cost[9]*pmax(0,train_D_mat[,9]-temp_level))+price_10_empirical(pmax(temp_level-train_D_mat[,9],best_decision_10) ))
              }
              
              price_9_var_est<-function(temp_level){
                return(h_cost[9]*pmax(0,temp_level-train_D_mat[,9])+b_cost[9]*pmax(0,train_D_mat[,9]-temp_level)+price_10_var_est(pmax(temp_level-train_D_mat[,9],best_decision_10) ))
              }
              
              res_9<-solnp(pars = c(my_initial,rep(1,m_n+4)),fun = price_9,LB=c(0,rep(0,m_n+4)),UB=c(2,rep(Inf,m_n+4)),control = list(trace=0,outer.iter=10,inner.iter=20))
              
              basis_mat_9<-bSpline(x=train_D_mat[,8]+res_9$par[1]*train_X_mat[,9],knots = knots_9,degree = 3,Boundary.knots = c(lower_knots[9],upper_knots[9]),intercept = TRUE )
              best_decision_9<-c(basis_mat_9%*%res_9$par[-1])
              beta_9<-res_9$par[1]
              
              
              
              
              knots_8<-seq(from=lower_knots[8],to=upper_knots[8],length.out=m_n+2)[-c(1,m_n+2)]
              
              price_8<-function(alpha){
                basis_mat_8<-bSpline(x=train_D_mat[,7]+alpha[1]*train_X_mat[,8],knots = knots_8,degree = 3,Boundary.knots = c(lower_knots[8],upper_knots[8]),intercept = TRUE )
                temp_level<-c(basis_mat_8%*%alpha[-1])
                return(sum(h_cost[8]*pmax(0,temp_level-train_D_mat[,8])+b_cost[8]*pmax(0,train_D_mat[,8]-temp_level))+price_9_empirical(pmax(temp_level-train_D_mat[,8],best_decision_9) ))
              }
              
              price_8_empirical<-function(temp_level){
                return(sum(h_cost[8]*pmax(0,temp_level-train_D_mat[,8])+b_cost[8]*pmax(0,train_D_mat[,8]-temp_level))+price_9_empirical(pmax(temp_level-train_D_mat[,8],best_decision_9) ))
              }
              
              price_8_var_est<-function(temp_level){
                return(h_cost[8]*pmax(0,temp_level-train_D_mat[,8])+b_cost[8]*pmax(0,train_D_mat[,8]-temp_level)+price_9_var_est(pmax(temp_level-train_D_mat[,8],best_decision_9) ))
              }
              
              res_8<-solnp(pars = c(my_initial,rep(1,m_n+4)),fun = price_8,LB=c(0,rep(0,m_n+4)),UB=c(2,rep(Inf,m_n+4)),control = list(trace=0,outer.iter=10,inner.iter=20))
              
              basis_mat_8<-bSpline(x=train_D_mat[,7]+res_8$par[1]*train_X_mat[,8],knots = knots_8,degree = 3,Boundary.knots = c(lower_knots[8],upper_knots[8]),intercept = TRUE )
              best_decision_8<-c(basis_mat_8%*%res_8$par[-1])
              beta_8<-res_8$par[1]
              
              
              
              knots_7<-seq(from=lower_knots[7],to=upper_knots[7],length.out=m_n+2)[-c(1,m_n+2)]
              
              price_7<-function(alpha){
                basis_mat_7<-bSpline(x=train_D_mat[,6]+alpha[1]*train_X_mat[,7],knots = knots_7,degree = 3,Boundary.knots = c(lower_knots[7],upper_knots[7]),intercept = TRUE )
                temp_level<-c(basis_mat_7%*%alpha[-1])
                return(sum(h_cost[7]*pmax(0,temp_level-train_D_mat[,7])+b_cost[7]*pmax(0,train_D_mat[,7]-temp_level))+price_8_empirical(pmax(temp_level-train_D_mat[,7],best_decision_8) ))
              }
              
              price_7_empirical<-function(temp_level){
                return(sum(h_cost[7]*pmax(0,temp_level-train_D_mat[,7])+b_cost[7]*pmax(0,train_D_mat[,7]-temp_level))+price_8_empirical(pmax(temp_level-train_D_mat[,7],best_decision_8) ))
              }
              
              price_7_var_est<-function(temp_level){
                return(h_cost[7]*pmax(0,temp_level-train_D_mat[,7])+b_cost[7]*pmax(0,train_D_mat[,7]-temp_level)+price_8_var_est(pmax(temp_level-train_D_mat[,7],best_decision_8) ))
              }
              
              res_7<-solnp(pars = c(my_initial,rep(1,m_n+4)),fun = price_7,LB=c(0,rep(0,m_n+4)),UB=c(2,rep(Inf,m_n+4)),control = list(trace=0,outer.iter=10,inner.iter=20))
              
              basis_mat_7<-bSpline(x=train_D_mat[,6]+res_7$par[1]*train_X_mat[,7],knots = knots_7,degree = 3,Boundary.knots = c(lower_knots[7],upper_knots[7]),intercept = TRUE )
              best_decision_7<-c(basis_mat_7%*%res_7$par[-1])
              beta_7<-res_7$par[1]
              
              
              knots_6<-seq(from=lower_knots[6],to=upper_knots[6],length.out=m_n+2)[-c(1,m_n+2)]
              
              price_6<-function(alpha){
                basis_mat_6<-bSpline(x=train_D_mat[,5]+alpha[1]*train_X_mat[,6],knots = knots_6,degree = 3,Boundary.knots = c(lower_knots[6],upper_knots[6]),intercept = TRUE )
                temp_level<-c(basis_mat_6%*%alpha[-1])
                return(sum(h_cost[6]*pmax(0,temp_level-train_D_mat[,6])+b_cost[6]*pmax(0,train_D_mat[,6]-temp_level))+price_7_empirical(pmax(temp_level-train_D_mat[,6],best_decision_7) ))
              }
              
              price_6_empirical<-function(temp_level){
                return(sum(h_cost[6]*pmax(0,temp_level-train_D_mat[,6])+b_cost[6]*pmax(0,train_D_mat[,6]-temp_level))+price_7_empirical(pmax(temp_level-train_D_mat[,6],best_decision_7) ))
              }
              
              price_6_var_est<-function(temp_level){
                return(h_cost[6]*pmax(0,temp_level-train_D_mat[,6])+b_cost[6]*pmax(0,train_D_mat[,6]-temp_level)+price_7_var_est(pmax(temp_level-train_D_mat[,6],best_decision_7) ))
              }
              
              res_6<-solnp(pars = c(my_initial,rep(1,m_n+4)),fun = price_6,LB=c(0,rep(0,m_n+4)),UB=c(2,rep(Inf,m_n+4)),control = list(trace=0,outer.iter=10,inner.iter=20))
              
              basis_mat_6<-bSpline(x=train_D_mat[,5]+res_6$par[1]*train_X_mat[,6],knots = knots_6,degree = 3,Boundary.knots = c(lower_knots[6],upper_knots[6]),intercept = TRUE )
              best_decision_6<-c(basis_mat_6%*%res_6$par[-1])
              beta_6<-res_6$par[1]
              
              
              knots_5<-seq(from=lower_knots[5],to=upper_knots[5],length.out=m_n+2)[-c(1,m_n+2)]
              
              price_5<-function(alpha){
                basis_mat_5<-bSpline(x=train_D_mat[,4]+alpha[1]*train_X_mat[,5],knots = knots_5,degree = 3,Boundary.knots = c(lower_knots[5],upper_knots[5]),intercept = TRUE )
                temp_level<-c(basis_mat_5%*%alpha[-1])
                return(sum(h_cost[5]*pmax(0,temp_level-train_D_mat[,5])+b_cost[5]*pmax(0,train_D_mat[,5]-temp_level))+price_6_empirical(pmax(temp_level-train_D_mat[,5],best_decision_6) ))
              }
              
              price_5_empirical<-function(temp_level){
                return(sum(h_cost[5]*pmax(0,temp_level-train_D_mat[,5])+b_cost[5]*pmax(0,train_D_mat[,5]-temp_level))+price_6_empirical(pmax(temp_level-train_D_mat[,5],best_decision_6) ))
              }
              
              price_5_var_est<-function(temp_level){
                return(h_cost[5]*pmax(0,temp_level-train_D_mat[,5])+b_cost[5]*pmax(0,train_D_mat[,5]-temp_level)+price_6_var_est(pmax(temp_level-train_D_mat[,5],best_decision_6) ))
              }
              
              res_5<-solnp(pars = c(my_initial,rep(1,m_n+4)),fun = price_5,LB=c(0,rep(0,m_n+4)),UB=c(2,rep(Inf,m_n+4)),control = list(trace=0,outer.iter=10,inner.iter=20))
              
              basis_mat_5<-bSpline(x=train_D_mat[,4]+res_5$par[1]*train_X_mat[,5],knots = knots_5,degree = 3,Boundary.knots = c(lower_knots[5],upper_knots[5]),intercept = TRUE )
              best_decision_5<-c(basis_mat_5%*%res_5$par[-1])
              beta_5<-res_5$par[1]
              
              knots_4<-seq(from=lower_knots[4],to=upper_knots[4],length.out=m_n+2)[-c(1,m_n+2)]
              
              price_4<-function(alpha){
                basis_mat_4<-bSpline(x=train_D_mat[,3]+alpha[1]*train_X_mat[,4],knots = knots_4,degree = 3,Boundary.knots = c(lower_knots[4],upper_knots[4]),intercept = TRUE )
                temp_level<-c(basis_mat_4%*%alpha[-1])
                return(sum(h_cost[4]*pmax(0,temp_level-train_D_mat[,4])+b_cost[4]*pmax(0,train_D_mat[,4]-temp_level))+price_5_empirical(pmax(temp_level-train_D_mat[,4],best_decision_5) ))
              }
              
              price_4_empirical<-function(temp_level){
                return(sum(h_cost[4]*pmax(0,temp_level-train_D_mat[,4])+b_cost[4]*pmax(0,train_D_mat[,4]-temp_level))+price_5_empirical(pmax(temp_level-train_D_mat[,4],best_decision_5) ))
              }
              
              price_4_var_est<-function(temp_level){
                return(h_cost[4]*pmax(0,temp_level-train_D_mat[,4])+b_cost[4]*pmax(0,train_D_mat[,4]-temp_level)+price_5_var_est(pmax(temp_level-train_D_mat[,4],best_decision_5) ))
              }
              
              res_4<-solnp(pars = c(my_initial,rep(1,m_n+4)),fun = price_4,LB=c(0,rep(0,m_n+4)),UB=c(2,rep(Inf,m_n+4)),control = list(trace=0,outer.iter=10,inner.iter=20))
              
              basis_mat_4<-bSpline(x=train_D_mat[,3]+res_4$par[1]*train_X_mat[,4],knots = knots_4,degree = 3,Boundary.knots = c(lower_knots[4],upper_knots[4]),intercept = TRUE )
              best_decision_4<-c(basis_mat_4%*%res_4$par[-1])
              beta_4<-res_4$par[1]
              
              knots_3<-seq(from=lower_knots[3],to=upper_knots[3],length.out=m_n+2)[-c(1,m_n+2)]
              
              price_3<-function(alpha){
                basis_mat_3<-bSpline(x=train_D_mat[,2]+alpha[1]*train_X_mat[,3],knots = knots_3,degree = 3,Boundary.knots = c(lower_knots[3],upper_knots[3]),intercept = TRUE )
                temp_level<-c(basis_mat_3%*%alpha[-1])
                return(sum(h_cost[3]*pmax(0,temp_level-train_D_mat[,3])+b_cost[3]*pmax(0,train_D_mat[,3]-temp_level))+price_4_empirical(pmax(temp_level-train_D_mat[,3],best_decision_4) ))
              }
              
              price_3_empirical<-function(temp_level){
                return(sum(h_cost[3]*pmax(0,temp_level-train_D_mat[,3])+b_cost[3]*pmax(0,train_D_mat[,3]-temp_level))+price_4_empirical(pmax(temp_level-train_D_mat[,3],best_decision_4) ))
              }
              
              price_3_var_est<-function(temp_level){
                return(h_cost[3]*pmax(0,temp_level-train_D_mat[,3])+b_cost[3]*pmax(0,train_D_mat[,3]-temp_level)+price_4_var_est(pmax(temp_level-train_D_mat[,3],best_decision_4) ))
              }
              
              res_3<-solnp(pars = c(my_initial,rep(1,m_n+4)),fun = price_3,LB=c(0,rep(0,m_n+4)),UB=c(2,rep(Inf,m_n+4)),control = list(trace=0,outer.iter=10,inner.iter=20))
              
              basis_mat_3<-bSpline(x=train_D_mat[,2]+res_3$par[1]*train_X_mat[,3],knots = knots_3,degree = 3,Boundary.knots = c(lower_knots[3],upper_knots[3]),intercept = TRUE )
              best_decision_3<-c(basis_mat_3%*%res_3$par[-1])
              beta_3<-res_3$par[1]
              
              
              
              knots_2<-seq(from=lower_knots[2],to=upper_knots[2],length.out=m_n+2)[-c(1,m_n+2)]
              
              price_2<-function(alpha){
                basis_mat_2<-bSpline(x=train_D_mat[,1]+alpha[1]*train_X_mat[,2],knots = knots_2,degree = 3,Boundary.knots = c(lower_knots[2],upper_knots[2]),intercept = TRUE )
                temp_level<-c(basis_mat_2%*%alpha[-1])
                return(sum(h_cost[2]*pmax(0,temp_level-train_D_mat[,2])+b_cost[2]*pmax(0,train_D_mat[,2]-temp_level))+price_3_empirical(pmax(temp_level-train_D_mat[,2],best_decision_3) ))
              }
              
              price_2_empirical<-function(temp_level_2){
                return(sum(h_cost[2]*pmax(0,temp_level_2-train_D_mat[,2])+b_cost[2]*pmax(0,train_D_mat[,2]-temp_level_2))+price_3_empirical(pmax(temp_level_2-train_D_mat[,2],best_decision_3) ))
              }
              
              
              price_2_var_est<-function(temp_level_2){
                return(h_cost[2]*pmax(0,temp_level_2-train_D_mat[,2])+b_cost[2]*pmax(0,train_D_mat[,2]-temp_level_2)+price_3_var_est(pmax(temp_level_2-train_D_mat[,2],best_decision_3) ))
              }
              
              res_2<-solnp(pars = c(my_initial,rep(1,m_n+4)),fun = price_2,LB=c(0,rep(0,m_n+4)),UB=c(2,rep(Inf,m_n+4)),control = list(trace=0,outer.iter=10,inner.iter=20))
              
              basis_mat_2<-bSpline(x=train_D_mat[,1]+res_2$par[1]*train_X_mat[,2],knots = knots_2,degree = 3,Boundary.knots = c(lower_knots[2],upper_knots[2]),intercept = TRUE )
              best_decision_2<-c(basis_mat_2%*%res_2$par[-1])
              beta_2<-res_2$par[1]
              
              
              
              knots_1<-seq(from=lower_knots[1],to=upper_knots[1],length.out=m_n+2)[-c(1,m_n+2)]
              basis_mat_1<-bSpline(x=train_X_mat[,1],knots = knots_1,degree = 3,Boundary.knots = c(lower_knots[1],upper_knots[1]),intercept = TRUE )
              price_1<-function(alpha){
                temp_level<-c(basis_mat_1%*%alpha)
                return(sum(h_cost[1]*pmax(0,temp_level-train_D_mat[,1])+b_cost[1]*pmax(0,train_D_mat[,1]-temp_level))+price_2_empirical(pmax(temp_level-train_D_mat[,1],best_decision_2)))
              }
              
              price_1_var_est<-function(alpha){
                temp_level<-c(basis_mat_1%*%alpha)
                return(h_cost[1]*pmax(0,temp_level-train_D_mat[,1])+b_cost[1]*pmax(0,train_D_mat[,1]-temp_level)+price_2_var_est(pmax(temp_level-train_D_mat[,1],best_decision_2)))
              }
              
              res_1<-solnp(pars = rep(1,m_n+4),fun = price_1,LB=rep(0,m_n+4),control = list(trace=0,outer.iter=10,inner.iter=20))
              
              if(is.null(nrow(test_D_mat))){
                basis_test_1<-bSpline(x=test_X_mat[1],knots = knots_1,degree = 3,Boundary.knots = c(lower_knots[1],upper_knots[1]),intercept = TRUE )
                test_decision_1<-c(basis_test_1%*%res_1$par)
                
                temp_x<-c(test_X_mat[2]*res_2$par[1]+test_D_mat[1])
                temp_base<-bSpline(x=temp_x,knots = knots_2,degree = 3,Boundary.knots = c(lower_knots[2],upper_knots[2]),intercept = TRUE )
                test_decision_2<-c(temp_base%*%res_2$par[-1])
                
                temp_x<-c(test_X_mat[3]*res_3$par[1]+test_D_mat[2])
                temp_base<-bSpline(x=temp_x,knots = knots_3,degree = 3,Boundary.knots = c(lower_knots[3],upper_knots[3]),intercept = TRUE )
                test_decision_3<-c(temp_base%*%res_3$par[-1])
                
                temp_x<-c(test_X_mat[4]*res_4$par[1]+test_D_mat[3])
                temp_base<-bSpline(x=temp_x,knots = knots_4,degree = 3,Boundary.knots = c(lower_knots[4],upper_knots[4]),intercept = TRUE )
                test_decision_4<-c(temp_base%*%res_4$par[-1])
                
                temp_x<-c(test_X_mat[5]*res_5$par[1]+test_D_mat[4])
                temp_base<-bSpline(x=temp_x,knots = knots_5,degree = 3,Boundary.knots = c(lower_knots[5],upper_knots[5]),intercept = TRUE )
                test_decision_5<-c(temp_base%*%res_5$par[-1])
                
                temp_x<-c(test_X_mat[6]*res_6$par[1]+test_D_mat[5])
                temp_base<-bSpline(x=temp_x,knots = knots_6,degree = 3,Boundary.knots = c(lower_knots[6],upper_knots[6]),intercept = TRUE )
                test_decision_6<-c(temp_base%*%res_6$par[-1])
                
                temp_x<-c(test_X_mat[7]*res_7$par[1]+test_D_mat[6])
                temp_base<-bSpline(x=temp_x,knots = knots_7,degree = 3,Boundary.knots = c(lower_knots[7],upper_knots[7]),intercept = TRUE )
                test_decision_7<-c(temp_base%*%res_7$par[-1])
                
                temp_x<-c(test_X_mat[8]*res_8$par[1]+test_D_mat[7])
                temp_base<-bSpline(x=temp_x,knots = knots_8,degree = 3,Boundary.knots = c(lower_knots[8],upper_knots[8]),intercept = TRUE )
                test_decision_8<-c(temp_base%*%res_8$par[-1])
                
                temp_x<-c(test_X_mat[9]*res_9$par[1]+test_D_mat[8])
                temp_base<-bSpline(x=temp_x,knots = knots_9,degree = 3,Boundary.knots = c(lower_knots[9],upper_knots[9]),intercept = TRUE )
                test_decision_9<-c(temp_base%*%res_9$par[-1])
                
                temp_x<-c(test_X_mat[10]*res_10$par[1]+test_D_mat[9])
                temp_base<-bSpline(x=temp_x,knots = knots_10,degree = 3,Boundary.knots = c(lower_knots[10],upper_knots[10]),intercept = TRUE )
                test_decision_10<-c(temp_base%*%res_10$par[-1])
                
                test_decision <- c(test_decision_2,test_decision_3,test_decision_4,test_decision_5,test_decision_6,test_decision_7,test_decision_8,test_decision_9,test_decision_10)
                
                temp_res<-0
                temp_level<-test_decision_1
                for(qqqq in 1:T){
                  temp_res<-temp_res+h_cost[qqqq]*pmax(0,temp_level-test_D_mat[qqqq])+b_cost[qqqq]*pmax(0,test_D_mat[qqqq]-temp_level)
                  # print(temp_res)
                  if(qqqq<T){
                    temp_level<-pmax(temp_level-test_D_mat[qqqq],test_decision[qqqq])
                  }
                }
                CV_mat[fold,ppp]<-temp_res
                
              }else{
                basis_test_1<-bSpline(x=test_X_mat[,1],knots = knots_1,degree = 3,Boundary.knots = c(lower_knots[1],upper_knots[1]),intercept = TRUE )
                test_decision_1<-c(basis_test_1%*%res_1$par)
                
                temp_x<-c(test_X_mat[,2]*res_2$par[1]+test_D_mat[,1])
                temp_base<-bSpline(x=temp_x,knots = knots_2,degree = 3,Boundary.knots = c(lower_knots[2],upper_knots[2]),intercept = TRUE )
                test_decision_2<-c(temp_base%*%res_2$par[-1])
                
                temp_x<-c(test_X_mat[,3]*res_3$par[1]+test_D_mat[,2])
                temp_base<-bSpline(x=temp_x,knots = knots_3,degree = 3,Boundary.knots = c(lower_knots[3],upper_knots[3]),intercept = TRUE )
                test_decision_3<-c(temp_base%*%res_3$par[-1])
                
                temp_x<-c(test_X_mat[,4]*res_4$par[1]+test_D_mat[,3])
                temp_base<-bSpline(x=temp_x,knots = knots_4,degree = 3,Boundary.knots = c(lower_knots[4],upper_knots[4]),intercept = TRUE )
                test_decision_4<-c(temp_base%*%res_4$par[-1])
                
                temp_x<-c(test_X_mat[,5]*res_5$par[1]+test_D_mat[,4])
                temp_base<-bSpline(x=temp_x,knots = knots_5,degree = 3,Boundary.knots = c(lower_knots[5],upper_knots[5]),intercept = TRUE )
                test_decision_5<-c(temp_base%*%res_5$par[-1])
                
                temp_x<-c(test_X_mat[,6]*res_6$par[1]+test_D_mat[,5])
                temp_base<-bSpline(x=temp_x,knots = knots_6,degree = 3,Boundary.knots = c(lower_knots[6],upper_knots[6]),intercept = TRUE )
                test_decision_6<-c(temp_base%*%res_6$par[-1])
                
                temp_x<-c(test_X_mat[,7]*res_7$par[1]+test_D_mat[,6])
                temp_base<-bSpline(x=temp_x,knots = knots_7,degree = 3,Boundary.knots = c(lower_knots[7],upper_knots[7]),intercept = TRUE )
                test_decision_7<-c(temp_base%*%res_7$par[-1])
                
                temp_x<-c(test_X_mat[,8]*res_8$par[1]+test_D_mat[,7])
                temp_base<-bSpline(x=temp_x,knots = knots_8,degree = 3,Boundary.knots = c(lower_knots[8],upper_knots[8]),intercept = TRUE )
                test_decision_8<-c(temp_base%*%res_8$par[-1])
                
                temp_x<-c(test_X_mat[,9]*res_9$par[1]+test_D_mat[,8])
                temp_base<-bSpline(x=temp_x,knots = knots_9,degree = 3,Boundary.knots = c(lower_knots[9],upper_knots[9]),intercept = TRUE )
                test_decision_9<-c(temp_base%*%res_9$par[-1])
                
                temp_x<-c(test_X_mat[,10]*res_10$par[1]+test_D_mat[,9])
                temp_base<-bSpline(x=temp_x,knots = knots_10,degree = 3,Boundary.knots = c(lower_knots[10],upper_knots[10]),intercept = TRUE )
                test_decision_10<-c(temp_base%*%res_10$par[-1])
                
                
                test_decision<-rbind(test_decision_2,test_decision_3,test_decision_4,test_decision_5,test_decision_6,test_decision_7,test_decision_8,test_decision_9,test_decision_10)
                
                temp_res<-0
                temp_level<-test_decision_1
                for(qqqq in 1:T){
                  temp_res<-temp_res+sum(h_cost[qqqq]*pmax(0,temp_level-test_D_mat[,qqqq])+b_cost[qqqq]*pmax(0,test_D_mat[,qqqq]-temp_level))
                  # print(temp_res)
                  if(qqqq<T){
                    temp_level<-pmax(temp_level-test_D_mat[,qqqq],test_decision[qqqq,]) 
                  }
                }
                CV_mat[fold,ppp]<-temp_res
              }
              # print(c(fold,ppp))
            }
          }
          CV_seq<-apply(CV_mat,MARGIN = 2,sum)
          
          m_n<-m_n_candidate[(which.min(CV_seq))[1]]
        }
        
        train_D_mat <- full_D_mat[1:cut_off, ]
        train_X_mat <- full_X_mat[1:cut_off, ]
        
        
        
        knots_10<-seq(from=lower_knots[10],to=upper_knots[10],length.out=m_n+2)[-c(1,m_n+2)]
        
        price_10<-function(alpha){
          basis_mat_10<-bSpline(x=train_D_mat[,9]+alpha[1]*train_X_mat[,10],knots = knots_10,degree = 3,Boundary.knots = c(lower_knots[10],upper_knots[10]),intercept = TRUE )
          temp_level<-c(basis_mat_10%*%alpha[-1])
          return(sum(h_cost[10]*pmax(0,temp_level-train_D_mat[,10])+b_cost[10]*pmax(0,train_D_mat[,10]-temp_level)))
        }
        
        price_10_empirical<-function(temp_level){
          return(sum(h_cost[10]*pmax(0,temp_level-train_D_mat[,10])+b_cost[10]*pmax(0,train_D_mat[,10]-temp_level)))
        }
        
        price_10_var_est<-function(temp_level){
          return(h_cost[10]*pmax(0,temp_level-train_D_mat[,10])+b_cost[10]*pmax(0,train_D_mat[,10]-temp_level))
        }
        
        res_10<-solnp(pars = c(my_initial,rep(1,m_n+4)),fun = price_10,LB=c(0,rep(0,m_n+4)),UB=c(2,rep(Inf,m_n+4)),control = list(trace=0,outer.iter=10,inner.iter=20))
        
        basis_mat_10<-bSpline(x=train_D_mat[,9]+res_10$par[1]*train_X_mat[,10],knots = knots_10,degree = 3,Boundary.knots = c(lower_knots[10],upper_knots[10]),intercept = TRUE )
        best_decision_10<-c(basis_mat_10%*%res_10$par[-1])
        beta_10<-res_10$par[1]
        
        
        knots_9<-seq(from=lower_knots[9],to=upper_knots[9],length.out=m_n+2)[-c(1,m_n+2)]
        
        price_9<-function(alpha){
          basis_mat_9<-bSpline(x=train_D_mat[,8]+alpha[1]*train_X_mat[,9],knots = knots_9,degree = 3,Boundary.knots = c(lower_knots[9],upper_knots[9]),intercept = TRUE )
          temp_level<-c(basis_mat_9%*%alpha[-1])
          return(sum(h_cost[9]*pmax(0,temp_level-train_D_mat[,9])+b_cost[9]*pmax(0,train_D_mat[,9]-temp_level))+price_10_empirical(pmax(temp_level-train_D_mat[,9],best_decision_10) ))
        }
        
        price_9_empirical<-function(temp_level){
          return(sum(h_cost[9]*pmax(0,temp_level-train_D_mat[,9])+b_cost[9]*pmax(0,train_D_mat[,9]-temp_level))+price_10_empirical(pmax(temp_level-train_D_mat[,9],best_decision_10) ))
        }
        
        price_9_var_est<-function(temp_level){
          return(h_cost[9]*pmax(0,temp_level-train_D_mat[,9])+b_cost[9]*pmax(0,train_D_mat[,9]-temp_level)+price_10_var_est(pmax(temp_level-train_D_mat[,9],best_decision_10) ))
        }
        
        res_9<-solnp(pars = c(my_initial,rep(1,m_n+4)),fun = price_9,LB=c(0,rep(0,m_n+4)),UB=c(2,rep(Inf,m_n+4)),control = list(trace=0,outer.iter=10,inner.iter=20))
        
        basis_mat_9<-bSpline(x=train_D_mat[,8]+res_9$par[1]*train_X_mat[,9],knots = knots_9,degree = 3,Boundary.knots = c(lower_knots[9],upper_knots[9]),intercept = TRUE )
        best_decision_9<-c(basis_mat_9%*%res_9$par[-1])
        beta_9<-res_9$par[1]
        
        
        
        
        knots_8<-seq(from=lower_knots[8],to=upper_knots[8],length.out=m_n+2)[-c(1,m_n+2)]
        
        price_8<-function(alpha){
          basis_mat_8<-bSpline(x=train_D_mat[,7]+alpha[1]*train_X_mat[,8],knots = knots_8,degree = 3,Boundary.knots = c(lower_knots[8],upper_knots[8]),intercept = TRUE )
          temp_level<-c(basis_mat_8%*%alpha[-1])
          return(sum(h_cost[8]*pmax(0,temp_level-train_D_mat[,8])+b_cost[8]*pmax(0,train_D_mat[,8]-temp_level))+price_9_empirical(pmax(temp_level-train_D_mat[,8],best_decision_9) ))
        }
        
        price_8_empirical<-function(temp_level){
          return(sum(h_cost[8]*pmax(0,temp_level-train_D_mat[,8])+b_cost[8]*pmax(0,train_D_mat[,8]-temp_level))+price_9_empirical(pmax(temp_level-train_D_mat[,8],best_decision_9) ))
        }
        
        price_8_var_est<-function(temp_level){
          return(h_cost[8]*pmax(0,temp_level-train_D_mat[,8])+b_cost[8]*pmax(0,train_D_mat[,8]-temp_level)+price_9_var_est(pmax(temp_level-train_D_mat[,8],best_decision_9) ))
        }
        
        res_8<-solnp(pars = c(my_initial,rep(1,m_n+4)),fun = price_8,LB=c(0,rep(0,m_n+4)),UB=c(2,rep(Inf,m_n+4)),control = list(trace=0,outer.iter=10,inner.iter=20))
        
        basis_mat_8<-bSpline(x=train_D_mat[,7]+res_8$par[1]*train_X_mat[,8],knots = knots_8,degree = 3,Boundary.knots = c(lower_knots[8],upper_knots[8]),intercept = TRUE )
        best_decision_8<-c(basis_mat_8%*%res_8$par[-1])
        beta_8<-res_8$par[1]
        
        
        
        knots_7<-seq(from=lower_knots[7],to=upper_knots[7],length.out=m_n+2)[-c(1,m_n+2)]
        
        price_7<-function(alpha){
          basis_mat_7<-bSpline(x=train_D_mat[,6]+alpha[1]*train_X_mat[,7],knots = knots_7,degree = 3,Boundary.knots = c(lower_knots[7],upper_knots[7]),intercept = TRUE )
          temp_level<-c(basis_mat_7%*%alpha[-1])
          return(sum(h_cost[7]*pmax(0,temp_level-train_D_mat[,7])+b_cost[7]*pmax(0,train_D_mat[,7]-temp_level))+price_8_empirical(pmax(temp_level-train_D_mat[,7],best_decision_8) ))
        }
        
        price_7_empirical<-function(temp_level){
          return(sum(h_cost[7]*pmax(0,temp_level-train_D_mat[,7])+b_cost[7]*pmax(0,train_D_mat[,7]-temp_level))+price_8_empirical(pmax(temp_level-train_D_mat[,7],best_decision_8) ))
        }
        
        price_7_var_est<-function(temp_level){
          return(h_cost[7]*pmax(0,temp_level-train_D_mat[,7])+b_cost[7]*pmax(0,train_D_mat[,7]-temp_level)+price_8_var_est(pmax(temp_level-train_D_mat[,7],best_decision_8) ))
        }
        
        res_7<-solnp(pars = c(my_initial,rep(1,m_n+4)),fun = price_7,LB=c(0,rep(0,m_n+4)),UB=c(2,rep(Inf,m_n+4)),control = list(trace=0,outer.iter=10,inner.iter=20))
        
        basis_mat_7<-bSpline(x=train_D_mat[,6]+res_7$par[1]*train_X_mat[,7],knots = knots_7,degree = 3,Boundary.knots = c(lower_knots[7],upper_knots[7]),intercept = TRUE )
        best_decision_7<-c(basis_mat_7%*%res_7$par[-1])
        beta_7<-res_7$par[1]
        
        
        knots_6<-seq(from=lower_knots[6],to=upper_knots[6],length.out=m_n+2)[-c(1,m_n+2)]
        
        price_6<-function(alpha){
          basis_mat_6<-bSpline(x=train_D_mat[,5]+alpha[1]*train_X_mat[,6],knots = knots_6,degree = 3,Boundary.knots = c(lower_knots[6],upper_knots[6]),intercept = TRUE )
          temp_level<-c(basis_mat_6%*%alpha[-1])
          return(sum(h_cost[6]*pmax(0,temp_level-train_D_mat[,6])+b_cost[6]*pmax(0,train_D_mat[,6]-temp_level))+price_7_empirical(pmax(temp_level-train_D_mat[,6],best_decision_7) ))
        }
        
        price_6_empirical<-function(temp_level){
          return(sum(h_cost[6]*pmax(0,temp_level-train_D_mat[,6])+b_cost[6]*pmax(0,train_D_mat[,6]-temp_level))+price_7_empirical(pmax(temp_level-train_D_mat[,6],best_decision_7) ))
        }
        
        price_6_var_est<-function(temp_level){
          return(h_cost[6]*pmax(0,temp_level-train_D_mat[,6])+b_cost[6]*pmax(0,train_D_mat[,6]-temp_level)+price_7_var_est(pmax(temp_level-train_D_mat[,6],best_decision_7) ))
        }
        
        res_6<-solnp(pars = c(my_initial,rep(1,m_n+4)),fun = price_6,LB=c(0,rep(0,m_n+4)),UB=c(2,rep(Inf,m_n+4)),control = list(trace=0,outer.iter=10,inner.iter=20))
        
        basis_mat_6<-bSpline(x=train_D_mat[,5]+res_6$par[1]*train_X_mat[,6],knots = knots_6,degree = 3,Boundary.knots = c(lower_knots[6],upper_knots[6]),intercept = TRUE )
        best_decision_6<-c(basis_mat_6%*%res_6$par[-1])
        beta_6<-res_6$par[1]
        
        
        knots_5<-seq(from=lower_knots[5],to=upper_knots[5],length.out=m_n+2)[-c(1,m_n+2)]
        
        price_5<-function(alpha){
          basis_mat_5<-bSpline(x=train_D_mat[,4]+alpha[1]*train_X_mat[,5],knots = knots_5,degree = 3,Boundary.knots = c(lower_knots[5],upper_knots[5]),intercept = TRUE )
          temp_level<-c(basis_mat_5%*%alpha[-1])
          return(sum(h_cost[5]*pmax(0,temp_level-train_D_mat[,5])+b_cost[5]*pmax(0,train_D_mat[,5]-temp_level))+price_6_empirical(pmax(temp_level-train_D_mat[,5],best_decision_6) ))
        }
        
        price_5_empirical<-function(temp_level){
          return(sum(h_cost[5]*pmax(0,temp_level-train_D_mat[,5])+b_cost[5]*pmax(0,train_D_mat[,5]-temp_level))+price_6_empirical(pmax(temp_level-train_D_mat[,5],best_decision_6) ))
        }
        
        price_5_var_est<-function(temp_level){
          return(h_cost[5]*pmax(0,temp_level-train_D_mat[,5])+b_cost[5]*pmax(0,train_D_mat[,5]-temp_level)+price_6_var_est(pmax(temp_level-train_D_mat[,5],best_decision_6) ))
        }
        
        res_5<-solnp(pars = c(my_initial,rep(1,m_n+4)),fun = price_5,LB=c(0,rep(0,m_n+4)),UB=c(2,rep(Inf,m_n+4)),control = list(trace=0,outer.iter=10,inner.iter=20))
        
        basis_mat_5<-bSpline(x=train_D_mat[,4]+res_5$par[1]*train_X_mat[,5],knots = knots_5,degree = 3,Boundary.knots = c(lower_knots[5],upper_knots[5]),intercept = TRUE )
        best_decision_5<-c(basis_mat_5%*%res_5$par[-1])
        beta_5<-res_5$par[1]
        
        knots_4<-seq(from=lower_knots[4],to=upper_knots[4],length.out=m_n+2)[-c(1,m_n+2)]
        
        price_4<-function(alpha){
          basis_mat_4<-bSpline(x=train_D_mat[,3]+alpha[1]*train_X_mat[,4],knots = knots_4,degree = 3,Boundary.knots = c(lower_knots[4],upper_knots[4]),intercept = TRUE )
          temp_level<-c(basis_mat_4%*%alpha[-1])
          return(sum(h_cost[4]*pmax(0,temp_level-train_D_mat[,4])+b_cost[4]*pmax(0,train_D_mat[,4]-temp_level))+price_5_empirical(pmax(temp_level-train_D_mat[,4],best_decision_5) ))
        }
        
        price_4_empirical<-function(temp_level){
          return(sum(h_cost[4]*pmax(0,temp_level-train_D_mat[,4])+b_cost[4]*pmax(0,train_D_mat[,4]-temp_level))+price_5_empirical(pmax(temp_level-train_D_mat[,4],best_decision_5) ))
        }
        
        price_4_var_est<-function(temp_level){
          return(h_cost[4]*pmax(0,temp_level-train_D_mat[,4])+b_cost[4]*pmax(0,train_D_mat[,4]-temp_level)+price_5_var_est(pmax(temp_level-train_D_mat[,4],best_decision_5) ))
        }
        
        res_4<-solnp(pars = c(my_initial,rep(1,m_n+4)),fun = price_4,LB=c(0,rep(0,m_n+4)),UB=c(2,rep(Inf,m_n+4)),control = list(trace=0,outer.iter=10,inner.iter=20))
        
        basis_mat_4<-bSpline(x=train_D_mat[,3]+res_4$par[1]*train_X_mat[,4],knots = knots_4,degree = 3,Boundary.knots = c(lower_knots[4],upper_knots[4]),intercept = TRUE )
        best_decision_4<-c(basis_mat_4%*%res_4$par[-1])
        beta_4<-res_4$par[1]
        
        knots_3<-seq(from=lower_knots[3],to=upper_knots[3],length.out=m_n+2)[-c(1,m_n+2)]
        
        price_3<-function(alpha){
          basis_mat_3<-bSpline(x=train_D_mat[,2]+alpha[1]*train_X_mat[,3],knots = knots_3,degree = 3,Boundary.knots = c(lower_knots[3],upper_knots[3]),intercept = TRUE )
          temp_level<-c(basis_mat_3%*%alpha[-1])
          return(sum(h_cost[3]*pmax(0,temp_level-train_D_mat[,3])+b_cost[3]*pmax(0,train_D_mat[,3]-temp_level))+price_4_empirical(pmax(temp_level-train_D_mat[,3],best_decision_4) ))
        }
        
        price_3_empirical<-function(temp_level){
          return(sum(h_cost[3]*pmax(0,temp_level-train_D_mat[,3])+b_cost[3]*pmax(0,train_D_mat[,3]-temp_level))+price_4_empirical(pmax(temp_level-train_D_mat[,3],best_decision_4) ))
        }
        
        price_3_var_est<-function(temp_level){
          return(h_cost[3]*pmax(0,temp_level-train_D_mat[,3])+b_cost[3]*pmax(0,train_D_mat[,3]-temp_level)+price_4_var_est(pmax(temp_level-train_D_mat[,3],best_decision_4) ))
        }
        
        res_3<-solnp(pars = c(my_initial,rep(1,m_n+4)),fun = price_3,LB=c(0,rep(0,m_n+4)),UB=c(2,rep(Inf,m_n+4)),control = list(trace=0,outer.iter=10,inner.iter=20))
        
        basis_mat_3<-bSpline(x=train_D_mat[,2]+res_3$par[1]*train_X_mat[,3],knots = knots_3,degree = 3,Boundary.knots = c(lower_knots[3],upper_knots[3]),intercept = TRUE )
        best_decision_3<-c(basis_mat_3%*%res_3$par[-1])
        beta_3<-res_3$par[1]
        
        
        
        knots_2<-seq(from=lower_knots[2],to=upper_knots[2],length.out=m_n+2)[-c(1,m_n+2)]
        
        price_2<-function(alpha){
          basis_mat_2<-bSpline(x=train_D_mat[,1]+alpha[1]*train_X_mat[,2],knots = knots_2,degree = 3,Boundary.knots = c(lower_knots[2],upper_knots[2]),intercept = TRUE )
          temp_level<-c(basis_mat_2%*%alpha[-1])
          return(sum(h_cost[2]*pmax(0,temp_level-train_D_mat[,2])+b_cost[2]*pmax(0,train_D_mat[,2]-temp_level))+price_3_empirical(pmax(temp_level-train_D_mat[,2],best_decision_3) ))
        }
        
        price_2_empirical<-function(temp_level_2){
          return(sum(h_cost[2]*pmax(0,temp_level_2-train_D_mat[,2])+b_cost[2]*pmax(0,train_D_mat[,2]-temp_level_2))+price_3_empirical(pmax(temp_level_2-train_D_mat[,2],best_decision_3) ))
        }
        
        
        price_2_var_est<-function(temp_level_2){
          return(h_cost[2]*pmax(0,temp_level_2-train_D_mat[,2])+b_cost[2]*pmax(0,train_D_mat[,2]-temp_level_2)+price_3_var_est(pmax(temp_level_2-train_D_mat[,2],best_decision_3) ))
        }
        
        res_2<-solnp(pars = c(my_initial,rep(1,m_n+4)),fun = price_2,LB=c(0,rep(0,m_n+4)),UB=c(2,rep(Inf,m_n+4)),control = list(trace=0,outer.iter=10,inner.iter=20))
        
        basis_mat_2<-bSpline(x=train_D_mat[,1]+res_2$par[1]*train_X_mat[,2],knots = knots_2,degree = 3,Boundary.knots = c(lower_knots[2],upper_knots[2]),intercept = TRUE )
        best_decision_2<-c(basis_mat_2%*%res_2$par[-1])
        beta_2<-res_2$par[1]
        
        
        
        knots_1<-seq(from=lower_knots[1],to=upper_knots[1],length.out=m_n+2)[-c(1,m_n+2)]
        basis_mat_1<-bSpline(x=train_X_mat[,1],knots = knots_1,degree = 3,Boundary.knots = c(lower_knots[1],upper_knots[1]),intercept = TRUE )
        price_1<-function(alpha){
          temp_level<-c(basis_mat_1%*%alpha)
          return(sum(h_cost[1]*pmax(0,temp_level-train_D_mat[,1])+b_cost[1]*pmax(0,train_D_mat[,1]-temp_level))+price_2_empirical(pmax(temp_level-train_D_mat[,1],best_decision_2)))
        }
        
        price_1_var_est<-function(alpha){
          temp_level<-c(basis_mat_1%*%alpha)
          return(h_cost[1]*pmax(0,temp_level-train_D_mat[,1])+b_cost[1]*pmax(0,train_D_mat[,1]-temp_level)+price_2_var_est(pmax(temp_level-train_D_mat[,1],best_decision_2)))
        }
        
        res_1<-solnp(pars = rep(1,m_n+4),fun = price_1,LB=rep(0,m_n+4),control = list(trace=0,outer.iter=10,inner.iter=20))
        
        
        
        test_D_mat <- full_D_mat[cut_off + 1, ]
        test_X_mat <- full_X_mat[cut_off + 1, ]
        
        basis_test_1<-bSpline(x=test_X_mat[1],knots = knots_1,degree = 3,Boundary.knots = c(lower_knots[1],upper_knots[1]),intercept = TRUE )
        test_decision_1<-c(basis_test_1%*%res_1$par)
        
        temp_x<-c(test_X_mat[2]*res_2$par[1]+test_D_mat[1])
        temp_base<-bSpline(x=temp_x,knots = knots_2,degree = 3,Boundary.knots = c(lower_knots[2],upper_knots[2]),intercept = TRUE )
        test_decision_2<-c(temp_base%*%res_2$par[-1])
        
        temp_x<-c(test_X_mat[3]*res_3$par[1]+test_D_mat[2])
        temp_base<-bSpline(x=temp_x,knots = knots_3,degree = 3,Boundary.knots = c(lower_knots[3],upper_knots[3]),intercept = TRUE )
        test_decision_3<-c(temp_base%*%res_3$par[-1])
        
        temp_x<-c(test_X_mat[4]*res_4$par[1]+test_D_mat[3])
        temp_base<-bSpline(x=temp_x,knots = knots_4,degree = 3,Boundary.knots = c(lower_knots[4],upper_knots[4]),intercept = TRUE )
        test_decision_4<-c(temp_base%*%res_4$par[-1])
        
        temp_x<-c(test_X_mat[5]*res_5$par[1]+test_D_mat[4])
        temp_base<-bSpline(x=temp_x,knots = knots_5,degree = 3,Boundary.knots = c(lower_knots[5],upper_knots[5]),intercept = TRUE )
        test_decision_5<-c(temp_base%*%res_5$par[-1])
        
        temp_x<-c(test_X_mat[6]*res_6$par[1]+test_D_mat[5])
        temp_base<-bSpline(x=temp_x,knots = knots_6,degree = 3,Boundary.knots = c(lower_knots[6],upper_knots[6]),intercept = TRUE )
        test_decision_6<-c(temp_base%*%res_6$par[-1])
        
        temp_x<-c(test_X_mat[7]*res_7$par[1]+test_D_mat[6])
        temp_base<-bSpline(x=temp_x,knots = knots_7,degree = 3,Boundary.knots = c(lower_knots[7],upper_knots[7]),intercept = TRUE )
        test_decision_7<-c(temp_base%*%res_7$par[-1])
        
        temp_x<-c(test_X_mat[8]*res_8$par[1]+test_D_mat[7])
        temp_base<-bSpline(x=temp_x,knots = knots_8,degree = 3,Boundary.knots = c(lower_knots[8],upper_knots[8]),intercept = TRUE )
        test_decision_8<-c(temp_base%*%res_8$par[-1])
        
        temp_x<-c(test_X_mat[9]*res_9$par[1]+test_D_mat[8])
        temp_base<-bSpline(x=temp_x,knots = knots_9,degree = 3,Boundary.knots = c(lower_knots[9],upper_knots[9]),intercept = TRUE )
        test_decision_9<-c(temp_base%*%res_9$par[-1])
        
        temp_x<-c(test_X_mat[10]*res_10$par[1]+test_D_mat[9])
        temp_base<-bSpline(x=temp_x,knots = knots_10,degree = 3,Boundary.knots = c(lower_knots[10],upper_knots[10]),intercept = TRUE )
        test_decision_10<-c(temp_base%*%res_10$par[-1])
        
        test_decision_our <- c(test_decision_1,test_decision_2,test_decision_3,test_decision_4,test_decision_5,test_decision_6,test_decision_7,test_decision_8,test_decision_9,test_decision_10)
      }
      
      test_cost_our <- 0
      test_cost_oracle <- 0
      temp_level_our <- test_decision_our[1]
      temp_level_oracle <- test_decision_oracle[1]
      for(qqqq in 1:T){
        test_cost_our<-test_cost_our+(h_cost[qqqq]*pmax(0,temp_level_our-test_D_mat[qqqq])+b_cost[qqqq]*pmax(0,test_D_mat[qqqq]-temp_level_our))
        test_cost_oracle<-test_cost_oracle+(h_cost[qqqq]*pmax(0,temp_level_oracle-test_D_mat[qqqq])+b_cost[qqqq]*pmax(0,test_D_mat[qqqq]-temp_level_oracle))
        
        if(qqqq<T){
          temp_level_our<-pmax(temp_level_our-test_D_mat[qqqq],test_decision_our[qqqq+1])
          temp_level_oracle<-pmax(temp_level_oracle-test_D_mat[qqqq],test_decision_oracle[qqqq+1])
        }
      }
      
      our_time <- (proc.time() - start_time)["elapsed"]
      
      final_res <- list(test_cost_our, test_cost_oracle, our_time)
      
      my_final_res[[cut_off + 1]] <- final_res
      print(cut_off)
    }
    
    return(my_final_res)
  }
  res<-try(my_simulation_procedure(iter),silent = TRUE)
  return(res)
}


clnum <- 20#detectCores()/2
cl <- makeCluster(getOption("cl.cores", clnum))
res<-parLapply(cl, 1:500,  simulation)
stopCluster(cl)

save(res,file = "T equals 10 with covariates online.rdata")