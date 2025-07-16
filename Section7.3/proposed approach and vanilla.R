rm(list=ls())
library(parallel)
iter <- 1
simulation<-function(iter){
  my_simulation_procedure<-function(iter){
    library(splines2)
    library(pracma)
    library(Rsolnp)
    # load("nocov_10regret_n200.rdata")
    n<-200
    T<-10
    D_mat<-matrix(0,n,T)
    for(i in 1:n){
      D_mat[i,1]<-runif(1,min = 1,5)
      D_mat[i,2]<-(D_mat[i,1])+runif(1,min = -0.5,0.5)
      for(j in 3:T){
        D_mat[i,j]<-D_mat[i,j-1]+0.2+runif(1,min = -0.2,0.2)
      }
    }
    
    load("optimal_policy_T10_nocov.rdata")
    b_cost<-rep(1,T)
    h_cost<-rep(1,T)
    b_cost[1]<-2
    full_D_mat <- D_mat
    
    my_final_res <- as.list(rep(NA, n))
    
    for(cut_off in 0 : (n - 1)){
      knots_2 <- res[[1]][[1]]
      knots_3 <- res[[1]][[2]]
      knots_4 <- res[[1]][[3]]
      knots_5 <- res[[1]][[4]]
      knots_6 <- res[[1]][[5]]
      knots_7 <- res[[1]][[6]]
      knots_8 <- res[[1]][[7]]
      knots_9 <- res[[1]][[8]]
      knots_10 <- res[[1]][[9]]
      
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
      
      test_decision_1<-pars1
      temp_base<-bSpline(x=test_D_mat[1],knots = knots_2,degree = 3,Boundary.knots = c(1,5),intercept = TRUE )
      test_decision_2<-c(temp_base%*%pars2)
      temp_base<-bSpline(x=test_D_mat[2],knots = knots_3,degree = 3,Boundary.knots = c(0.5,5.5),intercept = TRUE )
      test_decision_3<-c(temp_base%*%pars3)
      temp_base<-bSpline(x=test_D_mat[3],knots = knots_4,degree = 3,Boundary.knots = c(0.5,5.5+0.4*1),intercept = TRUE )
      test_decision_4<-c(temp_base%*%pars4)
      temp_base<-bSpline(x=test_D_mat[4],knots = knots_5,degree = 3,Boundary.knots = c(0.5,5.5+0.4*2),intercept = TRUE )
      test_decision_5<-c(temp_base%*%pars5)
      temp_base<-bSpline(x=test_D_mat[5],knots = knots_6,degree = 3,Boundary.knots = c(0.5,5.5+0.4*3),intercept = TRUE )
      test_decision_6<-c(temp_base%*%pars6)
      temp_base<-bSpline(x=test_D_mat[6],knots = knots_7,degree = 3,Boundary.knots = c(0.5,5.5+0.4*4),intercept = TRUE )
      test_decision_7<-c(temp_base%*%pars7)
      temp_base<-bSpline(x=test_D_mat[7],knots = knots_8,degree = 3,Boundary.knots = c(0.5,5.5+0.4*5),intercept = TRUE )
      test_decision_8<-c(temp_base%*%pars8)
      temp_base<-bSpline(x=test_D_mat[8],knots = knots_9,degree = 3,Boundary.knots = c(0.5,5.5+0.4*6),intercept = TRUE )
      test_decision_9<-c(temp_base%*%pars9)
      temp_base<-bSpline(x=test_D_mat[9],knots = knots_10,degree = 3,Boundary.knots = c(0.5,5.5+0.4*7),intercept = TRUE )
      test_decision_10<-c(temp_base%*%pars10)
      
      test_decision_oracle<- c(test_decision_1,test_decision_2,test_decision_3,test_decision_4,test_decision_5,test_decision_6,test_decision_7,test_decision_8,test_decision_9,test_decision_10)
      
      ### proposed approach
      if(cut_off == 0){
        start_time <- proc.time()
        test_decision_our <- rep(2, T)
        our_time <- (proc.time() - start_time)["elapsed"]
        our_CV_time <- 0
      }
      
      if(cut_off == 1){
        start_time <- proc.time()
        test_decision_our <- full_D_mat[1, ]
        our_time <- (proc.time() - start_time)["elapsed"]
        our_CV_time <- 0
      }
      
      if(cut_off > 1){
        m_n_candidate <- 1:10
        D_mat <- full_D_mat[1:cut_off, ]
        if(cut_off < 5){
          m_n <- sample(m_n_candidate, 1)
          our_CV_time <- 0
        }else{
          ########################################################################################################################################################################################################################################
          ########################################################################################################################################################################################################################################
          start_time <- proc.time()
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
            test_D_mat<-D_mat[test_index, ]
            
            train_index <- c(my_split[-fold,])
            train_index <- train_index[which(train_index < 5e3)]
            train_D_mat<-D_mat[train_index, ]
            
            for(ppp in 1:length(m_n_candidate)){
              m_n<-m_n_candidate[ppp]
              cut_quantile<-seq(from=0,to=1,length.out=m_n+2)[-c(1,m_n+2)]
              
              knots_10<-quantile(x=train_D_mat[,9],probs=cut_quantile)
              basis_mat_10<-bSpline(x=train_D_mat[,9],knots = knots_10,degree = 3,Boundary.knots = c(0.5,5.5+0.4*7),intercept = TRUE )
              
              price_10<-function(alpha){
                temp_level<-c(basis_mat_10%*%alpha)
                return(sum(h_cost[10]*pmax(0,temp_level-train_D_mat[,10])+b_cost[10]*pmax(0,train_D_mat[,10]-temp_level)))
              }
              
              price_10_empirical<-function(temp_level){
                return(sum(h_cost[10]*pmax(0,temp_level-train_D_mat[,10])+b_cost[10]*pmax(0,train_D_mat[,10]-temp_level)))
              }
              
              price_10_empirical_var_est<-function(temp_level){
                return(h_cost[10]*pmax(0,temp_level-train_D_mat[,10])+b_cost[10]*pmax(0,train_D_mat[,10]-temp_level))
              }
              
              res_10<-solnp(pars = rep(1,m_n+4),fun = price_10,LB=rep(0,m_n+4),control = list(trace=0))
              best_decision_10<-c(basis_mat_10%*%res_10$pars)
              
              
              
              
              knots_9<-quantile(x=train_D_mat[,8],probs=cut_quantile)
              basis_mat_9<-bSpline(x=train_D_mat[,8],knots = knots_9,degree = 3,Boundary.knots = c(0.5,5.5+0.4*6),intercept = TRUE )
              
              price_9<-function(alpha){
                temp_level<-c(basis_mat_9%*%alpha)
                return(sum(h_cost[9]*pmax(0,temp_level-train_D_mat[,9])+b_cost[9]*pmax(0,train_D_mat[,9]-temp_level))+price_10_empirical(pmax(best_decision_10,temp_level-train_D_mat[,9])) )
              }
              
              price_9_empirical<-function(temp_level){
                return(sum(h_cost[9]*pmax(0,temp_level-train_D_mat[,9])+b_cost[9]*pmax(0,train_D_mat[,9]-temp_level))+price_10_empirical(pmax(best_decision_10,temp_level-train_D_mat[,9])) )
              }
              
              price_9_empirical_var_est<-function(temp_level){
                return(h_cost[9]*pmax(0,temp_level-train_D_mat[,9])+b_cost[9]*pmax(0,train_D_mat[,9]-temp_level)+price_10_empirical_var_est(pmax(best_decision_10,temp_level-train_D_mat[,9])))
              }
              
              res_9<-solnp(pars = rep(1,m_n+4),fun = price_9,LB=rep(0,m_n+4),control = list(trace=0))
              best_decision_9<-c(basis_mat_9%*%res_9$pars)
              
              
              
              knots_8<-quantile(x=train_D_mat[,7],probs=cut_quantile)
              basis_mat_8<-bSpline(x=train_D_mat[,7],knots = knots_8,degree = 3,Boundary.knots = c(0.5,5.5+0.4*5),intercept = TRUE )
              
              price_8<-function(alpha){
                temp_level<-c(basis_mat_8%*%alpha)
                return(sum(h_cost[8]*pmax(0,temp_level-train_D_mat[,8])+b_cost[8]*pmax(0,train_D_mat[,8]-temp_level))+price_9_empirical(pmax(best_decision_9,temp_level-train_D_mat[,8])) )
              }
              
              price_8_empirical<-function(temp_level){
                return(sum(h_cost[8]*pmax(0,temp_level-train_D_mat[,8])+b_cost[8]*pmax(0,train_D_mat[,8]-temp_level))+price_9_empirical(pmax(best_decision_9,temp_level-train_D_mat[,8])) )
              }
              
              price_8_empirical_var_est<-function(temp_level){
                return(h_cost[8]*pmax(0,temp_level-train_D_mat[,8])+b_cost[8]*pmax(0,train_D_mat[,8]-temp_level)+price_9_empirical_var_est(pmax(best_decision_9,temp_level-train_D_mat[,8])))
              }
              
              res_8<-solnp(pars = rep(1,m_n+4),fun = price_8,LB=rep(0,m_n+4),control = list(trace=0))
              best_decision_8<-c(basis_mat_8%*%res_8$pars)
              
              
              
              knots_7<-quantile(x=train_D_mat[,6],probs=cut_quantile)
              basis_mat_7<-bSpline(x=train_D_mat[,6],knots = knots_7,degree = 3,Boundary.knots = c(0.5,5.5+0.4*4),intercept = TRUE )
              
              price_7<-function(alpha){
                temp_level<-c(basis_mat_7%*%alpha)
                return(sum(h_cost[7]*pmax(0,temp_level-train_D_mat[,7])+b_cost[7]*pmax(0,train_D_mat[,7]-temp_level))+price_8_empirical(pmax(best_decision_8,temp_level-train_D_mat[,7])) )
              }
              
              price_7_empirical<-function(temp_level){
                return(sum(h_cost[7]*pmax(0,temp_level-train_D_mat[,7])+b_cost[7]*pmax(0,train_D_mat[,7]-temp_level))+price_8_empirical(pmax(best_decision_8,temp_level-train_D_mat[,7])) )
              }
              
              price_7_empirical_var_est<-function(temp_level){
                return(h_cost[7]*pmax(0,temp_level-train_D_mat[,7])+b_cost[7]*pmax(0,train_D_mat[,7]-temp_level)+price_8_empirical_var_est(pmax(best_decision_8,temp_level-train_D_mat[,7])))
              }
              
              res_7<-solnp(pars = rep(1,m_n+4),fun = price_7,LB=rep(0,m_n+4),control = list(trace=0))
              best_decision_7<-c(basis_mat_7%*%res_7$pars)
              
              
              knots_6<-quantile(x=train_D_mat[,5],probs=cut_quantile)
              basis_mat_6<-bSpline(x=train_D_mat[,5],knots = knots_6,degree = 3,Boundary.knots = c(0.5,5.5+0.4*3),intercept = TRUE )
              
              price_6<-function(alpha){
                temp_level<-c(basis_mat_6%*%alpha)
                return(sum(h_cost[6]*pmax(0,temp_level-train_D_mat[,6])+b_cost[6]*pmax(0,train_D_mat[,6]-temp_level))+price_7_empirical(pmax(best_decision_7,temp_level-train_D_mat[,6])) )
              }
              
              price_6_empirical<-function(temp_level){
                return(sum(h_cost[6]*pmax(0,temp_level-train_D_mat[,6])+b_cost[6]*pmax(0,train_D_mat[,6]-temp_level))+price_7_empirical(pmax(best_decision_7,temp_level-train_D_mat[,6])) )
              }
              
              price_6_empirical_var_est<-function(temp_level){
                return(h_cost[6]*pmax(0,temp_level-train_D_mat[,6])+b_cost[6]*pmax(0,train_D_mat[,6]-temp_level)+price_7_empirical_var_est(pmax(best_decision_7,temp_level-train_D_mat[,6])))
              }
              
              res_6<-solnp(pars = rep(1,m_n+4),fun = price_6,LB=rep(0,m_n+4),control = list(trace=0))
              best_decision_6<-c(basis_mat_6%*%res_6$pars)
              
              
              
              knots_5<-quantile(x=train_D_mat[,4],probs=cut_quantile)
              basis_mat_5<-bSpline(x=train_D_mat[,4],knots = knots_5,degree = 3,Boundary.knots = c(0.5,5.5+0.4*2),intercept = TRUE )
              
              price_5<-function(alpha){
                temp_level<-c(basis_mat_5%*%alpha)
                return(sum(h_cost[5]*pmax(0,temp_level-train_D_mat[,5])+b_cost[5]*pmax(0,train_D_mat[,5]-temp_level))+price_6_empirical(pmax(best_decision_6,temp_level-train_D_mat[,5])) )
              }
              
              price_5_empirical<-function(temp_level){
                return(sum(h_cost[5]*pmax(0,temp_level-train_D_mat[,5])+b_cost[5]*pmax(0,train_D_mat[,5]-temp_level))+price_6_empirical(pmax(best_decision_6,temp_level-train_D_mat[,5])) )
              }
              
              price_5_empirical_var_est<-function(temp_level){
                return(h_cost[5]*pmax(0,temp_level-train_D_mat[,5])+b_cost[5]*pmax(0,train_D_mat[,5]-temp_level)+price_6_empirical_var_est(pmax(best_decision_6,temp_level-train_D_mat[,5])))
              }
              
              res_5<-solnp(pars = rep(1,m_n+4),fun = price_5,LB=rep(0,m_n+4),control = list(trace=0))
              best_decision_5<-c(basis_mat_5%*%res_5$pars)
              
              
              knots_4<-quantile(x=train_D_mat[,3],probs=cut_quantile)
              basis_mat_4<-bSpline(x=train_D_mat[,3],knots = knots_4,degree = 3,Boundary.knots = c(0.5,5.5+0.4*1),intercept = TRUE )
              
              price_4<-function(alpha){
                temp_level<-c(basis_mat_4%*%alpha)
                return(sum(h_cost[4]*pmax(0,temp_level-train_D_mat[,4])+b_cost[4]*pmax(0,train_D_mat[,4]-temp_level))+price_5_empirical(pmax(best_decision_5,temp_level-train_D_mat[,4])) )
              }
              
              price_4_empirical<-function(temp_level){
                return(sum(h_cost[4]*pmax(0,temp_level-train_D_mat[,4])+b_cost[4]*pmax(0,train_D_mat[,4]-temp_level))+price_5_empirical(pmax(best_decision_5,temp_level-train_D_mat[,4])) )
              }
              
              price_4_empirical_var_est<-function(temp_level){
                return(h_cost[4]*pmax(0,temp_level-train_D_mat[,4])+b_cost[4]*pmax(0,train_D_mat[,4]-temp_level)+price_5_empirical_var_est(pmax(best_decision_5,temp_level-train_D_mat[,4])))
              }
              
              res_4<-solnp(pars = rep(1,m_n+4),fun = price_4,LB=rep(0,m_n+4),control = list(trace=0))
              best_decision_4<-c(basis_mat_4%*%res_4$pars)
              
              
              
              
              
              knots_3<-quantile(x=train_D_mat[,2],probs=cut_quantile)
              basis_mat_3<-bSpline(x=train_D_mat[,2],knots = knots_3,degree = 3,Boundary.knots = c(0.5,5.5),intercept = TRUE )
              
              price_3<-function(alpha){
                temp_level<-c(basis_mat_3%*%alpha)
                return(sum(h_cost[3]*pmax(0,temp_level-train_D_mat[,3])+b_cost[3]*pmax(0,train_D_mat[,3]-temp_level))+price_4_empirical(pmax(best_decision_4,temp_level-train_D_mat[,3])) )
              }
              
              price_3_empirical<-function(temp_level){
                return(sum(h_cost[3]*pmax(0,temp_level-train_D_mat[,3])+b_cost[3]*pmax(0,train_D_mat[,3]-temp_level))+price_4_empirical(pmax(best_decision_4,temp_level-train_D_mat[,3])) )
              }
              
              price_3_empirical_var_est<-function(temp_level){
                return(h_cost[3]*pmax(0,temp_level-train_D_mat[,3])+b_cost[3]*pmax(0,train_D_mat[,3]-temp_level)+price_4_empirical_var_est(pmax(best_decision_4,temp_level-train_D_mat[,3])))
              }
              
              res_3<-solnp(pars = rep(1,m_n+4),fun = price_3,LB=rep(0,m_n+4),control = list(trace=0))
              best_decision_3<-c(basis_mat_3%*%res_3$pars)
              
              
              
              
              knots_2<-quantile(x=train_D_mat[,1],probs=cut_quantile)
              basis_mat_2<-bSpline(x=train_D_mat[,1],knots = knots_2,degree = 3,Boundary.knots = c(1,5),intercept = TRUE )
              
              price_2<-function(alpha){
                temp_level<-c(basis_mat_2%*%alpha)
                return(sum(h_cost[2]*pmax(0,temp_level-train_D_mat[,2])+b_cost[2]*pmax(0,train_D_mat[,2]-temp_level))+price_3_empirical(pmax(temp_level-train_D_mat[,2],best_decision_3)))
              }
              
              price_2_empirical<-function(temp_level){
                return(sum(h_cost[2]*pmax(0,temp_level-train_D_mat[,2])+b_cost[2]*pmax(0,train_D_mat[,2]-temp_level))+price_3_empirical(pmax(temp_level-train_D_mat[,2],best_decision_3)))
              }
              
              price_2_empirical_var_est<-function(temp_level){
                return(h_cost[2]*pmax(0,temp_level-train_D_mat[,2])+b_cost[2]*pmax(0,train_D_mat[,2]-temp_level)+price_3_empirical_var_est(pmax(temp_level-train_D_mat[,2],best_decision_3)))
              }
              
              res_2<-solnp(pars = rep(1,m_n+4),fun = price_2,LB=rep(0,m_n+4),control = list(trace=0))
              best_decision_2<-c(basis_mat_2%*%res_2$pars)
              
              temp_x<-seq(1,5,0.01)
              temp_base<-bSpline(x=temp_x,knots = knots_2,degree = 3,Boundary.knots = c(1,5),intercept = TRUE )
              S2_est<-c(temp_base%*%res_2$par)
              # plot(temp_x,S2_est,type = "l",lwd=3)
              # points(temp_x,temp_x,lwd=3,lty=2,type = "l",col="red")
              
              price_1<-function(temp_level){
                return(sum(h_cost[1]*pmax(0,temp_level-train_D_mat[,1])+b_cost[1]*pmax(0,train_D_mat[,1]-temp_level))+price_2_empirical(pmax(temp_level-train_D_mat[,1],best_decision_2)))
              }
              price_1_empirical_var_est<-function(temp_level){
                return((h_cost[1]*pmax(0,temp_level-train_D_mat[,1])+b_cost[1]*pmax(0,train_D_mat[,1]-temp_level))+price_2_empirical_var_est(pmax(temp_level-train_D_mat[,1],best_decision_2)))
              }
              res_1<-optimize(f=price_1,interval = c(0,50))
              
              if(is.null(nrow(test_D_mat))){
                test_decision_1<-res_1$minimum
                temp_base<-bSpline(x=test_D_mat[1],knots = knots_2,degree = 3,Boundary.knots = c(1,5),intercept = TRUE )
                test_decision_2<-c(temp_base%*%res_2$pars)
                temp_base<-bSpline(x=test_D_mat[2],knots = knots_3,degree = 3,Boundary.knots = c(0.5,5.5),intercept = TRUE )
                test_decision_3<-c(temp_base%*%res_3$pars)
                temp_base<-bSpline(x=test_D_mat[3],knots = knots_4,degree = 3,Boundary.knots = c(0.5,5.5+0.4*1),intercept = TRUE )
                test_decision_4<-c(temp_base%*%res_4$pars)
                temp_base<-bSpline(x=test_D_mat[4],knots = knots_5,degree = 3,Boundary.knots = c(0.5,5.5+0.4*2),intercept = TRUE )
                test_decision_5<-c(temp_base%*%res_5$pars)
                temp_base<-bSpline(x=test_D_mat[5],knots = knots_6,degree = 3,Boundary.knots = c(0.5,5.5+0.4*3),intercept = TRUE )
                test_decision_6<-c(temp_base%*%res_6$pars)
                temp_base<-bSpline(x=test_D_mat[6],knots = knots_7,degree = 3,Boundary.knots = c(0.5,5.5+0.4*4),intercept = TRUE )
                test_decision_7<-c(temp_base%*%res_7$pars)
                temp_base<-bSpline(x=test_D_mat[7],knots = knots_8,degree = 3,Boundary.knots = c(0.5,5.5+0.4*5),intercept = TRUE )
                test_decision_8<-c(temp_base%*%res_8$pars)
                temp_base<-bSpline(x=test_D_mat[8],knots = knots_9,degree = 3,Boundary.knots = c(0.5,5.5+0.4*6),intercept = TRUE )
                test_decision_9<-c(temp_base%*%res_9$pars)
                temp_base<-bSpline(x=test_D_mat[9],knots = knots_10,degree = 3,Boundary.knots = c(0.5,5.5+0.4*7),intercept = TRUE )
                test_decision_10<-c(temp_base%*%res_10$pars)
                
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
                test_decision_1<-res_1$minimum
                temp_base<-bSpline(x=test_D_mat[,1],knots = knots_2,degree = 3,Boundary.knots = c(1,5),intercept = TRUE )
                test_decision_2<-c(temp_base%*%res_2$pars)
                temp_base<-bSpline(x=test_D_mat[,2],knots = knots_3,degree = 3,Boundary.knots = c(0.5,5.5),intercept = TRUE )
                test_decision_3<-c(temp_base%*%res_3$pars)
                temp_base<-bSpline(x=test_D_mat[,3],knots = knots_4,degree = 3,Boundary.knots = c(0.5,5.5+0.4*1),intercept = TRUE )
                test_decision_4<-c(temp_base%*%res_4$pars)
                temp_base<-bSpline(x=test_D_mat[,4],knots = knots_5,degree = 3,Boundary.knots = c(0.5,5.5+0.4*2),intercept = TRUE )
                test_decision_5<-c(temp_base%*%res_5$pars)
                temp_base<-bSpline(x=test_D_mat[,5],knots = knots_6,degree = 3,Boundary.knots = c(0.5,5.5+0.4*3),intercept = TRUE )
                test_decision_6<-c(temp_base%*%res_6$pars)
                temp_base<-bSpline(x=test_D_mat[,6],knots = knots_7,degree = 3,Boundary.knots = c(0.5,5.5+0.4*4),intercept = TRUE )
                test_decision_7<-c(temp_base%*%res_7$pars)
                temp_base<-bSpline(x=test_D_mat[,7],knots = knots_8,degree = 3,Boundary.knots = c(0.5,5.5+0.4*5),intercept = TRUE )
                test_decision_8<-c(temp_base%*%res_8$pars)
                temp_base<-bSpline(x=test_D_mat[,8],knots = knots_9,degree = 3,Boundary.knots = c(0.5,5.5+0.4*6),intercept = TRUE )
                test_decision_9<-c(temp_base%*%res_9$pars)
                temp_base<-bSpline(x=test_D_mat[,9],knots = knots_10,degree = 3,Boundary.knots = c(0.5,5.5+0.4*7),intercept = TRUE )
                test_decision_10<-c(temp_base%*%res_10$pars)
                
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
          our_CV_time <- (proc.time() - start_time)["elapsed"]
        }
        
        start_time <- proc.time()
        cut_quantile<-seq(from=0,to=1,length.out=m_n+2)[-c(1,m_n+2)]
        
        train_D_mat <- full_D_mat[1:cut_off, ]
        
        
        knots_10<-quantile(x=train_D_mat[,9],probs=cut_quantile)
        basis_mat_10<-bSpline(x=train_D_mat[,9],knots = knots_10,degree = 3,Boundary.knots = c(0.5,5.5+0.4*7),intercept = TRUE )
        
        price_10<-function(alpha){
          temp_level<-c(basis_mat_10%*%alpha)
          return(sum(h_cost[10]*pmax(0,temp_level-train_D_mat[,10])+b_cost[10]*pmax(0,train_D_mat[,10]-temp_level)))
        }
        
        price_10_empirical<-function(temp_level){
          return(sum(h_cost[10]*pmax(0,temp_level-train_D_mat[,10])+b_cost[10]*pmax(0,train_D_mat[,10]-temp_level)))
        }
        
        price_10_empirical_var_est<-function(temp_level){
          return(h_cost[10]*pmax(0,temp_level-train_D_mat[,10])+b_cost[10]*pmax(0,train_D_mat[,10]-temp_level))
        }
        
        res_10<-solnp(pars = rep(1,m_n+4),fun = price_10,LB=rep(0,m_n+4),control = list(trace=0))
        best_decision_10<-c(basis_mat_10%*%res_10$pars)
        
        
        
        
        knots_9<-quantile(x=train_D_mat[,8],probs=cut_quantile)
        basis_mat_9<-bSpline(x=train_D_mat[,8],knots = knots_9,degree = 3,Boundary.knots = c(0.5,5.5+0.4*6),intercept = TRUE )
        
        price_9<-function(alpha){
          temp_level<-c(basis_mat_9%*%alpha)
          return(sum(h_cost[9]*pmax(0,temp_level-train_D_mat[,9])+b_cost[9]*pmax(0,train_D_mat[,9]-temp_level))+price_10_empirical(pmax(best_decision_10,temp_level-train_D_mat[,9])) )
        }
        
        price_9_empirical<-function(temp_level){
          return(sum(h_cost[9]*pmax(0,temp_level-train_D_mat[,9])+b_cost[9]*pmax(0,train_D_mat[,9]-temp_level))+price_10_empirical(pmax(best_decision_10,temp_level-train_D_mat[,9])) )
        }
        
        price_9_empirical_var_est<-function(temp_level){
          return(h_cost[9]*pmax(0,temp_level-train_D_mat[,9])+b_cost[9]*pmax(0,train_D_mat[,9]-temp_level)+price_10_empirical_var_est(pmax(best_decision_10,temp_level-train_D_mat[,9])))
        }
        
        res_9<-solnp(pars = rep(1,m_n+4),fun = price_9,LB=rep(0,m_n+4),control = list(trace=0))
        best_decision_9<-c(basis_mat_9%*%res_9$pars)
        
        
        
        knots_8<-quantile(x=train_D_mat[,7],probs=cut_quantile)
        basis_mat_8<-bSpline(x=train_D_mat[,7],knots = knots_8,degree = 3,Boundary.knots = c(0.5,5.5+0.4*5),intercept = TRUE )
        
        price_8<-function(alpha){
          temp_level<-c(basis_mat_8%*%alpha)
          return(sum(h_cost[8]*pmax(0,temp_level-train_D_mat[,8])+b_cost[8]*pmax(0,train_D_mat[,8]-temp_level))+price_9_empirical(pmax(best_decision_9,temp_level-train_D_mat[,8])) )
        }
        
        price_8_empirical<-function(temp_level){
          return(sum(h_cost[8]*pmax(0,temp_level-train_D_mat[,8])+b_cost[8]*pmax(0,train_D_mat[,8]-temp_level))+price_9_empirical(pmax(best_decision_9,temp_level-train_D_mat[,8])) )
        }
        
        price_8_empirical_var_est<-function(temp_level){
          return(h_cost[8]*pmax(0,temp_level-train_D_mat[,8])+b_cost[8]*pmax(0,train_D_mat[,8]-temp_level)+price_9_empirical_var_est(pmax(best_decision_9,temp_level-train_D_mat[,8])))
        }
        
        res_8<-solnp(pars = rep(1,m_n+4),fun = price_8,LB=rep(0,m_n+4),control = list(trace=0))
        best_decision_8<-c(basis_mat_8%*%res_8$pars)
        
        
        
        knots_7<-quantile(x=train_D_mat[,6],probs=cut_quantile)
        basis_mat_7<-bSpline(x=train_D_mat[,6],knots = knots_7,degree = 3,Boundary.knots = c(0.5,5.5+0.4*4),intercept = TRUE )
        
        price_7<-function(alpha){
          temp_level<-c(basis_mat_7%*%alpha)
          return(sum(h_cost[7]*pmax(0,temp_level-train_D_mat[,7])+b_cost[7]*pmax(0,train_D_mat[,7]-temp_level))+price_8_empirical(pmax(best_decision_8,temp_level-train_D_mat[,7])) )
        }
        
        price_7_empirical<-function(temp_level){
          return(sum(h_cost[7]*pmax(0,temp_level-train_D_mat[,7])+b_cost[7]*pmax(0,train_D_mat[,7]-temp_level))+price_8_empirical(pmax(best_decision_8,temp_level-train_D_mat[,7])) )
        }
        
        price_7_empirical_var_est<-function(temp_level){
          return(h_cost[7]*pmax(0,temp_level-train_D_mat[,7])+b_cost[7]*pmax(0,train_D_mat[,7]-temp_level)+price_8_empirical_var_est(pmax(best_decision_8,temp_level-train_D_mat[,7])))
        }
        
        res_7<-solnp(pars = rep(1,m_n+4),fun = price_7,LB=rep(0,m_n+4),control = list(trace=0))
        best_decision_7<-c(basis_mat_7%*%res_7$pars)
        
        
        knots_6<-quantile(x=train_D_mat[,5],probs=cut_quantile)
        basis_mat_6<-bSpline(x=train_D_mat[,5],knots = knots_6,degree = 3,Boundary.knots = c(0.5,5.5+0.4*3),intercept = TRUE )
        
        price_6<-function(alpha){
          temp_level<-c(basis_mat_6%*%alpha)
          return(sum(h_cost[6]*pmax(0,temp_level-train_D_mat[,6])+b_cost[6]*pmax(0,train_D_mat[,6]-temp_level))+price_7_empirical(pmax(best_decision_7,temp_level-train_D_mat[,6])) )
        }
        
        price_6_empirical<-function(temp_level){
          return(sum(h_cost[6]*pmax(0,temp_level-train_D_mat[,6])+b_cost[6]*pmax(0,train_D_mat[,6]-temp_level))+price_7_empirical(pmax(best_decision_7,temp_level-train_D_mat[,6])) )
        }
        
        price_6_empirical_var_est<-function(temp_level){
          return(h_cost[6]*pmax(0,temp_level-train_D_mat[,6])+b_cost[6]*pmax(0,train_D_mat[,6]-temp_level)+price_7_empirical_var_est(pmax(best_decision_7,temp_level-train_D_mat[,6])))
        }
        
        res_6<-solnp(pars = rep(1,m_n+4),fun = price_6,LB=rep(0,m_n+4),control = list(trace=0))
        best_decision_6<-c(basis_mat_6%*%res_6$pars)
        
        
        
        knots_5<-quantile(x=train_D_mat[,4],probs=cut_quantile)
        basis_mat_5<-bSpline(x=train_D_mat[,4],knots = knots_5,degree = 3,Boundary.knots = c(0.5,5.5+0.4*2),intercept = TRUE )
        
        price_5<-function(alpha){
          temp_level<-c(basis_mat_5%*%alpha)
          return(sum(h_cost[5]*pmax(0,temp_level-train_D_mat[,5])+b_cost[5]*pmax(0,train_D_mat[,5]-temp_level))+price_6_empirical(pmax(best_decision_6,temp_level-train_D_mat[,5])) )
        }
        
        price_5_empirical<-function(temp_level){
          return(sum(h_cost[5]*pmax(0,temp_level-train_D_mat[,5])+b_cost[5]*pmax(0,train_D_mat[,5]-temp_level))+price_6_empirical(pmax(best_decision_6,temp_level-train_D_mat[,5])) )
        }
        
        price_5_empirical_var_est<-function(temp_level){
          return(h_cost[5]*pmax(0,temp_level-train_D_mat[,5])+b_cost[5]*pmax(0,train_D_mat[,5]-temp_level)+price_6_empirical_var_est(pmax(best_decision_6,temp_level-train_D_mat[,5])))
        }
        
        res_5<-solnp(pars = rep(1,m_n+4),fun = price_5,LB=rep(0,m_n+4),control = list(trace=0))
        best_decision_5<-c(basis_mat_5%*%res_5$pars)
        
        
        knots_4<-quantile(x=train_D_mat[,3],probs=cut_quantile)
        basis_mat_4<-bSpline(x=train_D_mat[,3],knots = knots_4,degree = 3,Boundary.knots = c(0.5,5.5+0.4*1),intercept = TRUE )
        
        price_4<-function(alpha){
          temp_level<-c(basis_mat_4%*%alpha)
          return(sum(h_cost[4]*pmax(0,temp_level-train_D_mat[,4])+b_cost[4]*pmax(0,train_D_mat[,4]-temp_level))+price_5_empirical(pmax(best_decision_5,temp_level-train_D_mat[,4])) )
        }
        
        price_4_empirical<-function(temp_level){
          return(sum(h_cost[4]*pmax(0,temp_level-train_D_mat[,4])+b_cost[4]*pmax(0,train_D_mat[,4]-temp_level))+price_5_empirical(pmax(best_decision_5,temp_level-train_D_mat[,4])) )
        }
        
        price_4_empirical_var_est<-function(temp_level){
          return(h_cost[4]*pmax(0,temp_level-train_D_mat[,4])+b_cost[4]*pmax(0,train_D_mat[,4]-temp_level)+price_5_empirical_var_est(pmax(best_decision_5,temp_level-train_D_mat[,4])))
        }
        
        res_4<-solnp(pars = rep(1,m_n+4),fun = price_4,LB=rep(0,m_n+4),control = list(trace=0))
        best_decision_4<-c(basis_mat_4%*%res_4$pars)
        
        
        
        
        
        knots_3<-quantile(x=train_D_mat[,2],probs=cut_quantile)
        basis_mat_3<-bSpline(x=train_D_mat[,2],knots = knots_3,degree = 3,Boundary.knots = c(0.5,5.5),intercept = TRUE )
        
        price_3<-function(alpha){
          temp_level<-c(basis_mat_3%*%alpha)
          return(sum(h_cost[3]*pmax(0,temp_level-train_D_mat[,3])+b_cost[3]*pmax(0,train_D_mat[,3]-temp_level))+price_4_empirical(pmax(best_decision_4,temp_level-train_D_mat[,3])) )
        }
        
        price_3_empirical<-function(temp_level){
          return(sum(h_cost[3]*pmax(0,temp_level-train_D_mat[,3])+b_cost[3]*pmax(0,train_D_mat[,3]-temp_level))+price_4_empirical(pmax(best_decision_4,temp_level-train_D_mat[,3])) )
        }
        
        price_3_empirical_var_est<-function(temp_level){
          return(h_cost[3]*pmax(0,temp_level-train_D_mat[,3])+b_cost[3]*pmax(0,train_D_mat[,3]-temp_level)+price_4_empirical_var_est(pmax(best_decision_4,temp_level-train_D_mat[,3])))
        }
        
        res_3<-solnp(pars = rep(1,m_n+4),fun = price_3,LB=rep(0,m_n+4),control = list(trace=0))
        best_decision_3<-c(basis_mat_3%*%res_3$pars)
        
        
        
        
        knots_2<-quantile(x=train_D_mat[,1],probs=cut_quantile)
        basis_mat_2<-bSpline(x=train_D_mat[,1],knots = knots_2,degree = 3,Boundary.knots = c(1,5),intercept = TRUE )
        
        price_2<-function(alpha){
          temp_level<-c(basis_mat_2%*%alpha)
          return(sum(h_cost[2]*pmax(0,temp_level-train_D_mat[,2])+b_cost[2]*pmax(0,train_D_mat[,2]-temp_level))+price_3_empirical(pmax(temp_level-train_D_mat[,2],best_decision_3)))
        }
        
        price_2_empirical<-function(temp_level){
          return(sum(h_cost[2]*pmax(0,temp_level-train_D_mat[,2])+b_cost[2]*pmax(0,train_D_mat[,2]-temp_level))+price_3_empirical(pmax(temp_level-train_D_mat[,2],best_decision_3)))
        }
        
        price_2_empirical_var_est<-function(temp_level){
          return(h_cost[2]*pmax(0,temp_level-train_D_mat[,2])+b_cost[2]*pmax(0,train_D_mat[,2]-temp_level)+price_3_empirical_var_est(pmax(temp_level-train_D_mat[,2],best_decision_3)))
        }
        
        res_2<-solnp(pars = rep(1,m_n+4),fun = price_2,LB=rep(0,m_n+4),control = list(trace=0))
        best_decision_2<-c(basis_mat_2%*%res_2$pars)
        
        temp_x<-seq(1,5,0.01)
        temp_base<-bSpline(x=temp_x,knots = knots_2,degree = 3,Boundary.knots = c(1,5),intercept = TRUE )
        S2_est<-c(temp_base%*%res_2$par)
        # plot(temp_x,S2_est,type = "l",lwd=3)
        # points(temp_x,temp_x,lwd=3,lty=2,type = "l",col="red")
        
        price_1<-function(temp_level){
          return(sum(h_cost[1]*pmax(0,temp_level-train_D_mat[,1])+b_cost[1]*pmax(0,train_D_mat[,1]-temp_level))+price_2_empirical(pmax(temp_level-train_D_mat[,1],best_decision_2)))
        }
        price_1_empirical_var_est<-function(temp_level){
          return((h_cost[1]*pmax(0,temp_level-train_D_mat[,1])+b_cost[1]*pmax(0,train_D_mat[,1]-temp_level))+price_2_empirical_var_est(pmax(temp_level-train_D_mat[,1],best_decision_2)))
        }
        res_1<-optimize(f=price_1,interval = c(0,50))
        
        
        test_D_mat <- full_D_mat[cut_off + 1, ]
        test_decision_1<-res_1$minimum
        temp_base<-bSpline(x=test_D_mat[1],knots = knots_2,degree = 3,Boundary.knots = c(1,5),intercept = TRUE )
        test_decision_2<-c(temp_base%*%res_2$pars)
        temp_base<-bSpline(x=test_D_mat[2],knots = knots_3,degree = 3,Boundary.knots = c(0.5,5.5),intercept = TRUE )
        test_decision_3<-c(temp_base%*%res_3$pars)
        temp_base<-bSpline(x=test_D_mat[3],knots = knots_4,degree = 3,Boundary.knots = c(0.5,5.5+0.4*1),intercept = TRUE )
        test_decision_4<-c(temp_base%*%res_4$pars)
        temp_base<-bSpline(x=test_D_mat[4],knots = knots_5,degree = 3,Boundary.knots = c(0.5,5.5+0.4*2),intercept = TRUE )
        test_decision_5<-c(temp_base%*%res_5$pars)
        temp_base<-bSpline(x=test_D_mat[5],knots = knots_6,degree = 3,Boundary.knots = c(0.5,5.5+0.4*3),intercept = TRUE )
        test_decision_6<-c(temp_base%*%res_6$pars)
        temp_base<-bSpline(x=test_D_mat[6],knots = knots_7,degree = 3,Boundary.knots = c(0.5,5.5+0.4*4),intercept = TRUE )
        test_decision_7<-c(temp_base%*%res_7$pars)
        temp_base<-bSpline(x=test_D_mat[7],knots = knots_8,degree = 3,Boundary.knots = c(0.5,5.5+0.4*5),intercept = TRUE )
        test_decision_8<-c(temp_base%*%res_8$pars)
        temp_base<-bSpline(x=test_D_mat[8],knots = knots_9,degree = 3,Boundary.knots = c(0.5,5.5+0.4*6),intercept = TRUE )
        test_decision_9<-c(temp_base%*%res_9$pars)
        temp_base<-bSpline(x=test_D_mat[9],knots = knots_10,degree = 3,Boundary.knots = c(0.5,5.5+0.4*7),intercept = TRUE )
        test_decision_10<-c(temp_base%*%res_10$pars)
        
        test_decision_our <- c(test_decision_1,test_decision_2,test_decision_3,test_decision_4,test_decision_5,test_decision_6,test_decision_7,test_decision_8,test_decision_9,test_decision_10)
        our_time <- (proc.time() - start_time)["elapsed"]
      }
      
      
      ### vanilla approach
      start_time <- proc.time()
      if(cut_off == 0){
        test_decision_naive <- rep(2, T)
      }
      
      if(cut_off == 1){
        test_decision_naive <- full_D_mat[1, ]
      }
      
      if(cut_off > 1){
        D_mat <- full_D_mat[1:cut_off, ]
        train_D_mat <- full_D_mat[1:cut_off, ]
        ########################################################################################################################################################################################################################################
        ########################################################################################################################################################################################################################################
        ###naive policy
        
        price_10<-function(alpha){
          temp_level<-alpha
          return(sum(h_cost[10]*pmax(0,temp_level-train_D_mat[,10])+b_cost[10]*pmax(0,train_D_mat[,10]-temp_level)))
        }
        
        price_10_empirical<-function(temp_level){
          return(sum(h_cost[10]*pmax(0,temp_level-train_D_mat[,10])+b_cost[10]*pmax(0,train_D_mat[,10]-temp_level)))
        }
        
        price_10_empirical_var_est<-function(temp_level){
          return(h_cost[10]*pmax(0,temp_level-train_D_mat[,10])+b_cost[10]*pmax(0,train_D_mat[,10]-temp_level))
        }
        
        res_10<-solnp(pars = mean(train_D_mat[,10]),fun = price_10,LB = 0,control = list(trace=0))
        best_decision_10<-res_10$pars
        
        
        
        
        price_9<-function(alpha){
          temp_level<-alpha
          return(sum(h_cost[9]*pmax(0,temp_level-train_D_mat[,9])+b_cost[9]*pmax(0,train_D_mat[,9]-temp_level))+price_10_empirical(pmax(best_decision_10,temp_level-train_D_mat[,9])) )
        }
        
        price_9_empirical<-function(temp_level){
          return(sum(h_cost[9]*pmax(0,temp_level-train_D_mat[,9])+b_cost[9]*pmax(0,train_D_mat[,9]-temp_level))+price_10_empirical(pmax(best_decision_10,temp_level-train_D_mat[,9])) )
        }
        
        price_9_empirical_var_est<-function(temp_level){
          return(h_cost[9]*pmax(0,temp_level-train_D_mat[,9])+b_cost[9]*pmax(0,train_D_mat[,9]-temp_level)+price_10_empirical_var_est(pmax(best_decision_10,temp_level-train_D_mat[,9])))
        }
        
        res_9<-solnp(pars = mean(train_D_mat[,9]),fun = price_9,LB = 0,control = list(trace=0))
        best_decision_9<-res_9$pars
        
        
        price_8<-function(alpha){
          temp_level<-alpha
          return(sum(h_cost[8]*pmax(0,temp_level-train_D_mat[,8])+b_cost[8]*pmax(0,train_D_mat[,8]-temp_level))+price_9_empirical(pmax(best_decision_9,temp_level-train_D_mat[,8])) )
        }
        
        price_8_empirical<-function(temp_level){
          return(sum(h_cost[8]*pmax(0,temp_level-train_D_mat[,8])+b_cost[8]*pmax(0,train_D_mat[,8]-temp_level))+price_9_empirical(pmax(best_decision_9,temp_level-train_D_mat[,8])) )
        }
        
        price_8_empirical_var_est<-function(temp_level){
          return(h_cost[8]*pmax(0,temp_level-train_D_mat[,8])+b_cost[8]*pmax(0,train_D_mat[,8]-temp_level)+price_9_empirical_var_est(pmax(best_decision_9,temp_level-train_D_mat[,8])))
        }
        
        res_8<-solnp(pars = mean(train_D_mat[,8]),fun = price_8,LB=0,control = list(trace=0))
        best_decision_8<-res_8$pars
        
        
        
        price_7<-function(alpha){
          temp_level<-alpha
          return(sum(h_cost[7]*pmax(0,temp_level-train_D_mat[,7])+b_cost[7]*pmax(0,train_D_mat[,7]-temp_level))+price_8_empirical(pmax(best_decision_8,temp_level-train_D_mat[,7])) )
        }
        
        price_7_empirical<-function(temp_level){
          return(sum(h_cost[7]*pmax(0,temp_level-train_D_mat[,7])+b_cost[7]*pmax(0,train_D_mat[,7]-temp_level))+price_8_empirical(pmax(best_decision_8,temp_level-train_D_mat[,7])) )
        }
        
        price_7_empirical_var_est<-function(temp_level){
          return(h_cost[7]*pmax(0,temp_level-train_D_mat[,7])+b_cost[7]*pmax(0,train_D_mat[,7]-temp_level)+price_8_empirical_var_est(pmax(best_decision_8,temp_level-train_D_mat[,7])))
        }
        
        res_7<-solnp(pars = mean(train_D_mat[,7]),fun = price_7,LB=0,control = list(trace=0))
        best_decision_7<-res_7$pars
        
        
        price_6<-function(alpha){
          temp_level<-alpha
          return(sum(h_cost[6]*pmax(0,temp_level-train_D_mat[,6])+b_cost[6]*pmax(0,train_D_mat[,6]-temp_level))+price_7_empirical(pmax(best_decision_7,temp_level-train_D_mat[,6])) )
        }
        
        price_6_empirical<-function(temp_level){
          return(sum(h_cost[6]*pmax(0,temp_level-train_D_mat[,6])+b_cost[6]*pmax(0,train_D_mat[,6]-temp_level))+price_7_empirical(pmax(best_decision_7,temp_level-train_D_mat[,6])) )
        }
        
        price_6_empirical_var_est<-function(temp_level){
          return(h_cost[6]*pmax(0,temp_level-train_D_mat[,6])+b_cost[6]*pmax(0,train_D_mat[,6]-temp_level)+price_7_empirical_var_est(pmax(best_decision_7,temp_level-train_D_mat[,6])))
        }
        
        res_6<-solnp(pars = mean(train_D_mat[,6]),fun = price_6,LB=0,control = list(trace=0))
        best_decision_6<-res_6$pars
        
        
        
        price_5<-function(alpha){
          temp_level<-alpha
          return(sum(h_cost[5]*pmax(0,temp_level-train_D_mat[,5])+b_cost[5]*pmax(0,train_D_mat[,5]-temp_level))+price_6_empirical(pmax(best_decision_6,temp_level-train_D_mat[,5])) )
        }
        
        price_5_empirical<-function(temp_level){
          return(sum(h_cost[5]*pmax(0,temp_level-train_D_mat[,5])+b_cost[5]*pmax(0,train_D_mat[,5]-temp_level))+price_6_empirical(pmax(best_decision_6,temp_level-train_D_mat[,5])) )
        }
        
        price_5_empirical_var_est<-function(temp_level){
          return(h_cost[5]*pmax(0,temp_level-train_D_mat[,5])+b_cost[5]*pmax(0,train_D_mat[,5]-temp_level)+price_6_empirical_var_est(pmax(best_decision_6,temp_level-train_D_mat[,5])))
        }
        
        res_5<-solnp(pars = mean(train_D_mat[,5]),fun = price_5,LB=0,control = list(trace=0))
        best_decision_5<-res_5$pars
        
        
        price_4<-function(alpha){
          temp_level<-alpha
          return(sum(h_cost[4]*pmax(0,temp_level-train_D_mat[,4])+b_cost[4]*pmax(0,train_D_mat[,4]-temp_level))+price_5_empirical(pmax(best_decision_5,temp_level-train_D_mat[,4])) )
        }
        
        price_4_empirical<-function(temp_level){
          return(sum(h_cost[4]*pmax(0,temp_level-train_D_mat[,4])+b_cost[4]*pmax(0,train_D_mat[,4]-temp_level))+price_5_empirical(pmax(best_decision_5,temp_level-train_D_mat[,4])) )
        }
        
        price_4_empirical_var_est<-function(temp_level){
          return(h_cost[4]*pmax(0,temp_level-train_D_mat[,4])+b_cost[4]*pmax(0,train_D_mat[,4]-temp_level)+price_5_empirical_var_est(pmax(best_decision_5,temp_level-train_D_mat[,4])))
        }
        
        res_4<-solnp(pars = mean(train_D_mat[,4]),fun = price_4,LB=0,control = list(trace=0))
        best_decision_4<-res_4$pars
        
        
        
        
        
        price_3<-function(alpha){
          temp_level<-alpha
          return(sum(h_cost[3]*pmax(0,temp_level-train_D_mat[,3])+b_cost[3]*pmax(0,train_D_mat[,3]-temp_level))+price_4_empirical(pmax(best_decision_4,temp_level-train_D_mat[,3])) )
        }
        
        price_3_empirical<-function(temp_level){
          return(sum(h_cost[3]*pmax(0,temp_level-train_D_mat[,3])+b_cost[3]*pmax(0,train_D_mat[,3]-temp_level))+price_4_empirical(pmax(best_decision_4,temp_level-train_D_mat[,3])) )
        }
        
        price_3_empirical_var_est<-function(temp_level){
          return(h_cost[3]*pmax(0,temp_level-train_D_mat[,3])+b_cost[3]*pmax(0,train_D_mat[,3]-temp_level)+price_4_empirical_var_est(pmax(best_decision_4,temp_level-train_D_mat[,3])))
        }
        
        res_3<-solnp(pars = mean(train_D_mat[,3]),fun = price_3,LB=0,control = list(trace=0))
        best_decision_3<-res_3$pars
        
        
        
        
        price_2<-function(alpha){
          temp_level<-alpha
          return(sum(h_cost[2]*pmax(0,temp_level-train_D_mat[,2])+b_cost[2]*pmax(0,train_D_mat[,2]-temp_level))+price_3_empirical(pmax(temp_level-train_D_mat[,2],best_decision_3)))
        }
        
        price_2_empirical<-function(temp_level){
          return(sum(h_cost[2]*pmax(0,temp_level-train_D_mat[,2])+b_cost[2]*pmax(0,train_D_mat[,2]-temp_level))+price_3_empirical(pmax(temp_level-train_D_mat[,2],best_decision_3)))
        }
        
        price_2_empirical_var_est<-function(temp_level){
          return(h_cost[2]*pmax(0,temp_level-train_D_mat[,2])+b_cost[2]*pmax(0,train_D_mat[,2]-temp_level)+price_3_empirical_var_est(pmax(temp_level-train_D_mat[,2],best_decision_3)))
        }
        
        res_2<-solnp(pars = mean(train_D_mat[,2]),fun = price_2,LB=0,control = list(trace=0))
        best_decision_2<-res_2$pars
        
        
        price_1<-function(temp_level){
          return(sum(h_cost[1]*pmax(0,temp_level-train_D_mat[,1])+b_cost[1]*pmax(0,train_D_mat[,1]-temp_level))+price_2_empirical(pmax(temp_level-train_D_mat[,1],best_decision_2)))
        }
        price_1_empirical_var_est<-function(temp_level){
          return((h_cost[1]*pmax(0,temp_level-train_D_mat[,1])+b_cost[1]*pmax(0,train_D_mat[,1]-temp_level))+price_2_empirical_var_est(pmax(temp_level-train_D_mat[,1],best_decision_2)))
        }
        res_1<-optimize(f=price_1,interval = c(0,50))
        
        best_decision_1<-res_1$minimum
        
        test_decision_naive <- c(best_decision_1, best_decision_2, best_decision_3, best_decision_4, best_decision_5, best_decision_6, best_decision_7, best_decision_8, best_decision_9, best_decision_10)
      }
      vanilla_time <- (proc.time() - start_time)["elapsed"]
      
      test_cost_our <- 0
      test_cost_naive <- 0
      test_cost_oracle <- 0
      temp_level_our <- test_decision_our[1]
      temp_level_naive <- test_decision_naive[1]
      temp_level_oracle <- test_decision_oracle[1]
      for(qqqq in 1:T){
        test_cost_our<-test_cost_our+(h_cost[qqqq]*pmax(0,temp_level_our-test_D_mat[qqqq])+b_cost[qqqq]*pmax(0,test_D_mat[qqqq]-temp_level_our))
        test_cost_naive<-test_cost_naive+(h_cost[qqqq]*pmax(0,temp_level_naive-test_D_mat[qqqq])+b_cost[qqqq]*pmax(0,test_D_mat[qqqq]-temp_level_naive))
        test_cost_oracle<-test_cost_oracle+(h_cost[qqqq]*pmax(0,temp_level_oracle-test_D_mat[qqqq])+b_cost[qqqq]*pmax(0,test_D_mat[qqqq]-temp_level_oracle))
        
        if(qqqq<T){
          temp_level_our<-pmax(temp_level_our-test_D_mat[qqqq],test_decision_our[qqqq+1])
          temp_level_naive<-pmax(temp_level_naive-test_D_mat[qqqq],test_decision_naive[qqqq+1])
          temp_level_oracle<-pmax(temp_level_oracle-test_D_mat[qqqq],test_decision_oracle[qqqq+1])
        }
      }
      
      final_res <- list(test_cost_our, test_cost_naive, test_cost_oracle,
                        our_time, our_CV_time, vanilla_time)
      
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

save(res,file = "T equals 10 without covariates.rdata")