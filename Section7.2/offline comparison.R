rm(list=ls())
library(parallel)
simulation<-function(iter){
  my_simulation_procedure<-function(iter){
    load("optimal_policy_T10.rdata")
    library(splines2)
    library(pracma)
    library(Rsolnp)
    n<-50
    T<-10
    X_mat<-matrix(0,n,T)
    D_mat<-matrix(0,n,T)
    
    set.seed(iter+1000)
    X1=runif(n,1,5)
    D1=2*X1+runif(n,-0.5,0.5)
    X2=0.2*X1+runif(n,-0.5,0.5)
    D2=0.2+0.1*(D1+X2)+runif(n,-0.2,0.2)
    X3=0.1*(D1 + X2)+runif(n,-0.2,0.2)
    D3=1+D2+0.5*X3+runif(n,-0.2,0.2)
    X4=0.1*(D2 + 0.5 * X3)+runif(n,-0.2,0.2)
    D4=1+D3+0.5*X4+runif(n,-0.2,0.2)
    X5=0.1*(D3 + 0.5 * X4)+runif(n,-0.2,0.2)
    D5=1+D4+0.5*X5+runif(n,-0.2,0.2)
    X6=0.1*(D4 + 0.5 * X5)+runif(n,-0.2,0.2)
    D6=1+D5+0.5*X6+runif(n,-0.2,0.2)
    X7=0.1*(D5 + 0.5 * X6)+runif(n,-0.2,0.2)
    D7=1+D6+0.5*X7+runif(n,-0.2,0.2)
    X8=0.1*(D6 + 0.5 * X7)+runif(n,-0.2,0.2)
    D8=1+D7+0.5*X8+runif(n,-0.2,0.2)
    X9=0.1*(D7 + 0.5 * X8)+runif(n,-0.2,0.2)
    D9=1+D8+0.5*X9+runif(n,-0.2,0.2)
    X10=0.1*(D8 + 0.5 * X9)+runif(n,-0.2,0.2)
    D10=1+D9+0.5*X10+runif(n,-0.2,0.2)
    
    ###for X2~X10, I treat its range as [-0.5, 2]
    
    X_mat[,1]<-X1
    X_mat[,2]<-X2
    X_mat[,3]<-X3
    X_mat[,4]<-X4
    X_mat[,5]<-X5
    X_mat[,6]<-X6
    X_mat[,7]<-X7
    X_mat[,8]<-X8
    X_mat[,9]<-X9
    X_mat[,10]<-X10
    
    D_mat[,1]<-D1
    D_mat[,2]<-D2
    D_mat[,3]<-D3
    D_mat[,4]<-D4
    D_mat[,5]<-D5
    D_mat[,6]<-D6
    D_mat[,7]<-D7
    D_mat[,8]<-D8
    D_mat[,9]<-D9
    D_mat[,10]<-D10
    
    
    
    
    X_mat_test<-matrix(0,n,T)
    D_mat_test<-matrix(0,n,T)
    
    set.seed(iter+10000)
    X1=runif(n,1,5)
    D1=2*X1+runif(n,-0.5,0.5)
    X2=0.2*X1+runif(n,-0.5,0.5)
    D2=0.2+0.1*(D1+X2)+runif(n,-0.2,0.2)
    X3=0.1*(D1 + X2)+runif(n,-0.2,0.2)
    D3=1+D2+0.5*X3+runif(n,-0.2,0.2)
    X4=0.1*(D2 + 0.5 * X3)+runif(n,-0.2,0.2)
    D4=1+D3+0.5*X4+runif(n,-0.2,0.2)
    X5=0.1*(D3 + 0.5 * X4)+runif(n,-0.2,0.2)
    D5=1+D4+0.5*X5+runif(n,-0.2,0.2)
    X6=0.1*(D4 + 0.5 * X5)+runif(n,-0.2,0.2)
    D6=1+D5+0.5*X6+runif(n,-0.2,0.2)
    X7=0.1*(D5 + 0.5 * X6)+runif(n,-0.2,0.2)
    D7=1+D6+0.5*X7+runif(n,-0.2,0.2)
    X8=0.1*(D6 + 0.5 * X7)+runif(n,-0.2,0.2)
    D8=1+D7+0.5*X8+runif(n,-0.2,0.2)
    X9=0.1*(D7 + 0.5 * X8)+runif(n,-0.2,0.2)
    D9=1+D8+0.5*X9+runif(n,-0.2,0.2)
    X10=0.1*(D8 + 0.5 * X9)+runif(n,-0.2,0.2)
    D10=1+D9+0.5*X10+runif(n,-0.2,0.2)
    
    X_mat_test[,1]<-X1
    X_mat_test[,2]<-X2
    X_mat_test[,3]<-X3
    X_mat_test[,4]<-X4
    X_mat_test[,5]<-X5
    X_mat_test[,6]<-X6
    X_mat_test[,7]<-X7
    X_mat_test[,8]<-X8
    X_mat_test[,9]<-X9
    X_mat_test[,10]<-X10
    
    D_mat_test[,1]<-D1
    D_mat_test[,2]<-D2
    D_mat_test[,3]<-D3
    D_mat_test[,4]<-D4
    D_mat_test[,5]<-D5
    D_mat_test[,6]<-D6
    D_mat_test[,7]<-D7
    D_mat_test[,8]<-D8
    D_mat_test[,9]<-D9
    D_mat_test[,10]<-D10
    
    
    
    b_cost<-rep(1,T)
    b_cost[1]<-5
    h_cost<-rep(1,T)
    h_cost[2]<-5
    
    my_initial<-0.1
    
    lower_knots<-c(1,0.5,-0.9,0,0.8,1.5,2.5,3.5,4.5,5.5)
    upper_knots<-c(5,15,6,7.5,9,10,11.5,13,14.5,16)
    
    ########################################################################################################################################################################################################################################
    ########################################################################################################################################################################################################################################
    ptm<-proc.time()
    m_n_candidate<-1:10
    CV_mat<-matrix(0,5,length(m_n_candidate))
    my_split<-1:n
    my_split<-matrix(my_split,5,n/5,byrow = TRUE)
    for(fold in 1:5){
      test_D_mat<-D_mat[c(my_split[fold,]),]
      test_X_mat<-X_mat[c(my_split[fold,]),]
      
      train_D_mat<-D_mat[c(my_split[-fold,]),]
      train_X_mat<-X_mat[c(my_split[-fold,]),]
      
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
        
        # temp_x<-seq(lower_knots[1],upper_knots[1],0.01)
        # temp_base<-bSpline(x=temp_x,knots = knots_1,degree = 3,Boundary.knots = c(lower_knots[1],upper_knots[1]),intercept = TRUE )
        # S1_est<-c(temp_base%*%res_1$par)
        # plot(temp_x,S1_est,type = "l",lwd=3)
        # points(temp_x,temp_x,lwd=3,lty=2,type = "l",col="red")
        
        
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
        print(c(fold,ppp,beta_2,beta_3,beta_4,beta_5,beta_6,beta_7,beta_8,beta_9,beta_10))
      }
    }
    CVtime<-as.numeric((proc.time()-ptm)[3])
    CV_seq<-apply(CV_mat,MARGIN = 2,sum)
    
    ptm<-proc.time()
    m_n<-m_n_candidate[(which.min(CV_seq))[1]]
    
    train_D_mat<-D_mat
    train_X_mat<-X_mat
    
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
    
    
    test_D_mat <- D_mat_test
    test_X_mat <- X_mat_test
    
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
    
    test_decision_ZX<-rbind(test_decision_1,test_decision_2,test_decision_3,test_decision_4,test_decision_5,test_decision_6,test_decision_7,test_decision_8,test_decision_9,test_decision_10)
    
    
    Runtime<-as.numeric((proc.time()-ptm)[3])
    
    
    ####this is for the optimal policy
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
    
    lower_knots<-c(1,0.5,-0.9,0,0.8,1.5,2.5,3.5,4.5,5.5)
    upper_knots<-c(5,15,6,7.5,9,10,11.5,13,14.5,16)
    
    basis_test_1<-bSpline(x=test_X_mat[,1],knots = knots_1,degree = 3,Boundary.knots = c(lower_knots[1],upper_knots[1]),intercept = TRUE )
    test_decision_1<-c(basis_test_1%*%pars1)
    
    temp_x<-c(test_X_mat[,2]*pars2[1]+test_D_mat[,1])
    temp_base<-bSpline(x=temp_x,knots = knots_2,degree = 3,Boundary.knots = c(lower_knots[2],upper_knots[2]),intercept = TRUE )
    test_decision_2<-c(temp_base%*%pars2[-1])
    
    temp_x<-c(test_X_mat[,3]*pars3[1]+test_D_mat[,2])
    temp_base<-bSpline(x=temp_x,knots = knots_3,degree = 3,Boundary.knots = c(lower_knots[3],upper_knots[3]),intercept = TRUE )
    test_decision_3<-c(temp_base%*%pars3[-1])
    
    temp_x<-c(test_X_mat[,4]*pars4[1]+test_D_mat[,3])
    temp_base<-bSpline(x=temp_x,knots = knots_4,degree = 3,Boundary.knots = c(lower_knots[4],upper_knots[4]),intercept = TRUE )
    test_decision_4<-c(temp_base%*%pars4[-1])
    
    temp_x<-c(test_X_mat[,5]*pars5[1]+test_D_mat[,4])
    temp_base<-bSpline(x=temp_x,knots = knots_5,degree = 3,Boundary.knots = c(lower_knots[5],upper_knots[5]),intercept = TRUE )
    test_decision_5<-c(temp_base%*%pars5[-1])
    
    temp_x<-c(test_X_mat[,6]*pars6[1]+test_D_mat[,5])
    temp_base<-bSpline(x=temp_x,knots = knots_6,degree = 3,Boundary.knots = c(lower_knots[6],upper_knots[6]),intercept = TRUE )
    test_decision_6<-c(temp_base%*%pars6[-1])
    
    temp_x<-c(test_X_mat[,7]*pars7[1]+test_D_mat[,6])
    temp_base<-bSpline(x=temp_x,knots = knots_7,degree = 3,Boundary.knots = c(lower_knots[7],upper_knots[7]),intercept = TRUE )
    test_decision_7<-c(temp_base%*%pars7[-1])
    
    temp_x<-c(test_X_mat[,8]*pars8[1]+test_D_mat[,7])
    temp_base<-bSpline(x=temp_x,knots = knots_8,degree = 3,Boundary.knots = c(lower_knots[8],upper_knots[8]),intercept = TRUE )
    test_decision_8<-c(temp_base%*%pars8[-1])
    
    temp_x<-c(test_X_mat[,9]*pars9[1]+test_D_mat[,8])
    temp_base<-bSpline(x=temp_x,knots = knots_9,degree = 3,Boundary.knots = c(lower_knots[9],upper_knots[9]),intercept = TRUE )
    test_decision_9<-c(temp_base%*%pars9[-1])
    
    temp_x<-c(test_X_mat[,10]*pars10[1]+test_D_mat[,9])
    temp_base<-bSpline(x=temp_x,knots = knots_10,degree = 3,Boundary.knots = c(lower_knots[10],upper_knots[10]),intercept = TRUE )
    test_decision_10<-c(temp_base%*%pars10[-1])
    
    test_decision_oracle<-rbind(test_decision_1,test_decision_2,test_decision_3,test_decision_4,test_decision_5,test_decision_6,test_decision_7,test_decision_8,test_decision_9,test_decision_10)
    
    test_cost_ZX<-rep(0, n)
    test_cost_oracle<-rep(0, n)
    temp_level_ZX<-test_decision_ZX[1, ]
    temp_level_oracle <- test_decision_oracle[1, ]
    for(qqqq in 1:T){
      test_cost_ZX<-test_cost_ZX+(h_cost[qqqq]*pmax(0,temp_level_ZX-test_D_mat[,qqqq])+b_cost[qqqq]*pmax(0,test_D_mat[,qqqq]-temp_level_ZX))
      test_cost_oracle<-test_cost_oracle+(h_cost[qqqq]*pmax(0,temp_level_oracle-test_D_mat[,qqqq])+b_cost[qqqq]*pmax(0,test_D_mat[,qqqq]-temp_level_oracle))
      
      if(qqqq<T){
        temp_level_ZX<-pmax(temp_level_ZX-test_D_mat[,qqqq],test_decision_ZX[qqqq+1,])
        temp_level_oracle<-pmax(temp_level_oracle-test_D_mat[,qqqq],test_decision_oracle[qqqq+1,])
      }
    }
    
    sum(test_cost_ZX - test_cost_oracle)
    
    return(list(CVtime,Runtime,test_cost_ZX,test_cost_oracle,X_mat,D_mat,X_mat_test,D_mat_test))
  }
  res<-try(my_simulation_procedure(iter),silent = TRUE)
  return(res)
}


clnum<-20#detectCores()/2
cl <- makeCluster(getOption("cl.cores", clnum))
res<-parLapply(cl, 1:500,  simulation)
stopCluster(cl)

save(res,file = "T equals 10 with covariates offline.rdata")