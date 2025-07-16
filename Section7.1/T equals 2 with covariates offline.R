rm(list=ls())
library(parallel)
simulation<-function(iter){
  my_simulation_procedure<-function(iter){
    library(splines2)
    library(pracma)
    library(Rsolnp)
    library(rpart)
    library(rpart.plot)
    library(randomForest)
    n<-200
    T<-3
    X_mat<-cbind(runif(n),runif(n),runif(n))
    D_mat<-matrix(0,n,T)
    set.seed(iter)
    X1=runif(n,1,5)
    D1=2*X1+runif(n,-0.5,0.5)###1.5 10.5
    X2=0.2*X1+runif(n,-0.5,0.5)##-0.3  1.5
    D2=0.2+0.1*(D1+X2)+runif(n,-0.2,0.2)
    
    X_mat[,1]<-X1
    X_mat[,2]<-X2
    D_mat[,1]<-D1
    D_mat[,2]<-D2
    
    X_mat_test<-cbind(runif(n),runif(n),runif(n))
    D_mat_test<-matrix(0,n,T)
    set.seed(iter+1000)
    X1=runif(n,1,5)
    D1=2*X1+runif(n,-0.5,0.5)###1.5 10.5
    X2=0.2*X1+runif(n,-0.5,0.5)##-0.3  1.5
    D2=0.2+0.1*(D1+X2)+runif(n,-0.2,0.2)
    
    X_mat_test[,1]<-X1
    X_mat_test[,2]<-X2
    D_mat_test[,1]<-D1
    D_mat_test[,2]<-D2
    
    X_mat_test<-X_mat_test[1:10,]
    D_mat_test<-D_mat_test[1:10,]
    
    
    n<-50
    X_mat<-X_mat[1:n,]
    D_mat<-D_mat[1:n,]
    X_mat_test<-X_mat_test[1:10,]
    D_mat_test<-D_mat_test[1:10,]
    
    b_cost<-c(5,1,1)
    h_cost<-c(1,5,1)
    ########################################################################################################################################################################################################################################
    ########################################################################################################################################################################################################################################
    gaussian_kernel<-function(x1,x2,sigma=1){
      return(exp(-(x1-x2)^2/sigma^2))
    }
    gaussian_kernel_high<-function(x1,x2,sigma=1){
      return(exp(-sum((x1-x2)^2)/sigma^2))
    }
    ptm<-proc.time()
    sigma_candidate<-seq(2,20,2)
    lambda_candidate<-5^(-c(1:10))
    CV_array<-array(0,dim = c(length(sigma_candidate),length(lambda_candidate),5))
    my_split<-1:n
    my_split<-matrix(my_split,5,n/5,byrow = TRUE)
    for(fold in 1:5){
      test_X<-D_mat[c(my_split[fold,]),1]
      test_Y<-D_mat[c(my_split[fold,]),2]
      test_Z1<-X_mat[c(my_split[fold,]),1]
      test_Z2<-X_mat[c(my_split[fold,]),2]
      
      train_X<-D_mat[c(my_split[-fold,]),1]
      train_Y<-D_mat[c(my_split[-fold,]),2]
      train_Z1<-X_mat[c(my_split[-fold,]),1]
      train_Z2<-X_mat[c(my_split[-fold,]),2]
      for(ppp in 1:length(sigma_candidate)){
        for(qqq in 1:length(lambda_candidate)){
          k_mat<-matrix(0,0.8*n,0.8*n)
          for(i in 1:(0.8*n)){
            k_mat[i,]<-gaussian_kernel(x1=train_Z1[i],x2=train_Z1,sigma = sigma_candidate[ppp])
          }
          
          k_mat_high<-matrix(0,0.8*n,0.8*n)
          for(i in 1:(0.8*n)){
            for(j in 1:(0.8*n)){
              k_mat_high[i,j]<-gaussian_kernel_high(x1=c(train_Z2[i],train_X[i]),x2=c(train_Z2[j],train_X[j]),sigma = sigma_candidate[ppp])
            }
          }
          
          keropt_cost<-function(a_vec,mylambda=numeric(0)){
            initial_level<-c(k_mat%*%a_vec[1:(0.8*n)])
            temp_level<-c(k_mat_high%*%a_vec[-c(1:(0.8*n))])
            return(sum(h_cost[1]*pmax(0,initial_level-train_X)+b_cost[1]*pmax(0,train_X-initial_level)+h_cost[2]*pmax(0,pmax(initial_level-train_X,temp_level)-train_Y)+b_cost[2]*pmax(0,train_Y-pmax(initial_level-train_X,temp_level)))+c(0.8*n*mylambda*t(a_vec[1:(0.8*n)])%*%k_mat%*%a_vec[1:(0.8*n)])+c(0.8*n*mylambda*t(a_vec[-c(1:(0.8*n))])%*%k_mat_high%*%a_vec[-c(1:(0.8*n))])+0.8*n*1e4*max(0,-min(initial_level,temp_level)))
          }
          
          opt_res<-solnp(pars = rep(0.1,0.8*n*2),fun = keropt_cost,mylambda=lambda_candidate[qqq],control = list(trace=0))
          #,ineqfun = inequality_constraint,ineqLB = rep(0,2),ineqUB = rep(Inf,2)
          
          k_mat_test<-matrix(0,0.2*n,0.8*n)
          for(i in 1:(0.2*n)){
            k_mat_test[i,]<-gaussian_kernel(x1=test_Z1[i],x2=train_Z1,sigma = sigma_candidate[ppp])
          }
          test_decision_1<-c(k_mat_test%*%opt_res$pars[1:(0.8*n)])
          
          k_mat_test_high<-matrix(0,0.2*n,0.8*n)
          for(i in 1:(0.2*n)){
            for(j in 1:(0.8*n)){
              k_mat_test_high[i,j]<-gaussian_kernel_high(x1=c(test_Z2[i],test_X[i]),x2=c(train_Z2[j],train_X[j]),sigma = sigma_candidate[ppp])
            }
          }
          test_decision_2<-c(k_mat_test_high%*%opt_res$pars[-c(1:(0.8*n))])
          
          temp_res<-0
          for(ooo in 1:length(test_X)){
            temp_res<-temp_res+h_cost[1]*max(0,test_decision_1[ooo]-test_X[ooo])+b_cost[1]*max(0,test_X[ooo]-test_decision_1[ooo])
            temp_res<-temp_res+h_cost[2]*max(0,max(test_decision_1[ooo]-test_X[ooo],test_decision_2[ooo])-test_Y[ooo])+b_cost[2]*max(0,test_Y[ooo]-max(test_decision_1[ooo]-test_X[ooo],test_decision_2[ooo]))
          }
          
          CV_array[ppp,qqq,fold]<-temp_res
          print(c(fold,ppp,qqq,opt_res$convergence))
        }
      }
    }
    CV_mat<-apply(CV_array,MARGIN = c(1,2),sum)
    CVtime_RKHS<-as.numeric((proc.time()-ptm)[3])
    ptm<-proc.time()
    CV_mat_RKHS<-CV_mat
    ind_selected<-which(CV_mat == min(CV_mat), arr.ind = TRUE)[1,]
    sigma<-sigma_candidate[ind_selected[1]] 
    mylambda<-lambda_candidate[ind_selected[2]]
    ind_selected_RKHS<-ind_selected
    
    train_X<-D_mat[,1]
    train_Y<-D_mat[,2]
    train_Z1<-X_mat[,1]
    train_Z2<-X_mat[,2]
    
    
    test_X<-D_mat_test[,1]
    test_Y<-D_mat_test[,2]
    test_Z1<-X_mat_test[,1]
    test_Z2<-X_mat_test[,2]
    
    
    
    k_mat<-matrix(0,n,n)
    for(i in 1:(n)){
      k_mat[i,]<-gaussian_kernel(x1=train_Z1[i],x2=train_Z1,sigma = sigma)
    }
    
    k_mat_high<-matrix(0,n,n)
    for(i in 1:(n)){
      for(j in 1:(n)){
        k_mat_high[i,j]<-gaussian_kernel_high(x1=c(train_Z2[i],train_X[i]),x2=c(train_Z2[j],train_X[j]),sigma = sigma)
      }
    }
    
    keropt_cost<-function(a_vec,mylambda=numeric(0)){
      initial_level<-c(k_mat%*%a_vec[1:(n)])
      temp_level<-c(k_mat_high%*%a_vec[-c(1:(n))])
      return(sum(h_cost[1]*pmax(0,initial_level-train_X)+b_cost[1]*pmax(0,train_X-initial_level)+h_cost[2]*pmax(0,pmax(initial_level-train_X,temp_level)-train_Y)+b_cost[2]*pmax(0,train_Y-pmax(initial_level-train_X,temp_level)))+c(n*mylambda*t(a_vec[1:(n)])%*%k_mat%*%a_vec[1:(n)])+c(n*mylambda*t(a_vec[-c(1:(n))])%*%k_mat_high%*%a_vec[-c(1:(n))])+n*1e4*max(0,-min(initial_level,temp_level)))
    }
    
    opt_res<-solnp(pars = rep(0.1,n*2),fun = keropt_cost,mylambda=mylambda,control = list(trace=0))
    #,ineqfun = inequality_constraint,ineqLB = rep(0,2),ineqUB = rep(Inf,2)
    
    k_mat_test<-matrix(0,10,n)
    for(i in 1:(10)){
      k_mat_test[i,]<-gaussian_kernel(x1=test_Z1[i],x2=train_Z1,sigma = sigma)
    }
    test_decision_1<-c(k_mat_test%*%opt_res$pars[1:(n)])
    
    k_mat_test_high<-matrix(0,10,n)
    for(i in 1:(10)){
      for(j in 1:(n)){
        k_mat_test_high[i,j]<-gaussian_kernel_high(x1=c(test_Z2[i],test_X[i]),x2=c(train_Z2[j],train_X[j]),sigma = sigma)
      }
    }
    test_decision_2<-c(k_mat_test_high%*%opt_res$pars[-c(1:(n))])
    
    temp_res<-rep(0,10)
    for(ooo in 1:length(test_X)){
      temp_res[ooo]<-temp_res[ooo]+h_cost[1]*max(0,test_decision_1[ooo]-test_X[ooo])+b_cost[1]*max(0,test_X[ooo]-test_decision_1[ooo])
      temp_res[ooo]<-temp_res[ooo]+h_cost[2]*max(0,max(test_decision_1[ooo]-test_X[ooo],test_decision_2[ooo])-test_Y[ooo])+b_cost[2]*max(0,test_Y[ooo]-max(test_decision_1[ooo]-test_X[ooo],test_decision_2[ooo]))
    }
    cost_RKHS<-temp_res
    time_RKHS<-as.numeric((proc.time()-ptm)[3])
    
    
    
    ########################################################################################################################################################################################################################################
    ########################################################################################################################################################################################################################################
    ptm<-proc.time()
    m_n_candidate<-1:10
    CV_mat<-matrix(0,5,length(m_n_candidate))
    my_split<-1:n
    my_split<-matrix(my_split,5,n/5,byrow = TRUE)
    for(fold in 1:5){
      test_X<-D_mat[c(my_split[fold,]),1]
      test_Y<-D_mat[c(my_split[fold,]),2]
      test_Z1<-X_mat[c(my_split[fold,]),1]
      test_Z2<-X_mat[c(my_split[fold,]),2]
      
      train_X<-D_mat[c(my_split[-fold,]),1]
      train_Y<-D_mat[c(my_split[-fold,]),2]
      train_Z1<-X_mat[c(my_split[-fold,]),1]
      train_Z2<-X_mat[c(my_split[-fold,]),2]
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
        
        price_2_var_est<-function(temp_level_2){
          return(h_cost[2]*pmax(0,temp_level_2-train_Y)+b_cost[2]*pmax(0,train_Y-temp_level_2))
        }
        
        res_2<-solnp(pars = c(0.1,rep(0.1,m_n+4)),fun = price_2,LB=c(0,rep(0,m_n+4)),UB=c(2,rep(Inf,m_n+4)),control = list(trace=0))
        basis_mat_2<-bSpline(x=train_X+res_2$pars[1]*train_Z2,knots = knots_2,degree = 3,Boundary.knots = c(0.9,13.5),intercept = TRUE )
        best_decision_2<-c(basis_mat_2%*%res_2$pars[-1])
        
        knots_1<-seq(from=1,to=5,length.out=m_n+2)[-c(1,m_n+2)]
        basis_mat_1<-bSpline(x=train_Z1,knots = knots_1,degree = 3,Boundary.knots = c(1,5),intercept = TRUE )
        price_1<-function(alpha){
          temp_level<-c(basis_mat_1%*%alpha)
          return(sum(h_cost[1]*pmax(0,temp_level-train_X)+b_cost[1]*pmax(0,train_X-temp_level))+price_2_empirical(pmax(temp_level-train_X,best_decision_2)))
        }
        price_1_var_est<-function(alpha){
          temp_level<-c(basis_mat_1%*%alpha)
          return(h_cost[1]*pmax(0,temp_level-train_X)+b_cost[1]*pmax(0,train_X-temp_level)+price_2_var_est(pmax(temp_level-train_X,best_decision_2)))
        }
        
        res_1<-solnp(pars = rep(0.1,m_n+4),fun = price_1,LB=rep(0,m_n+4),control = list(trace=0))
        
        basis_test_1<-bSpline(x=test_Z1,knots = knots_1,degree = 3,Boundary.knots = c(1,5),intercept = TRUE )
        test_decision_1<-c(basis_test_1%*%res_1$pars)
        
        temp_x<-c(test_Z2*res_2$pars[1]+test_X)
        temp_base<-bSpline(x=temp_x,knots = knots_2,degree = 3,Boundary.knots = c(0.9,13.5),intercept = TRUE )
        test_decision_2<-c(temp_base%*%res_2$pars[-1])
        
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
    CVtime_ZX<-as.numeric((proc.time()-ptm)[3])
    ptm<-proc.time()
    CV_seq_ZX<-CV_seq
    m_n<-m_n_candidate[(which.min(CV_seq))[1]]
    
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
    
    price_2_var_est<-function(temp_level_2){
      return(h_cost[2]*pmax(0,temp_level_2-train_Y)+b_cost[2]*pmax(0,train_Y-temp_level_2))
    }
    
    res_2<-solnp(pars = c(0.1,rep(0.1,m_n+4)),fun = price_2,LB=c(0,rep(0,m_n+4)),UB=c(2,rep(Inf,m_n+4)),control = list(trace=0))
    basis_mat_2<-bSpline(x=train_X+res_2$pars[1]*train_Z2,knots = knots_2,degree = 3,Boundary.knots = c(0.9,13.5),intercept = TRUE )
    best_decision_2<-c(basis_mat_2%*%res_2$pars[-1])
    
    knots_1<-seq(from=1,to=5,length.out=m_n+2)[-c(1,m_n+2)]
    basis_mat_1<-bSpline(x=train_Z1,knots = knots_1,degree = 3,Boundary.knots = c(1,5),intercept = TRUE )
    price_1<-function(alpha){
      temp_level<-c(basis_mat_1%*%alpha)
      return(sum(h_cost[1]*pmax(0,temp_level-train_X)+b_cost[1]*pmax(0,train_X-temp_level))+price_2_empirical(pmax(temp_level-train_X,best_decision_2)))
    }
    price_1_var_est<-function(alpha){
      temp_level<-c(basis_mat_1%*%alpha)
      return(h_cost[1]*pmax(0,temp_level-train_X)+b_cost[1]*pmax(0,train_X-temp_level)+price_2_var_est(pmax(temp_level-train_X,best_decision_2)))
    }
    
    res_1<-solnp(pars = rep(0.1,m_n+4),fun = price_1,LB=rep(0,m_n+4),control = list(trace=0))
    
    basis_test_1<-bSpline(x=X_mat_test[,1],knots = knots_1,degree = 3,Boundary.knots = c(1,5),intercept = TRUE )
    test_decision_1<-c(basis_test_1%*%res_1$pars)
    
    temp_x<-c(X_mat_test[,2]*res_2$pars[1]+D_mat_test[,1])
    temp_base<-bSpline(x=temp_x,knots = knots_2,degree = 3,Boundary.knots = c(0.9,13.5),intercept = TRUE )
    test_decision_2<-c(temp_base%*%res_2$pars[-1])
    
    temp_res<-rep(0,10)
    for(ooo in 1:10){
      temp_res[ooo]<-temp_res[ooo]+h_cost[1]*max(0,test_decision_1[ooo]-D_mat_test[ooo,1])+b_cost[1]*max(0,D_mat_test[ooo,1]-test_decision_1[ooo])
      temp_res[ooo]<-temp_res[ooo]+h_cost[2]*max(0,max(test_decision_1[ooo]-D_mat_test[ooo,1],test_decision_2[ooo])-D_mat_test[ooo,2])+b_cost[2]*max(0,D_mat_test[ooo,2]-max(test_decision_1[ooo]-D_mat_test[ooo,1],test_decision_2[ooo]))
    }
    cost_ZX<-temp_res
    time_ZX<-as.numeric((proc.time()-ptm)[3])
    
    ########################################################################################################################################################################################################################################
    ########################################################################################################################################################################################################################################
    euclidean<-function(x,y){
      res<-rep(0,nrow(y))
      for(i in 1:nrow(y)){
        res[i]<-sqrt(sum((x-y[i,])^2)) 
      }
      return(res)
    }
    ptm<-proc.time()
    neighbor_candidate<-round(seq(2,n/3,length.out=5))
    CV_mat<-matrix(0,5,length(neighbor_candidate))
    my_split<-1:n
    my_split<-matrix(my_split,5,n/5,byrow = TRUE)
    for(fold in 1:5){
      test_X<-D_mat[c(my_split[fold,]),1]
      test_Y<-D_mat[c(my_split[fold,]),2]
      test_Z1<-X_mat[c(my_split[fold,]),1]
      test_Z2<-X_mat[c(my_split[fold,]),2]
      
      train_X<-D_mat[c(my_split[-fold,]),1]
      train_Y<-D_mat[c(my_split[-fold,]),2]
      train_Z1<-X_mat[c(my_split[-fold,]),1]
      train_Z2<-X_mat[c(my_split[-fold,]),2]
      
      dist_mat_train<-matrix(0,length(train_X),length(train_X))
      for(xxx in 1:length(train_X)){
        dist_mat_train[xxx,]<-order(euclidean(x=c(train_Z2[xxx],train_X[xxx]),y=cbind(train_Z2,train_X)))
      }
      
      dist_mat_test<-matrix(0,length(test_X),length(train_X))
      for(xxx in 1:length(test_X)){
        dist_mat_test[xxx,]<-order(euclidean(x=c(test_Z2[xxx],test_X[xxx]),y=cbind(train_Z2,train_X)))
      }
      
      dist_mat_test1<-matrix(0,length(test_X),length(train_X))
      for(xxx in 1:length(test_X)){
        dist_mat_test1[xxx,]<-order(abs(test_Z1[xxx]-train_Z1))
      }
      
      for(ppp in 1:length(neighbor_candidate)){
        k_num<-neighbor_candidate[ppp]
        
        price_2<-function(z,x=numeric(0),starting_level=numeric(0)){
          temp_ind<-dist_mat_train[x,]
          temp_trainY<-train_Y[temp_ind[1:k_num]]
          return(mean(h_cost[2]*pmax(0,pmax(starting_level,z)-temp_trainY)+b_cost[2]*pmax(0,temp_trainY-pmax(starting_level,z))))
        }
        
        price_1<-function(z0,x0=numeric(0)){
          temp_ind<-dist_mat_test1[x0,]
          temp_ind<-temp_ind[1:k_num]
          next_cost<-rep(0,k_num)
          for(www in 1:k_num){
            temp_www<-solnp(pars = mean(train_Y),fun = price_2,control = list(trace=0),x=temp_ind[www],starting_level=z0-train_X[temp_ind[www]],LB=0)$values
            next_cost[www]<-temp_www[length(temp_www)]
          }
          temp_trainX<-train_X[temp_ind[1:k_num]]
          return(sum(h_cost[1]*pmax(z0-temp_trainX,0)+b_cost[1]*pmax(temp_trainX-z0,0)+next_cost)) 
        }
        
        price_2_test<-function(z,x=numeric(0),starting_level=numeric(0)){
          temp_ind<-dist_mat_test[x,]
          temp_trainY<-train_Y[temp_ind[1:k_num]]
          return(mean(h_cost[2]*pmax(0,pmax(starting_level,z)-temp_trainY)+b_cost[2]*pmax(0,temp_trainY-pmax(starting_level,z))))
        }
        
        res<-0
        for(rrr in 1:length(test_Z1)){
          res_period1<-solnp(pars = mean(train_X),fun = price_1,LB=0,control = list(trace=0),x0=rrr)
          res_period1<-res_period1$pars
          z_level<-solnp(pars = mean(train_Y),fun = price_2_test,control = list(trace=0),x=rrr,starting_level=res_period1-test_X[rrr],LB=0)$pars
          z_level<-max(z_level,res_period1-test_X[rrr])
          res<-res+h_cost[1]*max(0,res_period1-test_X[rrr])+b_cost[1]*max(0,test_X[rrr]-res_period1)+h_cost[2]*max(0,z_level-test_Y[rrr])+b_cost[2]*max(0,test_Y[rrr]-z_level)
        }
        
        CV_mat[fold,ppp]<-res
        print(c(fold,ppp))
      }
    }
    CV_seq<-apply(CV_mat,MARGIN = 2,sum)
    CVtime_knn<-as.numeric((proc.time()-ptm)[3])
    
    
    CV_seq_knn<-CV_seq
    k_num<-neighbor_candidate[(which.min(CV_seq))[1]]
    
    train_X<-D_mat[,1]
    train_Y<-D_mat[,2]
    train_Z1<-X_mat[,1]
    train_Z2<-X_mat[,2]
    
    
    test_X<-D_mat_test[,1]
    test_Y<-D_mat_test[,2]
    test_Z1<-X_mat_test[,1]
    test_Z2<-X_mat_test[,2]
    
    dist_mat_train<-matrix(0,length(train_X),length(train_X))
    for(xxx in 1:length(train_X)){
      dist_mat_train[xxx,]<-order(euclidean(x=c(train_Z2[xxx],train_X[xxx]),y=cbind(train_Z2,train_X)))
    }
    
    dist_mat_test<-matrix(0,length(test_X),length(train_X))
    for(xxx in 1:length(test_X)){
      dist_mat_test[xxx,]<-order(euclidean(x=c(test_Z2[xxx],test_X[xxx]),y=cbind(train_Z2,train_X)))
    }
    
    dist_mat_test1<-matrix(0,length(test_X),length(train_X))
    for(xxx in 1:length(test_X)){
      dist_mat_test1[xxx,]<-order(abs(test_Z1[xxx]-train_Z1))
    }
    
    price_2<-function(z,x=numeric(0),starting_level=numeric(0)){
      temp_ind<-dist_mat_train[x,]
      temp_trainY<-train_Y[temp_ind[1:k_num]]
      return(mean(h_cost[2]*pmax(0,pmax(starting_level,z)-temp_trainY)+b_cost[2]*pmax(0,temp_trainY-pmax(starting_level,z))))
    }
    
    price_1<-function(z0,x0=numeric(0)){
      temp_ind<-dist_mat_test1[x0,]
      temp_ind<-temp_ind[1:k_num]
      next_cost<-rep(0,k_num)
      for(www in 1:k_num){
        temp_www<-solnp(pars = mean(train_Y),fun = price_2,control = list(trace=0),x=temp_ind[www],starting_level=z0-train_X[temp_ind[www]],LB=0)$values
        next_cost[www]<-temp_www[length(temp_www)]
      }
      temp_trainX<-train_X[temp_ind[1:k_num]]
      return(sum(h_cost[1]*pmax(z0-temp_trainX,0)+b_cost[1]*pmax(temp_trainX-z0,0)+next_cost)) 
    }
    
    price_2_test<-function(z,x=numeric(0),starting_level=numeric(0)){
      temp_ind<-dist_mat_test[x,]
      temp_trainY<-train_Y[temp_ind[1:k_num]]
      return(mean(h_cost[2]*pmax(0,pmax(starting_level,z)-temp_trainY)+b_cost[2]*pmax(0,temp_trainY-pmax(starting_level,z))))
    }
    
    
    
    res<-rep(0,length(test_Z1))
    time_knn<-rep(0,length(test_Z1))
    for(rrr in 1:length(test_Z1)){
      ptm<-proc.time()
      res_period1<-solnp(pars = mean(train_X),fun = price_1,LB=0,control = list(trace=0),x0=rrr)
      res_period1<-res_period1$pars
      z_level<-solnp(pars = mean(train_Y),fun = price_2_test,control = list(trace=0),x=rrr,starting_level=res_period1-test_X[rrr],LB=0)$pars
      z_level<-max(z_level,res_period1-test_X[rrr])
      res[rrr]<-res[rrr]+h_cost[1]*max(0,res_period1-test_X[rrr])+b_cost[1]*max(0,test_X[rrr]-res_period1)+h_cost[2]*max(0,z_level-test_Y[rrr])+b_cost[2]*max(0,test_Y[rrr]-z_level)
      time_knn[rrr]<-as.numeric((proc.time()-ptm)[3])
    }
    cost_knn<-res
    
    
    
    ########################################################################################################################################################################################################################################
    ########################################################################################################################################################################################################################################
    ptm<-proc.time()
    minbucket_candidate<-rev(unique(round(seq(2,n/3,length.out=5)) )) 
    CV_array<-matrix(0,length(minbucket_candidate),5)
    my_split<-1:n
    my_split<-matrix(my_split,5,n/5,byrow = TRUE)
    for(fold in 1:5){
      test_X<-D_mat[c(my_split[fold,]),1]
      test_Y<-D_mat[c(my_split[fold,]),2]
      test_Z1<-X_mat[c(my_split[fold,]),1]
      test_Z2<-X_mat[c(my_split[fold,]),2]
      
      train_X<-D_mat[c(my_split[-fold,]),1]
      train_Y<-D_mat[c(my_split[-fold,]),2]
      train_Z1<-X_mat[c(my_split[-fold,]),1]
      train_Z2<-X_mat[c(my_split[-fold,]),2]
      
      my_data<-as.data.frame(cbind(train_Z1,train_X))
      test_data<-as.data.frame(cbind(test_Z1,test_X))
      names(test_data)<-c("train_Z1","train_X")
      
      
      my_data1<-as.data.frame(cbind(train_Z2,train_X,train_Y))
      test_data1<-as.data.frame(cbind(test_Z2,test_X,test_Y))
      names(test_data1)<-c("train_Z2","train_X","train_Y")
      for(qqq in 1:length(minbucket_candidate)){
        m1<-rpart(train_Y~train_Z2+train_X,data = my_data1,method  = "anova",control = list(xval=5,cp=0,minsplit=2,minbucket=minbucket_candidate[qqq]))
        # rpart.plot(m1)
        
        
        pred_value_train<-predict(m1, newdata = my_data1)
        
        pred_value_test<-predict(m1, newdata = test_data1)
        
        
        price_2<-function(z,fff=numeric(0),starting_level=numeric(0)){
          aaa<-which(pred_value_train==pred_value_train[fff]) 
          temp_trainY<-train_Y[aaa]
          return(mean(h_cost[2]*pmax(0,pmax(z,starting_level)-temp_trainY)+b_cost[2]*pmax(0,temp_trainY-pmax(z,starting_level))))
        }
        
        price_1<-function(z0,weight_vector=numeric(0)){
          next_cost<-rep(0,length(train_X))
          for(www in 1:length(train_X)){
            if(weight_vector[www]>0){
              temp_www<-solnp(pars = mean(train_Y),fun = price_2,control = list(trace=0),fff=www,starting_level=z0-train_X[www],LB=0)$values
              next_cost[www]<-temp_www[length(temp_www)]
            }
          }
          return(sum(weight_vector*(h_cost[1]*pmax(z0-train_X,0)+b_cost[1]*pmax(train_X-z0,0)+next_cost))) 
        }
        
        m2<-rpart(train_X~train_Z1,data = my_data,method  = "anova",control = list(xval=5,cp=0,minsplit=2,minbucket=minbucket_candidate[qqq]))
        # rpart.plot(m2)
        
        
        pred_value_train_initial<-predict(m2, newdata = my_data)
        
        pred_value_test_initial<-predict(m2, newdata = test_data)
        
        # res_period1<-solnp(pars = mean(train_X),fun = price_1,LB=0,control = list(trace=0))
        # res_period1<-res_period1$pars
        
        price_2_test<-function(z,fff=numeric(0),starting_level=numeric(0)){
          aaa<-which(pred_value_train==pred_value_test[fff]) 
          temp_trainY<-train_Y[aaa]
          return(mean(h_cost[2]*pmax(0,pmax(z,starting_level)-temp_trainY)+b_cost[2]*pmax(0,temp_trainY-pmax(z,starting_level))))
        }
        
        res<-0
        for(rrr in 1:length(test_X)){
          my_qqq<-which(pred_value_test_initial[rrr]==pred_value_train_initial)
          weight_vector<-rep(0,length(train_X))
          weight_vector[my_qqq]<-1
          res_period1<-solnp(pars = mean(train_X),fun = price_1,LB=0,control = list(trace=0),weight_vector=weight_vector)
          res_period1<-res_period1$pars
          z_level<-solnp(pars = mean(train_Y),fun = price_2_test,control = list(trace=0),fff=rrr,starting_level=res_period1-test_X[rrr],LB=0)$pars
          z_level<-max(z_level,res_period1-test_X[rrr])
          res<-res+h_cost[1]*max(0,res_period1-test_X[rrr])+b_cost[1]*max(0,test_X[rrr]-res_period1)+h_cost[2]*max(0,z_level-test_Y[rrr])+b_cost[2]*max(0,test_Y[rrr]-z_level)
        }
        CV_array[qqq,fold]<-res
        print(c(qqq,fold))
      }
    }
    CV_seq<-apply(CV_array,MARGIN = 1,sum)
    CVtime_cart<-as.numeric((proc.time()-ptm)[3])
    CV_seq_cart<-CV_seq
    ind_selected<-which.min(CV_seq)[1]
    minbucket<-minbucket_candidate[ind_selected] 
    ind_selected_cart<-ind_selected
    
    
    train_X<-D_mat[,1]
    train_Y<-D_mat[,2]
    train_Z1<-X_mat[,1]
    train_Z2<-X_mat[,2]
    
    
    test_X<-D_mat_test[,1]
    test_Y<-D_mat_test[,2]
    test_Z1<-X_mat_test[,1]
    test_Z2<-X_mat_test[,2]
    
    my_data<-as.data.frame(cbind(train_Z1,train_X))
    test_data<-as.data.frame(cbind(test_Z1,test_X))
    names(test_data)<-c("train_Z1","train_X")
    
    
    my_data1<-as.data.frame(cbind(train_Z2,train_X,train_Y))
    test_data1<-as.data.frame(cbind(test_Z2,test_X,test_Y))
    names(test_data1)<-c("train_Z2","train_X","train_Y")
    
    
    m1<-rpart(train_Y~train_Z2+train_X,data = my_data1,method  = "anova",control = list(xval=5,cp=0,minsplit=2,minbucket=minbucket))
    # rpart.plot(m1)
    
    
    pred_value_train<-predict(m1, newdata = my_data1)
    
    pred_value_test<-predict(m1, newdata = test_data1)
    
    
    price_2<-function(z,fff=numeric(0),starting_level=numeric(0)){
      aaa<-which(pred_value_train==pred_value_train[fff]) 
      temp_trainY<-train_Y[aaa]
      return(mean(h_cost[2]*pmax(0,pmax(z,starting_level)-temp_trainY)+b_cost[2]*pmax(0,temp_trainY-pmax(z,starting_level))))
    }
    
    price_1<-function(z0,weight_vector=numeric(0)){
      next_cost<-rep(0,length(train_X))
      for(www in 1:length(train_X)){
        if(weight_vector[www]>0){
          temp_www<-solnp(pars = mean(train_Y),fun = price_2,control = list(trace=0),fff=www,starting_level=z0-train_X[www],LB=0)$values
          next_cost[www]<-temp_www[length(temp_www)]
        }
      }
      return(sum(weight_vector*(h_cost[1]*pmax(z0-train_X,0)+b_cost[1]*pmax(train_X-z0,0)+next_cost))) 
    }
    
    m2<-rpart(train_X~train_Z1,data = my_data,method  = "anova",control = list(xval=5,cp=0,minsplit=2,minbucket=minbucket))
    
    pred_value_train_initial<-predict(m2, newdata = my_data)
    
    pred_value_test_initial<-predict(m2, newdata = test_data)
    
    price_2_test<-function(z,fff=numeric(0),starting_level=numeric(0)){
      aaa<-which(pred_value_train==pred_value_test[fff]) 
      temp_trainY<-train_Y[aaa]
      return(mean(h_cost[2]*pmax(0,pmax(z,starting_level)-temp_trainY)+b_cost[2]*pmax(0,temp_trainY-pmax(z,starting_level))))
    }
    
    res<-rep(0,length(test_X))
    time_cart<-rep(0,length(test_X))
    for(rrr in 1:length(test_X)){
      ptm<-proc.time()
      my_qqq<-which(pred_value_test_initial[rrr]==pred_value_train_initial)
      weight_vector<-rep(0,length(train_X))
      weight_vector[my_qqq]<-1
      res_period1<-solnp(pars = mean(train_X),fun = price_1,LB=0,control = list(trace=0),weight_vector=weight_vector)
      res_period1<-res_period1$pars
      z_level<-solnp(pars = mean(train_Y),fun = price_2_test,control = list(trace=0),fff=rrr,starting_level=res_period1-test_X[rrr],LB=0)$pars
      z_level<-max(z_level,res_period1-test_X[rrr])
      res[rrr]<-res[rrr]+h_cost[1]*max(0,res_period1-test_X[rrr])+b_cost[1]*max(0,test_X[rrr]-res_period1)+h_cost[2]*max(0,z_level-test_Y[rrr])+b_cost[2]*max(0,test_Y[rrr]-z_level)
      time_cart[rrr]<-as.numeric((proc.time()-ptm)[3])
    }
    cost_cart<-res
    
    
    
    
    
    ########################################################################################################################################################################################################################################
    ########################################################################################################################################################################################################################################
    ptm<-proc.time()
    nodesize_candidate<-rev(unique(round(seq(2,n/3,length.out=5)) )) 
    CV_mat<-matrix(0,5,length(nodesize_candidate))
    my_split<-1:n
    my_split<-matrix(my_split,5,n/5,byrow = TRUE)
    for(fold in 1:5){
      test_X<-D_mat[c(my_split[fold,]),1]
      test_Y<-D_mat[c(my_split[fold,]),2]
      test_Z1<-X_mat[c(my_split[fold,]),1]
      test_Z2<-X_mat[c(my_split[fold,]),2]
      
      train_X<-D_mat[c(my_split[-fold,]),1]
      train_Y<-D_mat[c(my_split[-fold,]),2]
      train_Z1<-X_mat[c(my_split[-fold,]),1]
      train_Z2<-X_mat[c(my_split[-fold,]),2]
      
      my_data<-as.data.frame(cbind(train_Z1,train_X))
      test_data<-as.data.frame(cbind(test_Z1,test_X))
      names(test_data)<-c("train_Z1","train_X")
      
      
      my_data1<-as.data.frame(cbind(train_Z2,train_X,train_Y))
      test_data1<-as.data.frame(cbind(test_Z2,test_X,test_Y))
      names(test_data1)<-c("train_Z2","train_X","train_Y")
      for(ppp in 1:length(nodesize_candidate)){
        aaa<-randomForest(train_Y~train_Z2+train_X,mtry=1,nodesize=nodesize_candidate[ppp],data=my_data1)
        ccc<-predict(aaa,my_data1,predict.all = TRUE)
        ddd_train<-ccc$individual
        
        ccc<-predict(aaa,test_data1,predict.all = TRUE)
        ddd_test<-ccc$individual
        
        temp_w_mat<-matrix(0,length(train_X),length(train_X))
        temp_weight_matrix<-array(0,dim = c(length(train_X),500,length(train_X)))
        for(hhh in 1:length(train_X)){
          for(kkk in 1:500){
            temp_weight_matrix[which(ddd_train[hhh,kkk]==ddd_train[,kkk]),kkk,hhh]<-1/length(which(ddd_train[hhh,kkk]==ddd_train[,kkk])) 
          }
          temp_w_mat[,hhh]<-apply(temp_weight_matrix[,,hhh],MARGIN = 1,mean) 
        }
        
        
        price_2<-function(z,fff=numeric(0),starting_level=numeric(0)){
          temp_w<-temp_w_mat[,fff]
          return(sum((h_cost[2]*pmax(0,pmax(z,starting_level)-train_Y)+b_cost[2]*pmax(0,train_Y-pmax(z,starting_level)))*temp_w))
        }
        
        price_1<-function(z0,weight_vector=numeric(0)){
          next_cost<-rep(0,length(train_X))
          for(www in 1:length(train_X)){
            if(weight_vector[www]>0){
              temp_www<-solnp(pars = mean(train_Y),fun = price_2,control = list(trace=0),fff=www,starting_level=z0-train_X[www],LB=0)$values
              next_cost[www]<-temp_www[length(temp_www)]
            }
          }
          return(sum(weight_vector*(h_cost[1]*pmax(z0-train_X,0)+b_cost[1]*pmax(train_X-z0,0)+next_cost))) 
        }
        
        
        aaa1<-randomForest(train_X~train_Z1,mtry=1,nodesize=nodesize_candidate[ppp],data=my_data)
        ccc1<-predict(aaa1,my_data,predict.all = TRUE)
        ddd_train1<-ccc1$individual
        
        ccc1<-predict(aaa1,test_data,predict.all = TRUE)
        ddd_test1<-ccc1$individual
        
        old_temp_w_mat<-temp_w_mat
        
        temp_w_mat<-matrix(0,length(train_X),length(test_X))
        temp_weight_matrix<-array(0,dim = c(length(train_X),500,length(test_X)))
        for(hhh in 1:length(test_X)){
          for(kkk in 1:500){
            temp_weight_matrix[which(ddd_test[hhh,kkk]==ddd_train[,kkk]),kkk,hhh]<-1/length(which(ddd_test[hhh,kkk]==ddd_train[,kkk])) 
          }
          temp_w_mat[,hhh]<-apply(temp_weight_matrix[,,hhh],MARGIN = 1,mean) 
        }
        
        temp_w_mat1<-temp_w_mat
        temp_w_mat<-old_temp_w_mat
        
        price_2_test<-function(z,fff=numeric(0),starting_level=numeric(0)){
          temp_w<-temp_w_mat1[,fff]
          return(sum((h_cost[2]*pmax(0,pmax(z,starting_level)-train_Y)+b_cost[2]*pmax(0,train_Y-pmax(z,starting_level)))*temp_w))
        }
        
        res<-0
        for(rrr in 1:length(test_X)){
          my_qqq<-matrix(0,500,length(train_X))
          for(ggg in 1:500){
            my_qqq[ggg,which(ddd_test1[rrr,ggg]==ddd_train1[,ggg])]<-1/length(which(ddd_test1[rrr,ggg]==ddd_train1[,ggg]))
          }
          weight_vector<-apply(my_qqq,MARGIN = 2,mean)
          res_period1<-solnp(pars = mean(train_X),fun = price_1,LB=0,control = list(trace=0),weight_vector=weight_vector)
          res_period1<-res_period1$pars
          z_level<-solnp(pars = mean(train_Y),fun = price_2_test,control = list(trace=0),fff=rrr,starting_level=res_period1-test_X[rrr],LB=0)$pars
          z_level<-max(z_level,res_period1-test_X[rrr])
          res<-res+h_cost[1]*max(0,res_period1-test_X[rrr])+b_cost[1]*max(0,test_X[rrr]-res_period1)+h_cost[2]*max(0,z_level-test_Y[rrr])+b_cost[2]*max(0,test_Y[rrr]-z_level)
        }
        CV_mat[fold,ppp]<-res
        print(c(fold,ppp))
      }
    }
    CV_seq<-apply(CV_mat,MARGIN = 2,sum)
    CVtime_rf<-as.numeric((proc.time()-ptm)[3])
    
    
    CV_seq_rf<-CV_seq
    ind_selected<-which(CV_seq == min(CV_seq), arr.ind = TRUE)[1]
    nodesize<-nodesize_candidate[ind_selected] 
    ind_selected_rf<-ind_selected
    
    
    train_X<-D_mat[,1]
    train_Y<-D_mat[,2]
    train_Z1<-X_mat[,1]
    train_Z2<-X_mat[,2]
    
    
    test_X<-D_mat_test[,1]
    test_Y<-D_mat_test[,2]
    test_Z1<-X_mat_test[,1]
    test_Z2<-X_mat_test[,2]
    
    my_data<-as.data.frame(cbind(train_Z1,train_X))
    test_data<-as.data.frame(cbind(test_Z1,test_X))
    names(test_data)<-c("train_Z1","train_X")
    
    
    my_data1<-as.data.frame(cbind(train_Z2,train_X,train_Y))
    test_data1<-as.data.frame(cbind(test_Z2,test_X,test_Y))
    names(test_data1)<-c("train_Z2","train_X","train_Y")
    
    
    aaa<-randomForest(train_Y~train_Z2+train_X,mtry=1,nodesize=nodesize,data=my_data1)
    ccc<-predict(aaa,my_data1,predict.all = TRUE)
    ddd_train<-ccc$individual
    
    ccc<-predict(aaa,test_data1,predict.all = TRUE)
    ddd_test<-ccc$individual
    
    temp_w_mat<-matrix(0,length(train_X),length(train_X))
    temp_weight_matrix<-array(0,dim = c(length(train_X),500,length(train_X)))
    for(hhh in 1:length(train_X)){
      for(kkk in 1:500){
        temp_weight_matrix[which(ddd_train[hhh,kkk]==ddd_train[,kkk]),kkk,hhh]<-1/length(which(ddd_train[hhh,kkk]==ddd_train[,kkk])) 
      }
      temp_w_mat[,hhh]<-apply(temp_weight_matrix[,,hhh],MARGIN = 1,mean) 
    }
    
    
    price_2<-function(z,fff=numeric(0),starting_level=numeric(0)){
      temp_w<-temp_w_mat[,fff]
      return(sum((h_cost[2]*pmax(0,pmax(z,starting_level)-train_Y)+b_cost[2]*pmax(0,train_Y-pmax(z,starting_level)))*temp_w))
    }
    
    price_1<-function(z0,weight_vector=numeric(0)){
      next_cost<-rep(0,length(train_X))
      for(www in 1:length(train_X)){
        if(weight_vector[www]>0){
          temp_www<-solnp(pars = mean(train_Y),fun = price_2,control = list(trace=0),fff=www,starting_level=z0-train_X[www],LB=0)$values
          next_cost[www]<-temp_www[length(temp_www)]
        }
      }
      return(sum(weight_vector*(h_cost[1]*pmax(z0-train_X,0)+b_cost[1]*pmax(train_X-z0,0)+next_cost))) 
    }
    
    
    aaa1<-randomForest(train_X~train_Z1,mtry=1,nodesize=nodesize,data=my_data)
    ccc1<-predict(aaa1,my_data,predict.all = TRUE)
    ddd_train1<-ccc1$individual
    
    ccc1<-predict(aaa1,test_data,predict.all = TRUE)
    ddd_test1<-ccc1$individual
    
    old_temp_w_mat<-temp_w_mat
    
    temp_w_mat<-matrix(0,length(train_X),length(test_X))
    temp_weight_matrix<-array(0,dim = c(length(train_X),500,length(test_X)))
    for(hhh in 1:length(test_X)){
      for(kkk in 1:500){
        temp_weight_matrix[which(ddd_test[hhh,kkk]==ddd_train[,kkk]),kkk,hhh]<-1/length(which(ddd_test[hhh,kkk]==ddd_train[,kkk])) 
      }
      temp_w_mat[,hhh]<-apply(temp_weight_matrix[,,hhh],MARGIN = 1,mean) 
    }
    
    temp_w_mat1<-temp_w_mat
    temp_w_mat<-old_temp_w_mat
    
    price_2_test<-function(z,fff=numeric(0),starting_level=numeric(0)){
      temp_w<-temp_w_mat1[,fff]
      return(sum((h_cost[2]*pmax(0,pmax(z,starting_level)-train_Y)+b_cost[2]*pmax(0,train_Y-pmax(z,starting_level)))*temp_w))
    }
    
    res<-rep(0,length(test_X))
    time_rf<-rep(0,length(test_X))
    for(rrr in 1:length(test_X)){
      ptm<-proc.time()
      my_qqq<-matrix(0,500,length(train_X))
      for(ggg in 1:500){
        my_qqq[ggg,which(ddd_test1[rrr,ggg]==ddd_train1[,ggg])]<-1/length(which(ddd_test1[rrr,ggg]==ddd_train1[,ggg]))
      }
      weight_vector<-apply(my_qqq,MARGIN = 2,mean)
      res_period1<-solnp(pars = mean(train_X),fun = price_1,LB=0,control = list(trace=0),weight_vector=weight_vector)
      res_period1<-res_period1$pars
      z_level<-solnp(pars = mean(train_Y),fun = price_2_test,control = list(trace=0),fff=rrr,starting_level=res_period1-test_X[rrr],LB=0)$pars
      z_level<-max(z_level,res_period1-test_X[rrr])
      res[rrr]<-res[rrr]+h_cost[1]*max(0,res_period1-test_X[rrr])+b_cost[1]*max(0,test_X[rrr]-res_period1)+h_cost[2]*max(0,z_level-test_Y[rrr])+b_cost[2]*max(0,test_Y[rrr]-z_level)
      time_rf[rrr]<-as.numeric((proc.time()-ptm)[3])
    }
    cost_rf<-res
    
    ########################################################################################################################################################################################################################################
    ########################################################################################################################################################################################################################################
    ###calculate the cost of the optimal policy
    load("cov_2_true1.rdata")
    
    knots_1 <- res[[1]][[1]]
    knots_2 <- res[[1]][[2]]
    
    
    pars1 <- res[[2]][[1]]
    pars2 <- res[[2]][[2]]
    
    basis_test_1<-bSpline(x=X_mat_test[,1],knots = knots_1,degree = 3,Boundary.knots = c(1,5),intercept = TRUE )
    test_decision_1<-c(basis_test_1%*%pars1)
    
    temp_x<-c(X_mat_test[,2]*pars2[1]+D_mat_test[,1])
    temp_base<-bSpline(x=temp_x,knots = knots_2,degree = 3,Boundary.knots = c(0.9,13.5),intercept = TRUE )
    test_decision_2<-c(temp_base%*%pars2[-1])
    
    temp_res<-rep(0,10)
    for(ooo in 1:10){
      temp_res[ooo]<-temp_res[ooo]+h_cost[1]*max(0,test_decision_1[ooo]-D_mat_test[ooo,1])+b_cost[1]*max(0,D_mat_test[ooo,1]-test_decision_1[ooo])
      temp_res[ooo]<-temp_res[ooo]+h_cost[2]*max(0,max(test_decision_1[ooo]-D_mat_test[ooo,1],test_decision_2[ooo])-D_mat_test[ooo,2])+b_cost[2]*max(0,D_mat_test[ooo,2]-max(test_decision_1[ooo]-D_mat_test[ooo,1],test_decision_2[ooo]))
    }
    cost_oracle <- temp_res
    
    return(list(cost_ZX,cost_RKHS,cost_knn,cost_cart,cost_rf,time_ZX,time_RKHS,time_knn,time_cart,time_rf,CVtime_ZX,CVtime_RKHS,CVtime_knn,CVtime_cart,CVtime_rf,m_n,ind_selected_RKHS,k_num,ind_selected_cart,ind_selected_rf,CV_seq_ZX,CV_mat_RKHS,CV_seq_knn,CV_seq_cart,CV_seq_rf,D_mat,X_mat,D_mat_test,X_mat_test,cost_oracle))
  }
  res<-try(my_simulation_procedure(iter),silent = TRUE)
  return(res)
}


clnum<-20#detectCores()/2
cl <- makeCluster(getOption("cl.cores", clnum))
res<-parLapply(cl, 1:100,  simulation)
stopCluster(cl)

save(res,file = "T equals 2 with covariates.rdata")
