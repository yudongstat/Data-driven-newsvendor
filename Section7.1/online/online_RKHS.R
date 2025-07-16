rm(list=ls())
library(parallel)
cores_called <- 50
simulation<-function(iter){
  my_simulation_procedure<-function(iter){
    library(splines2)
    library(pracma)
    library(Rsolnp)
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
    
    full_D_mat <- D_mat
    full_X_mat <- X_mat
    
    b_cost<-c(5,1,1)
    h_cost<-c(1,5,1)
    
    gaussian_kernel<-function(x1,x2,sigma=1){
      return(exp(-(x1-x2)^2/sigma^2))
    }
    gaussian_kernel_high<-function(x1,x2,sigma=1){
      return(exp(-sum((x1-x2)^2)/sigma^2))
    }
    
    tuning_length <- 10
    
    window_length <- 100
    my_online_res <- matrix(NA, window_length, 3)
    
    for(cut_off in 0:(window_length - 1)){
      n <- cut_off
      if(cut_off > 0){
        D_mat <- full_D_mat[1:cut_off, ]
        X_mat <- full_X_mat[1:cut_off, ]
      }
      
      X_test <- full_X_mat[cut_off + 1, ]
      D_test <- full_D_mat[cut_off + 1, ]
      
      ptm <- proc.time()
      # we update the tuning parameters every 10 selling seasons
      if((cut_off + 1) %% tuning_length == 1){

        sigma_candidate <- seq(1,20,1)
        lambda_candidate <- 5^(-c(1:10))
        
        if(cut_off == 0){
          # at the first selling season, we randomly choose the tuning parameters
          sigma <- sample(x = sigma_candidate, size = 1)
          mylambda <- sample(x = lambda_candidate, size = 1)
          
        }else{
          # update the tuning parameters based on cross validation
          
          # RKHS
          CV_array<-array(0,dim = c(length(sigma_candidate),length(lambda_candidate),5))
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
          ind_selected<-which(CV_mat == min(CV_mat), arr.ind = TRUE)[1,]
          sigma<-sigma_candidate[ind_selected[1]] 
          mylambda<-lambda_candidate[ind_selected[2]]
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
        
        k_mat_test<-gaussian_kernel(x1=X_test[1],x2=train_Z1,sigma = sigma)
        test_decision_1 <- sum(k_mat_test * opt_res$pars[1:(n)])
        
        k_mat_test_high <- rep(0, n)
        for(j in 1:(n)){
          k_mat_test_high[j]<-gaussian_kernel_high(x1=c(X_test[2], D_test[1]),x2=c(train_Z2[j],train_X[j]),sigma = sigma)
        }
        test_decision_2<- sum(k_mat_test_high * opt_res$pars[-c(1:(n))])
      }
      
      current_level <- test_decision_1
      temp_res <- h_cost[1]*max(0,current_level - D_test[1]) + b_cost[1]*max(0,D_test[1]-current_level)
      current_level <- max(current_level - D_test[1], test_decision_2)
      temp_res <- temp_res + h_cost[2]*max(0, current_level - D_test[2])+b_cost[2]*max(0,D_test[2] - current_level)
      cost <- temp_res
      time <- as.numeric((proc.time()-ptm)[3])
      my_online_res[cut_off + 1, ] <- c(cost, time, CVtime)
      print(c(cut_off + 1, cost, time, CVtime))
    }
    
    return(my_online_res)
  }
  res<-try(my_simulation_procedure(iter),silent = TRUE)
  return(res)
}


clnum <- cores_called#detectCores()/2
cl <- makeCluster(getOption("cl.cores", clnum))
res<-parLapply(cl, 1:100,  simulation)
stopCluster(cl)

save(res,file = "online_RKHS.rdata")
