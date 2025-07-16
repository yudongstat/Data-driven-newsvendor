rm(list=ls())
library(parallel)
cores_called <- 50
simulation<-function(iter){
  my_simulation_procedure<-function(iter){
    library(splines2)
    library(pracma)
    library(Rsolnp)
    load("ye_cov_n200.rdata")
    full_D_mat <- res[[iter]][[26]]
    full_X_mat <- res[[iter]][[27]]
    
    b_cost<-c(5,1,1)
    h_cost<-c(1,5,1)
    
    euclidean<-function(x,y){
      res<-rep(0,nrow(y))
      for(i in 1:nrow(y)){
        res[i]<-sqrt(sum((x-y[i,])^2)) 
      }
      return(res)
    }
    
    tuning_length <- 10
    
    window_length <- 100
    my_online_res <- matrix(NA, window_length, 3)
    
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
        
        neighbor_candidate <- unique(round(seq(2, max(cut_off, 10)/3, length.out = 5)))
        
        if(cut_off == 0){
          # at the first selling season, we randomly choose the tuning parameters
          k_num <- sample(x = neighbor_candidate, size = 1)
          
        }else{
          # update the tuning parameters based on cross validation
          # knn
          CV_mat<-matrix(0,5,length(neighbor_candidate))
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
          k_num<-neighbor_candidate[(which.min(CV_seq))[1]]
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
        
        
        test_X<-D_test[1]
        test_Y<-D_test[2]
        test_Z1<-X_test[1]
        test_Z2<-X_test[2]
        
        dist_mat_train<-matrix(0,length(train_X),length(train_X))
        for(xxx in 1:length(train_X)){
          dist_mat_train[xxx,]<-order(euclidean(x=c(train_Z2[xxx],train_X[xxx]),y=cbind(train_Z2,train_X)))
        }
        
        dist_mat_test <- order(euclidean(x=c(test_Z2,test_X),y=cbind(train_Z2,train_X)))
        
        dist_mat_test1<- order(abs(test_Z1 - train_Z1))
        
        price_2<-function(z,x=numeric(0),starting_level=numeric(0)){
          temp_ind<-dist_mat_train[x,]
          temp_trainY<-train_Y[temp_ind[1:k_num]]
          return(mean(h_cost[2]*pmax(0,pmax(starting_level,z)-temp_trainY)+b_cost[2]*pmax(0,temp_trainY-pmax(starting_level,z))))
        }
        
        price_1<-function(z0){
          temp_ind<-dist_mat_test1
          temp_ind<-temp_ind[1:k_num]
          next_cost<-rep(0,k_num)
          for(www in 1:k_num){
            temp_www<-solnp(pars = mean(train_Y),fun = price_2,control = list(trace=0),x=temp_ind[www],starting_level=z0-train_X[temp_ind[www]],LB=0)$values
            next_cost[www]<-temp_www[length(temp_www)]
          }
          temp_trainX<-train_X[temp_ind[1:k_num]]
          return(sum(h_cost[1]*pmax(z0-temp_trainX,0)+b_cost[1]*pmax(temp_trainX-z0,0)+next_cost)) 
        }
        
        price_2_test<-function(z,starting_level=numeric(0)){
          temp_ind<-dist_mat_test
          temp_trainY<-train_Y[temp_ind[1:k_num]]
          return(mean(h_cost[2]*pmax(0,pmax(starting_level,z)-temp_trainY)+b_cost[2]*pmax(0,temp_trainY-pmax(starting_level,z))))
        }
        
        res_period1<-solnp(pars = mean(train_X),fun = price_1,LB=0,control = list(trace=0))
        test_decision_1 <- res_period1$pars
        test_decision_2<-solnp(pars = mean(train_Y),fun = price_2_test,control = list(trace=0), starting_level=test_decision_1-test_X,LB=0)$pars
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

save(res,file = "online_knn.rdata")