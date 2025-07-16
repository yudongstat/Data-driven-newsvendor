rm(list=ls())
library(parallel)
cores_called <- 50
simulation<-function(iter){
  my_simulation_procedure<-function(iter){
    library(splines2)
    library(pracma)
    library(Rsolnp)
    library(rpart)
    library(rpart.plot)
    library(randomForest)
    load("ye_cov_n200.rdata")
    full_D_mat <- res[[iter]][[26]]
    full_X_mat <- res[[iter]][[27]]
    
    b_cost<-c(5,1,1)
    h_cost<-c(1,5,1)
    
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
        
        nodesize_candidate <- rev(unique(round(seq(2, max(cut_off, 10)/3, length.out = 5))))
        
        if(cut_off == 0){
          # at the first selling season, we randomly choose the tuning parameters
          nodesize <- sample(x = nodesize_candidate, size = 1)
          
        }else{
          # update the tuning parameters based on cross validation
          # RF
          CV_mat<-matrix(0,5,length(nodesize_candidate))
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
          ind_selected<-which(CV_seq == min(CV_seq), arr.ind = TRUE)[1]
          nodesize<-nodesize_candidate[ind_selected] 
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
        
        my_data<-as.data.frame(cbind(train_Z1,train_X))
        
        
        my_data1<-as.data.frame(cbind(train_Z2,train_X,train_Y))
        
        
        aaa<-randomForest(train_Y~train_Z2+train_X,mtry=1,nodesize=nodesize,data=my_data1)
        ccc<-predict(aaa,my_data1,predict.all = TRUE)
        ddd_train<-ccc$individual
        
        ccc<-predict(aaa,data.frame(train_Z2 = test_Z2,
                                    train_X = test_X,
                                    train_Y = test_Y),
                     predict.all = TRUE)
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
        
        ccc1<-predict(aaa1, data.frame(train_Z1 = test_Z1,
                                       train_X = test_X),
                      predict.all = TRUE)
        ddd_test1<-ccc1$individual
        
        old_temp_w_mat<-temp_w_mat
        
        temp_weight_matrix <- matrix(0, length(train_X), 500)
        for(kkk in 1:500){
          temp_weight_matrix[which(ddd_test[1, kkk]==ddd_train[,kkk]),kkk] <- 1/length(which(ddd_test[1, kkk]==ddd_train[,kkk])) 
        }
        temp_w_mat <- apply(temp_weight_matrix, MARGIN = 1, mean) 
        
        temp_w_mat1<-temp_w_mat
        temp_w_mat<-old_temp_w_mat
        
        price_2_test<-function(z, starting_level=numeric(0)){
          temp_w<-temp_w_mat1
          return(sum((h_cost[2]*pmax(0,pmax(z,starting_level)-train_Y)+b_cost[2]*pmax(0,train_Y-pmax(z,starting_level)))*temp_w))
        }
        
        my_qqq<-matrix(0,500,length(train_X))
        for(ggg in 1:500){
          my_qqq[ggg,which(ddd_test1[1, ggg]==ddd_train1[,ggg])]<-1/length(which(ddd_test1[1, ggg]==ddd_train1[,ggg]))
        }
        weight_vector<-apply(my_qqq,MARGIN = 2,mean)
        res_period1<-solnp(pars = mean(train_X),fun = price_1,LB=0,control = list(trace=0),weight_vector=weight_vector)
        test_decision_1 <- res_period1$pars
        test_decision_2 <-solnp(pars = mean(train_Y),fun = price_2_test,control = list(trace=0), starting_level=test_decision_1-test_X,LB=0)$pars
        
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

save(res,file = "online_RF.rdata")