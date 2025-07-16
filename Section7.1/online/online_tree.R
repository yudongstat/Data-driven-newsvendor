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
        
        minbucket_candidate <- rev(unique(round(seq(2, max(cut_off, 10)/3, length.out = 5))))

        if(cut_off == 0){
          # at the first selling season, we randomly choose the tuning parameters
          minbucket <- sample(x = minbucket_candidate, size = 1)
          
        }else{
          # update the tuning parameters based on cross validation
          # CART
          CV_array<-matrix(0,length(minbucket_candidate),5)
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
          ind_selected<-which.min(CV_seq)[1]
          minbucket<-minbucket_candidate[ind_selected] 
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
        
        
        m1<-rpart(train_Y~train_Z2+train_X,data = my_data1,method  = "anova",control = list(xval=5,cp=0,minsplit=2,minbucket=minbucket))
        # rpart.plot(m1)
        
        
        pred_value_train<-predict(m1, newdata = my_data1)
        
        pred_value_test<-predict(m1, newdata = data.frame(train_Z2 = test_Z2, 
                                                          train_X = test_X,
                                                          train_Y = test_Y))
        
        
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
        
        pred_value_test_initial<-predict(m2, newdata = data.frame(train_Z1 = test_Z1,
                                                                  train_X = test_X))
        
        price_2_test<-function(z, starting_level=numeric(0)){
          aaa<-which(pred_value_train==pred_value_test) 
          temp_trainY<-train_Y[aaa]
          return(mean(h_cost[2]*pmax(0,pmax(z,starting_level)-temp_trainY)+b_cost[2]*pmax(0,temp_trainY-pmax(z,starting_level))))
        }
        
        
        my_qqq<-which(pred_value_test_initial==pred_value_train_initial)
        weight_vector<-rep(0,length(train_X))
        weight_vector[my_qqq]<-1
        res_period1<-solnp(pars = mean(train_X),fun = price_1,LB=0,control = list(trace=0),weight_vector=weight_vector)
        test_decision_1 <- res_period1$pars
        test_decision_2 <-solnp(pars = mean(train_Y),fun = price_2_test,control = list(trace=0), starting_level= test_decision_1 - test_X,LB=0)$pars
        
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

save(res,file = "online_tree.rdata")