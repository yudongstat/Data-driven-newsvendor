rm(list = ls())
library(splines2)
library(pracma)
library(Rsolnp)
chosen_index <- 2
load("bakery_data.rdata")

T <- 6
h_cost <- rep(1, T)
b_cost <- rep(1, T)
n <- nrow(my_data)/T

full_D_mat <- matrix(my_data[[products_all[indexing_order[chosen_index]]]],
                     n, T, byrow = TRUE)

cost_list_vanilla <- rep(NA, n)
time_vanilla <- rep(NA, n)

for(cut_off in 0:(n - 1)){
  start_time <- proc.time()
  
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
    
    train_D_mat <- D_mat
    test_D_mat <- full_D_mat[cut_off + 1, ]
    
    price_6<-function(alpha){
      temp_level<-alpha
      return(sum(h_cost[6]*pmax(0,temp_level-train_D_mat[,6])+b_cost[6]*pmax(0,train_D_mat[,6]-temp_level)) )
    }
    
    price_6_empirical<-function(temp_level){
      return(sum(h_cost[6]*pmax(0,temp_level-train_D_mat[,6])+b_cost[6]*pmax(0,train_D_mat[,6]-temp_level)) )
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
    
    res_5<-solnp(pars = mean(train_D_mat[,5]),fun = price_5,LB=0,control = list(trace=0))
    best_decision_5<-res_5$pars
    
    
    price_4<-function(alpha){
      temp_level<-alpha
      return(sum(h_cost[4]*pmax(0,temp_level-train_D_mat[,4])+b_cost[4]*pmax(0,train_D_mat[,4]-temp_level))+price_5_empirical(pmax(best_decision_5,temp_level-train_D_mat[,4])) )
    }
    
    price_4_empirical<-function(temp_level){
      return(sum(h_cost[4]*pmax(0,temp_level-train_D_mat[,4])+b_cost[4]*pmax(0,train_D_mat[,4]-temp_level))+price_5_empirical(pmax(best_decision_5,temp_level-train_D_mat[,4])) )
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
    
    res_3<-solnp(pars = mean(train_D_mat[,3]),fun = price_3,LB=0,control = list(trace=0))
    best_decision_3<-res_3$pars
    
    
    
    
    price_2<-function(alpha){
      temp_level<-alpha
      return(sum(h_cost[2]*pmax(0,temp_level-train_D_mat[,2])+b_cost[2]*pmax(0,train_D_mat[,2]-temp_level))+price_3_empirical(pmax(temp_level-train_D_mat[,2],best_decision_3)))
    }
    
    price_2_empirical<-function(temp_level){
      return(sum(h_cost[2]*pmax(0,temp_level-train_D_mat[,2])+b_cost[2]*pmax(0,train_D_mat[,2]-temp_level))+price_3_empirical(pmax(temp_level-train_D_mat[,2],best_decision_3)))
    }
    
    res_2<-solnp(pars = mean(train_D_mat[,2]),fun = price_2,LB=0,control = list(trace=0))
    best_decision_2<-res_2$pars
    
    
    price_1<-function(temp_level){
      return(sum(h_cost[1]*pmax(0,temp_level-train_D_mat[,1])+b_cost[1]*pmax(0,train_D_mat[,1]-temp_level))+price_2_empirical(pmax(temp_level-train_D_mat[,1],best_decision_2)))
    }
    price_1_empirical_var_est<-function(temp_level){
      return((h_cost[1]*pmax(0,temp_level-train_D_mat[,1])+b_cost[1]*pmax(0,train_D_mat[,1]-temp_level))+price_2_empirical_var_est(pmax(temp_level-train_D_mat[,1],best_decision_2)))
    }
    res_1<-optimize(f=price_1,interval = c(50,500))
    
    best_decision_1<-res_1$minimum
    
    test_decision_our <- c(best_decision_1,
                           best_decision_2,
                           best_decision_3,
                           best_decision_4,
                           best_decision_5,
                           best_decision_6)
  }
  
  time_vanilla[cut_off + 1] <- (proc.time() - start_time)["elapsed"]
  
  
  test_cost_our <- 0
  temp_level_our <- test_decision_our[1]
  for(qqqq in 1:T){
    test_cost_our<-test_cost_our+(h_cost[qqqq]*pmax(0,temp_level_our-test_D_mat[qqqq])+b_cost[qqqq]*pmax(0,test_D_mat[qqqq]-temp_level_our))
    
    if(qqqq<T){
      temp_level_our<-pmax(temp_level_our-test_D_mat[qqqq],test_decision_our[qqqq+1])
    }
  }
  cost_list_vanilla[cut_off + 1] <- test_cost_our
  
  print(c(cut_off, test_cost_our))
}


cost_vanilla <- cost_list_vanilla

save(cost_vanilla, time_vanilla, file = "vanilla.rdata")
