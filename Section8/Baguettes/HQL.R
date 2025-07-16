rm(list = ls())
library(splines2)
library(pracma)
library(Rsolnp)
chosen_index <- 1
load("bakery_data.rdata")

T <- 6
H <- T
h_cost <- rep(1, T)
b_cost <- rep(1, T)
n <- nrow(my_data)/T

full_D_mat <- matrix(my_data[[products_all[indexing_order[chosen_index]]]],
                     n, T, byrow = TRUE)


upper_bound <- 710
action_space <- seq(0, upper_bound, 0.2)
best_base_level <- rep(200, T)

train_D_mat <- full_D_mat

# initialization
Q_table <- matrix(0, H, length(action_space))
V_table <- matrix(0, H, length(action_space))

feasible_actions <- as.list(rep(NA, H + 1))
feasible_actions[[1]] <- 1:length(action_space)
feasible_actions[[2]] <- 1:length(action_space)
feasible_actions[[3]] <- 1:length(action_space)
feasible_actions[[4]] <- 1:length(action_space)
feasible_actions[[5]] <- 1:length(action_space)
feasible_actions[[6]] <- 1:length(action_space)
feasible_actions[[7]] <- 1:length(action_space)



test_cost_HQL <- 0
temp_level_HQL <- best_base_level[1]
for(qqqq in 1:T){
  test_cost_HQL <- test_cost_HQL + (h_cost[qqqq]*pmax(0, temp_level_HQL - train_D_mat[1, qqqq]) + b_cost[qqqq]*pmax(0, train_D_mat[1, qqqq] - temp_level_HQL))
  if(qqqq<T){
    temp_level_HQL <- pmax(temp_level_HQL - train_D_mat[1, qqqq], best_base_level[qqqq + 1])
  }
}

cost_HQL_list <- rep(NA, n)
cost_HQL_list[1] <- test_cost_HQL

time_HQL_list <- rep(NA, n)

for(epoch in 1:(n-1)){
  start_time <- proc.time()
  alpha <- (H + 1)/(H + epoch*H)
  # CB1 <- sqrt(H * log(length(action_space) * H * n)) / sqrt(epoch)
  CB1 <- sqrt(H^5 * log(length(action_space) * H * n)) / sqrt(epoch)
  D_list <- c()
  current_stock <- 0
  for(h in 1:H){
    base_stock <- max(feasible_actions[[h]])
    current_stock <- max(base_stock, current_stock)
    D_list <- c(D_list, min(current_stock, train_D_mat[epoch, h]))
    current_stock <- max(0, current_stock - train_D_mat[epoch, h])
  }
  
  for(h in H:1){
    temp_feasible <- feasible_actions[[h]]
    for(a_idx in 1:length(temp_feasible)){
      # determine next stopping time
      temp_leftover <- action_space[temp_feasible[a_idx]]
      temp_leftover <- temp_leftover - D_list[h]
      moving_idx <- h + 1
      next_feasible <- feasible_actions[[moving_idx]]
      while(temp_leftover > max(action_space[next_feasible])){
        temp_leftover <- temp_leftover - D_list[[moving_idx]]
        moving_idx <- moving_idx + 1
        next_feasible <- feasible_actions[[moving_idx]]
      }
      r_term <- 0
      starting_level <- action_space[temp_feasible[a_idx]]
      for(temp_h in h:(moving_idx - 1)){
        r_term <- r_term + b_cost[temp_h] * min(starting_level, D_list[temp_h]) - h_cost[temp_h] * max(0, starting_level - D_list[temp_h]) 
        starting_level <- starting_level - D_list[temp_h]
      }
      
      aaaa <- which.min(abs(starting_level - action_space)) 
      
      if(moving_idx > H){
        Q_table[h, temp_feasible[a_idx]] <- (1 - alpha) * Q_table[h, temp_feasible[a_idx]] + alpha * r_term
      }else{
        Q_table[h, temp_feasible[a_idx]] <- (1 - alpha) * Q_table[h, temp_feasible[a_idx]] + alpha * (r_term + V_table[moving_idx, aaaa])
      }
    }
    
    posit <- which.max(Q_table[h, feasible_actions[[h]]])
    
    best_base_level[h] <- action_space[feasible_actions[[h]][posit]]  #action_space[which.max(Q_table[h, ])] 
    
    for(temp_ii in 1:length(action_space)){
      V_table[h, temp_ii] <- max(Q_table[h, temp_ii:length(action_space)])
    }
  }
  
  for(h in 1:H){
    old_feasible_action <- feasible_actions[[h]]
    new_feasible_action <- c()
    for(temp_jj in 1:length(old_feasible_action)){
      if(abs(Q_table[h, old_feasible_action[temp_jj]] - max(Q_table[h, feasible_actions[[h]]])) < CB1){
        new_feasible_action <- c(new_feasible_action, old_feasible_action[temp_jj])
      }
    }
    # print(c(h, new_feasible_action))
    feasible_actions[[h]] <- new_feasible_action
  }
  test_cost_HQL <- 0
  temp_level_HQL <- best_base_level[1]
  for(qqqq in 1:T){
    test_cost_HQL <- test_cost_HQL + (h_cost[qqqq]*pmax(0, temp_level_HQL - train_D_mat[epoch+1, qqqq]) + b_cost[qqqq]*pmax(0, train_D_mat[epoch+1, qqqq] - temp_level_HQL))
    if(qqqq<T){
      temp_level_HQL <- pmax(temp_level_HQL - train_D_mat[epoch+1, qqqq], best_base_level[qqqq + 1])
    }
  }
  
  cost_HQL_list[epoch+1] <- test_cost_HQL
  time_HQL_list[epoch+1] <- (proc.time() - start_time)["elapsed"]
  
  print(epoch)
}

cost_HQL <- cost_HQL_list

cumsum(cost_HQL[21:84])

time_HQL <- time_HQL_list

save(cost_HQL, time_HQL, file = "HQL_p1.rdata")
