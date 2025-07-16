rm(list = ls())
bakery_data <- read.csv(file = "Bakery sales.csv", header = TRUE)
dates_all <- unique(bakery_data$date) 
products_all <- unique(bakery_data$article)

# # Count occurrences of each level
# article_counts <- table(bakery_data$article)
# 
# # Filter levels that appear more than 1000 times
# frequent_articles <- names(article_counts[article_counts > 5000])
# 
# # View the result
# print(frequent_articles)

my_data <- as.data.frame(matrix(0, 600, length(products_all) + 1))

for(i in 1:length(dates_all)){
  temp_data <- bakery_data[which(bakery_data$date == dates_all[i]), ]
  for(j in 1:nrow(temp_data)){
    temp_a <- which(temp_data$article[j] == products_all) 
    if(length(temp_a) > 0){
      my_data[i, 1+temp_a] <- my_data[i, 1+temp_a] + temp_data$Quantity[j]
    }
  }
  my_data[i, 1] <- dates_all[i]
}

my_data$V1 <- as.Date(my_data$V1)
my_data$dow <- weekdays(my_data$V1)

table(my_data$dow)

my_data <- my_data[my_data$dow != "Wednesday", ]
my_data <- my_data[-c(1:4), ]

dow_mat <- matrix(my_data$dow, 89, 6, byrow = TRUE)
date_mat <- matrix(as.character(my_data$V1), 89, 6, byrow = TRUE)

my_data <- my_data[-which((my_data$V1 >= as.Date("2021-05-13"))&(my_data$V1 <= as.Date("2021-05-25"))), ]
my_data <- my_data[-c(524, 525), ]

dow_mat <- matrix(my_data$dow, 88, 6, byrow = TRUE)
date_mat <- matrix(as.character(my_data$V1), 88, 6, byrow = TRUE)

my_data <- my_data[-which((my_data$V1 >= as.Date("2021-12-02"))&(my_data$V1 <= as.Date("2021-12-05"))), ]

dow_mat <- matrix(my_data$dow, 87, 6, byrow = TRUE)
date_mat <- matrix(as.character(my_data$V1), 87, 6, byrow = TRUE)

my_data <- my_data[-which((my_data$V1 >= as.Date("2021-12-30"))&(my_data$V1 <= as.Date("2022-01-04"))), ]

dow_mat <- matrix(my_data$dow, 86, 6, byrow = TRUE)
date_mat <- matrix(as.character(my_data$V1), 86, 6, byrow = TRUE)


my_data <- my_data[-which((my_data$V1 >= as.Date("2022-09-02"))&(my_data$V1 <= as.Date("2022-09-06"))), ]

dow_mat <- matrix(my_data$dow, 85, 6, byrow = TRUE)
date_mat <- matrix(as.character(my_data$V1), 85, 6, byrow = TRUE)


my_data <- my_data[-which((my_data$V1 >= as.Date("2022-09-15"))&(my_data$V1 <= as.Date("2022-09-20"))), ]

my_data <- my_data[, c(1, 151, 2:150)]

names(my_data) <- c("Date", "DoW", products_all)

indexing_order <- order(apply(my_data[, -c(1:2)], MARGIN = 2, sum), decreasing = TRUE)

save(my_data, indexing_order, products_all, file = "bakery_data.rdata")