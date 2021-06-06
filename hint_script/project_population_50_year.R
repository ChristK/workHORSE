## project population to 2071 ( 2021 + 50)
## adapting file ./ONS_data/pop_size/transform_pops.R
## Vincy Huang
## Jun2021

## set-up ####
if (!require(CKutils)) {
  if (!require(remotes)) install.packages("remotes")
  remotes::install_github("ChristK/CKutils")
}
dependencies(c("Rcpp", "dqrng", "fst", "qs", "gamlss", "reldist", "future", "future.apply", "data.table"))
options(future.fork.enable = TRUE) # enable fork in Rstudio
plan(multicore)

# load data 
population_before_2041 <- read.fst("./ONS_data/pop_size/pop_proj.fst", as.data.table = TRUE)

# train model ####
population_train <- copy(population_before_2041)
population_train[, sex := factor(sex, levels = c("men", "women"), labels = c("men", "women"))]

model_population<- glm(pops ~ age + sex + year + LAD17CD , family = gaussian, data = population_train)
qsave(model_population, "./hint_script/population_prediction.qs")

# pedict till 2071 ####
model <- qread("./hint_script/population_prediction.qs")
population_predict <- CJ(age = unique(population_train$age), sex = unique(population_train$sex), 
              year = 2042: 2071)
test <- copy(population_predict)
test <- test[, LAD17CD := "E06000001"]
test <- test[, pops := predict.glm(model, newdata = test)]
test <- test[, .(LAD17CD, sex, age, year, pops)]
test_all <- rbind(population_train[, c(1, 3:6)], test)
data_plot <- test_all[pop_sum := summary()]
data_plot <- test_all[LAD17CD == "E06000001"][, .(pop_sum = sum(pops)), by = .(year, sex)]
ggplot() + 
  geom_point(data = data_plot, aes(x = year, y = pop_sum, group = sex))
