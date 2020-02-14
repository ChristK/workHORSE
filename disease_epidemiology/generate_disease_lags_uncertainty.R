cat("Initialising workHORSE model...\n\n")
if (!require(CKutils)) {
  if (!require(remotes)) install.packages("remotes")
  remotes::install_github("ChristK/CKutils")
  library(CKutils)
}
dependencies(c("data.table", "fst"))

set.seed(2506120L)

mc_max  <- 2e3
mc_max2 <- 10L # To allow 2 or more diseases have the same lag in years but be independent

dt <- CJ(mc = 1:mc_max, disease_enum = 1:mc_max2)
for (i in 3:9) {
  set(dt, NULL, paste0("lag_", i), 2L + rbinom(nrow(dt), 8L, (i - 2) / 8L))
}
write_fst(dt, "./disease_epidemiology/disease_lags_l.fst", 100L)

# create a table with row numbers for each mc
dt[, rn := .I]
tt <- dt[, .(from = min(rn), to = max(rn)), keyby = mc]
write_fst(tt, "./disease_epidemiology/disease_lags_indx.fst", 100L)

# tt <- 2L + rbinom(1e5, 8L, (9L - 2L) / 8L)
# mean(tt)
# hist(tt)
