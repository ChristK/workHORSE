
# Generate uncertainty (2e3 mc iterations)
# asumming 97.5 precentile is up to 20% higher than the median for incid, prev, fatal
# No uncertainty is assumed for disease duration
# NOTE Currently we assume that uncertainty covaries for all disease, age, sex, qimd,
# so when i.e. chd incid is low then lung ca inc is also low
mc_max <- 2e3
uncertainty_space <- 1.2


library(data.table)
library(fst)
library(future.apply)
library(SimJoint)
library(dqrng)
library(rriskDistributions)
if (!require(workHORSEmisc)) {
  if (!require(remotes))
    install.packages("remotes")
  remotes::install_local("./Rpackage/workHORSE_model_pkg/")
  library(workHORSEmisc)
}

options(
  future.fork.enable = TRUE, # enable fork in Rstudio
  future.globals.maxSize = +Inf,
  future.globals.onReference = "ignore"
)
plan(multiprocess, workers = 20L)


disease_epi <- fread("./disease_epidemiology/disease_epi.csv", stringsAsFactors = TRUE,
                     select = c("age", "sex", "qimd", "disease", "out_incidence_rates",
                                "out_prevalence_rates", "out_case_fatality_rates",
                                "out_duration_years"))
disease_epi <- melt(disease_epi, c("age", "sex", "qimd", "disease"), variable.name = "type")

disease_epi[, disease := gsub(" ", "_", disease)]
disease_epi[, disease := gsub("cancer", "ca", disease)]
disease_epi[, type := gsub("out_", "", type)]

get_beta_par <- function(mu, scale = 1.2, tol = 0.001) { # scale of 1.5 is 50% higher
  stopifnot(scale > 1)
  suppressMessages(rriskDistributions::get.beta.par(c(0.5, 0.975), c(mu, mu * scale), FALSE, FALSE, tol))
}

# too small values do not converge I need to multiply scaling factor here and
# divide just before the end
disease_epi[, mltp_factor := 1]
disease_epi[value * uncertainty_space < 1e-2, mltp_factor := 1e2]
disease_epi[value * uncertainty_space < 1e-3, mltp_factor := 1e3]
disease_epi[value * uncertainty_space < 1e-4, mltp_factor := 1e4]
disease_epi[value * uncertainty_space < 1e-5, mltp_factor := 1e5]

disease_epi[type == "incidence_rates", c("a", "b") := {
  l <- future_lapply(value * mltp_factor, get_beta_par, uncertainty_space, 0.001)
  list(sapply(l, `[[`, 1), sapply(l, `[[`, 2))
}]
disease_epi[type == "prevalence_rates", c("a", "b") := {
  l <- future_lapply(value * mltp_factor, get_beta_par, uncertainty_space, 0.001)
  list(sapply(l, `[[`, 1), sapply(l, `[[`, 2))
}]
disease_epi[type == "case_fatality_rates" & value == 0, `:=`(value = 1e-8, mltp_factor = 1e6)] # replaces 0 fatality
disease_epi[type == "case_fatality_rates", c("a", "b") := {
  l <- future_lapply(value * mltp_factor, get_beta_par, uncertainty_space, 0.001)
  list(sapply(l, `[[`, 1), sapply(l, `[[`, 2))
}]


disease_epi_l <- CJ(mc = 1:mc_max, age = disease_epi$age, sex = disease_epi$sex,
                    qimd = disease_epi$qimd, type = disease_epi$type,
                    disease = disease_epi$disease, unique = TRUE)
disease_epi_l <- disease_epi_l[!(sex == "men" & disease == "breast_ca")]

absorb_dt(disease_epi_l, disease_epi)

# Generate swarm of correlated random numbers
dqRNGkind("pcg64") # dqRNGkind("Xoroshiro128+") ~10% faster
SEED <- 8003283L # sample(1e7, 1)
set.seed(SEED) # Current is to ensure reproducibility.
dqset.seed(SEED)# Ensure that seed hasn't been used elsewher in the model


# assume 0.6 correlation between prev and mort, 0.5 between incid and prev. 0.4 cor between incid and mort
# 1, 2, 3 is incid, prev, mort respectively
cm <- matrix(c(1, 0.5, 0.4, 0.5, 1, 0.6, 0.4, 0.6, 1), 3)
# cm <- matrix(c(1, 0.1, 0.1, 0.1, 1, 0.1, 0.1, 0.1, 1), 3)

colnames(cm) <- c("incd", "prvl", "mrtl")
rownames(cm) <- c("incd", "prvl", "mrtl")
tt <- generate_corr_unifs(mc_max, cm)
# cor(tt)
rank_mtx <- data.table(mc = 1:mc_max, tt)

absorb_dt(disease_epi_l, rank_mtx)
disease_epi_l[type == "incidence_rates", p := incd]
disease_epi_l[type == "prevalence_rates", p := prvl]
disease_epi_l[type == "case_fatality_rates", p := mrtl] # despite I store it with fatality this is really used in mortality
disease_epi_l[, c("incd", "prvl", "mrtl") := NULL]
setnames(disease_epi_l, "value", "median_value")
disease_epi_l[, value := qbeta(p, a, b)/mltp_factor, by = .(age, qimd, sex)]
# disease_epi_l[type != "duration_years", anyNA(value)]

disease_epi_l[, c("a", "b", "mltp_factor") := NULL]
disease_epi_l <- dcast(disease_epi_l, mc+age+sex+qimd~disease+type, drop = c(TRUE, TRUE),
                       value.var = c("value", "median_value", "p"), fill = 0)
tt <- disease_epi_l[, transpose(lapply(.SD, function(x) all(is.na(x))), keep.names = "rn")][(V2)]
disease_epi_l[, (tt$V1) := NULL]



write_fst(disease_epi_l, "./disease_epidemiology/disease_epi_l.fst", 100L)
# create a table with row numbers for each mc
disease_epi_l[, rn := .I]
tt <- disease_epi_l[, .(from = min(rn), to = max(rn)), keyby = mc]
write_fst(tt, "./disease_epidemiology/disease_epi_indx.fst", 100L)


