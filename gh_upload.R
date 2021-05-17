library(piggyback)
tag <- "v0.0.2"
repo <- "ChristK/workHORSE"
# pb_new_release(repo, tag)

pb_track(c("disease_epidemiology/disease_epi_l.fst",
  "lifecourse_models/*.fst",
  "simulation/health_econ/informal_care_costs_l.fst",
  "simulation/health_econ/productivity_costs_l.fst"
))

pb_upload(pb_track(), repo = repo, tag = tag)

pb_list(repo = repo, tag = tag)
