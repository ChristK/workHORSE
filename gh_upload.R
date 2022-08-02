library(piggyback)
tag <- "v0.0.2"
repo <- "ChristK/workHORSE"
# pb_release_create(repo, tag)

pb_upload("disease_epi_l.fst", repo = repo, tag = tag, dir = "disease_epidemiology/")
pb_upload("informal_care_costs_l.fst", repo = repo, tag = tag, dir = "simulation/health_econ/")
pb_upload("productivity_costs_l.fst", repo = repo, tag = tag, dir = "simulation/health_econ/")
pb_upload(list.files("lifecourse_models", ".fst$"), repo = repo, tag = tag, dir = "lifecourse_models/")

pb_list(repo = repo, tag = tag)
