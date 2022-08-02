library(piggyback)
tag <- "v0.0.2"
repo <- "ChristK/workHORSE"



if (interactive()) {
  print("interactive")
  dir.create("./disease_epidemiology", showWarnings = FALSE, recursive = TRUE)
  dir.create("./simulation/health_econ", showWarnings = FALSE, recursive = TRUE)
  dir.create("./lifecourse_models", showWarnings = FALSE, recursive = TRUE)


  pb_download("disease_epi_l.fst",
              "disease_epidemiology/",
              repo = repo,
              tag = tag)
  pb_download(
    "informal_care_costs_l.fst",
    "simulation/health_econ/",
    repo = repo,
    tag = tag
  )
  pb_download(
    "productivity_costs_l.fst",
    "simulation/health_econ/",
    repo = repo,
    tag = tag
  )
  pb_download(
    NULL,
    "lifecourse_models/",
    ignore = c(
      "disease_epi_l.fst",
      "informal_care_costs_l.fst",
      "productivity_costs_l.fst"
    ),
    repo = repo,
    tag = tag
  )
} else { # used with Rscript
  # i.e. Rscript /root/workHORSE/gh_deploy.R "/root/workHORSE/"
  args <- commandArgs(TRUE)

  dir.create(file.path(args, "disease_epidemiology"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(args, "simulation", "health_econ"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(args, "lifecourse_models"), showWarnings = FALSE, recursive = TRUE)

  pb_download("disease_epi_l.fst",
              file.path(args, "disease_epidemiology"),
              repo = repo,
              tag = tag)
  pb_download(
    "informal_care_costs_l.fst",
    file.path(args, "simulation", "health_econ"),
    repo = repo,
    tag = tag
  )
  pb_download(
    "productivity_costs_l.fst",
    file.path(args, "simulation", "health_econ"),
    repo = repo,
    tag = tag
  )
  pb_download(
    NULL,
    file.path(args, "lifecourse_models"),
    ignore = c(
      "disease_epi_l.fst",
      "informal_care_costs_l.fst",
      "productivity_costs_l.fst"
    ),
    repo = repo,
    tag = tag
  )
}


#
