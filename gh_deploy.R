library(piggyback)
tag <- "v0.0.2"
repo <- "ChristK/workHORSE"

if (interactive()) {
  print("interactive")
  pb_download(repo = repo, tag = tag, show_progress = TRUE)
} else { # used with Rscript
  # i.e. Rscript /root/workHORSE/gh_deploy.R "/root/workHORSE/"
  args <- commandArgs(TRUE)
  print(args)
  pb_download(dest = args, repo = repo, tag = tag, show_progress = FALSE)
}


#
