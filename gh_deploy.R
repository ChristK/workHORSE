library(piggyback)
tag <- "v0.0.2"
repo <- "ChristK/workHORSE"
# pb_new_release(repo, tag)

pb_download(repo = repo, tag = tag, show_progress = FALSE)
