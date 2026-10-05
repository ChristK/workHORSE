# Downloads the large data files that are not kept in git (see .gitignore) from
# the GitHub release below, and checks each one against the size and md5 listed
# in gh_deploy_files.csv. Files that are already present and intact are skipped,
# so it is safe to run again.
#
# It uses the public release download links, not the GitHub API, so the API
# rate limit (60 requests an hour without a token) does not apply. It stops with
# an error if any file cannot be downloaded intact, so the Docker build, which
# runs it, fails instead of producing an image with broken data.
#
# Usage, from the repo root in R: source("gh_deploy.R")
#   or from a shell:              Rscript gh_deploy.R /path/to/workHORSE
#
# After uploading new files with gh_upload.R, update release_url and regenerate
# gh_deploy_files.csv from the repo root:
#   f <- c("disease_epidemiology/disease_epi_l.fst",
#          "simulation/health_econ/informal_care_costs_l.fst",
#          "simulation/health_econ/productivity_costs_l.fst",
#          sort(list.files("lifecourse_models", "fst$", full.names = TRUE),
#               method = "radix"))
#   write.csv(data.frame(file = basename(f), dir = dirname(f),
#                        size = file.size(f), md5 = unname(tools::md5sum(f))),
#             "gh_deploy_files.csv", row.names = FALSE, quote = FALSE)

local({
  release_url <- "https://github.com/ChristK/workHORSE/releases/download/v0.0.2/"

  root <- if (interactive()) "." else commandArgs(trailingOnly = TRUE)[1L]
  if (is.na(root) || !nzchar(root)) root <- "."

  manifest <- file.path(root, "gh_deploy_files.csv")
  if (!file.exists(manifest)) {
    stop("Cannot find ", manifest, ". Run this from the workHORSE repo root, ",
         "or pass the repo path: Rscript gh_deploy.R /path/to/workHORSE",
         call. = FALSE)
  }
  files <- utils::read.csv(manifest, colClasses = "character")

  intact <- function(path, size, md5) {
    file.exists(path) &&
      file.size(path) == as.numeric(size) &&
      unname(tools::md5sum(path)) == md5
  }

  # download.file() gives up after getOption("timeout") seconds (default 60)
  op <- options(timeout = max(3600, getOption("timeout")))
  on.exit(options(op), add = TRUE)

  failed <- character(0)
  for (i in seq_len(nrow(files))) {
    f <- files[i, ]
    dest <- file.path(root, f$dir, f$file)
    if (intact(dest, f$size, f$md5)) next

    dir.create(dirname(dest), showWarnings = FALSE, recursive = TRUE)
    tmp <- paste0(dest, ".part")
    size <- as.numeric(f$size)
    ok <- FALSE
    for (attempt in 1:3) {
      if (attempt > 1L) Sys.sleep(10)
      message("Downloading ", f$file, " (",
              if (size < 1e6) sprintf("%.1f kB", size / 1e3)
              else sprintf("%.1f MB", size / 1e6), ")",
              if (attempt > 1L) paste0(", attempt ", attempt))
      tryCatch(
        utils::download.file(paste0(release_url, f$file), tmp,
                             mode = "wb", quiet = TRUE),
        error = function(e) message("  ", conditionMessage(e))
      )
      if (intact(tmp, f$size, f$md5)) {
        ok <- file.rename(tmp, dest)
        break
      }
      message("  ", f$file, " is missing or corrupt after the download")
    }
    if (!ok) {
      unlink(tmp)
      failed <- c(failed, f$file)
    }
  }

  if (length(failed) > 0L) {
    stop(length(failed), " of ", nrow(files),
         " data files could not be downloaded intact: ",
         paste(failed, collapse = ", "), call. = FALSE)
  }
  message("All ", nrow(files), " data files are present and intact.")
})
