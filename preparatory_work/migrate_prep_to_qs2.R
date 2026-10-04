# Migrate preparatory_work/*.R from qs to qs2 (run from the repo root of a LOCAL clone).
# Usage: Rscript preparatory_work/migrate_prep_to_qs2.R [--dry-run]
#
# The qs package is archived and cannot be installed on R >= 4.6. The model
# fitting scripts only need the full R serialisation of the fitted objects, so
# they move to qs2::qs_save()/qs_read() (never qd_save(): it cannot store model
# objects). Each file keeps its own line endings (fit_ckd_model.R is CRLF).
# Commit the migration together with "*.qs2" in .gitignore (the fitted models
# contain HSE microdata).
dry <- "--dry-run" %in% commandArgs(TRUE)
files <- list.files("preparatory_work", "\\.R$", full.names = TRUE)
files <- setdiff(files, file.path("preparatory_work", "migrate_prep_to_qs2.R")) # not itself
rules <- list(                                   # applied in this order
  c('qs::qread\\(',                     'qs2::qs_read('),
  c('\\bqsave\\(',                      'qs_save('),
  c('\\bqread\\(',                      'qs_read('),
  c(',\\s*preset\\s*=\\s*"high"',       ', compress_level = 4L'),   # qs "high" = zstd level 4
  c('(\\.qs)"',                         '\\12"'),                   # "....qs" -> "....qs2"
  c('"qs"',                             '"qs2"')                    # dependencies(c(..., "qs", ...))
)
# Line ending of a file: "\r\n" if every line ends in CRLF, "\n" if none does
line_ending <- function(f) {
  b <- readBin(f, "raw", file.size(f))
  n_lf <- sum(b == as.raw(10L))
  n_crlf <- sum(b[-1L] == as.raw(10L) & b[-length(b)] == as.raw(13L))
  if (n_crlf > 0L && n_crlf < n_lf) stop(f, ": mixed LF/CRLF line endings, fix by hand")
  if (n_crlf > 0L) "\r\n" else "\n"
}
tot <- 0L
for (f in files) {
  x <- readLines(f, warn = FALSE); y <- x       # readLines drops the CR of CRLF
  for (r in rules) y <- gsub(r[1], r[2], y, perl = TRUE)
  n <- sum(x != y)
  if (n == 0L) next
  tot <- tot + n
  eol <- line_ending(f)                         # before the file is overwritten
  invisible(parse(text = y, keep.source = FALSE))           # fails loudly on a syntax error
  left <- grep('\\bq(save|read)\\(|"qs"|\\.qs"|preset\\s*=', y, value = TRUE, perl = TRUE)
  if (length(left)) stop(f, ": unconverted lines remain:\n", paste(left, collapse = "\n"))
  cat(sprintf("%-55s %3d lines changed\n", f, n))
  if (!dry) {
    con <- file(f, "wb")                        # binary: no \n -> \r\n translation on Windows
    writeLines(y, con, sep = eol)
    close(con)
  }
}
cat("Total lines changed:", tot, if (dry) "(dry run, nothing written)" else "", "\n")
