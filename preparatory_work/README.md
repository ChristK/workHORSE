# preparatory_work

Offline scripts that fit the lifecourse models (`fit_*_model.R`), extract the HSE
correlation structure and build other model inputs. They are run by hand, not by
the app.

## The qs package

`qs` is archived on CRAN and cannot be installed on R >= 4.6. The model and the
Shiny app no longer use it (they use base R `saveRDS()`/`readRDS()`), but these
scripts still call `qsave()`/`qread()` and list `"qs"` in `dependencies()`.

Before re-fitting any model, migrate the scripts to `qs2`, from the repo root:

```
Rscript preparatory_work/migrate_prep_to_qs2.R --dry-run   # preview, writes nothing
Rscript preparatory_work/migrate_prep_to_qs2.R
```

The migration script rewrites every script that uses `qs`, keeping each file's
line endings (`fit_ckd_model.R` is CRLF):

- `qsave()`/`qread()` become `qs_save()`/`qs_read()`. Never use `qd_save()` for the
  fitted model objects: it replaces formulas, functions and environments with
  `NULL`.
- `preset = "high"` becomes `compress_level = 4L`.
- `.qs` file names become `.qs2`, and `"qs"` becomes `"qs2"` in `dependencies()`.

Add `*.qs2` to `.gitignore` in the same commit as the migration. The fitted model
objects contain HSE microdata (`model$data`), which is why `*.qs` is ignored today.

## Existing name mismatches in fit_af_diag_model.R

The migration does not touch these. Fix them when re-fitting.

- Line 57 checks for `marginal_distr_af.qs`, but lines 58 and 61 read and write
  `marginal_distr_af_diag.qs`.
- Line 165 saves `af_dgn_model.qs`, but line 190 reads `af_diag_model.qs`.
