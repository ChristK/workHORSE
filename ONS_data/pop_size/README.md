# Population size inputs

workHORSE scales the synthetic population to ONS population counts by year,
single year of age and sex, using two files built here:

| File | Contents |
|---|---|
| `pop_proj.fst` | `year, age, sex, LAD23CD, LAD23NM, pops`; 296 local authorities (April 2023 boundaries), ages 0–90 (90 = 90+), 2003–2047. Used for regions and local authorities (summed over the selected LAs). |
| `pop_proj_england.fst` | `year, age, sex, pops`; England, ages 0–90, 2003–2047. Used whenever "England" is selected. |

## Sources

| Years | Local authorities (`pop_proj.fst`) | England (`pop_proj_england.fst`) |
|---|---|---|
| 2003–2010 | LSOA 2011 mid-year estimates (`synthpop/lsoa_mid_year_population_estimates.fst`) summed to LAD23 | Sum of LAs |
| 2011–2025 | ONS mid-year estimates, *mid-2011 to mid-2025 detailed time series* (MYEB1) | Sum of LAs (equals the ONS England total) |
| 2026–2047 | ONS **2022-based** subnational population projections (SNPP), **migration category variant**, April 2023 LA boundaries | ONS **2024-based** national population projections (NPP), principal projection |

- **SNPP 2022-based**, released 24 June 2025, mid-2022 to mid-2047. This
  release has no principal projection; ONS uses the migration category variant
  in its place as the headline projection.
  [Dataset (Z1)](https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationprojections/datasets/localauthoritiesinenglandz1),
  file `2022snpppopulationsyoamigcat23.zip`.
- **Mid-year estimates**, released 29 July 2026; mid-2022 to mid-2024 were
  revised for updated international migration.
  [Dataset](https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationestimates/datasets/estimatesofthepopulationforenglandandwales),
  file `myebtablesenglandwales20112025.xlsx`, sheet MYEB1.
- **NPP 2024-based**, released 28 April 2026, England, principal projection
  (`en_ppp_machine_readable.xlsx`, sheet Population).
  [Dataset (Z3)](https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationprojections/datasets/z3zippedpopulationprojectionsdatafilesengland).

Because England comes from the 2024-based NPP and local authorities from the
2022-based SNPP, from 2026 onwards the population of England is **not** the
sum of its regions or local authorities.

All ONS data are © Crown copyright, reused under the
[Open Government Licence v3.0](http://www.nationalarchives.gov.uk/doc/open-government-licence).

## Regenerating

Run from the project root, in this order:

1. `ONS_data/pop_size/download_ons_pop.R`: downloads the ONS files into
   `ONS_data/pop_size/raw/` (gitignored) and writes the compact extracts
   `ons_snpp2022_migcat_lad23.fst`, `ons_mye_2011_2025_lad23.fst` and
   `ons_npp2024_ppp_england.fst`.
2. `ONS_data/pop_size/update_locality_indx_lad23.R`: adds `LAD23CD`/`LAD23NM`
   to `synthpop/lsoa_to_locality_indx.fst` and writes
   `lad17_to_lad23_lookup.csv`.
3. `ONS_data/pop_size/transform_pops.R`: builds `pop_proj.fst` and
   `pop_proj_england.fst`, checks them, and writes junction diagnostics to
   `ONS_data/pop_size/diagnostics/`.

Cached synthetic populations store their weights, so they are regenerated
automatically when either population file changes (the files' md5 is part of
the synthpop checksum).
