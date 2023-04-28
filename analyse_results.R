library(fst)
library(data.table)
library(ggplot2)

# make the smoke prevalence thingy - 1hr
dt = read.fst("./results.fst")
dt = setDT(dt)
smoke_prev = fread("./smk_output.csv", data.table = TRUE)

# analyse data
colnames = c('mc', 'friendly_name', 'year', 'sex', 'agegrp', 'qimd', 'ethnicity', 
             'pops',
             'all_cause_mrtl')

names = grep(c("^cpp|^cypp|^dpp"), names(dt), value = TRUE)
colnames = c(colnames, names)

tt = setDT(dt[, ..colnames])
setkey(tt, year)
# sum all cpp, cypp
tt[, `:=`(
  cpp_all_cause = 
    cpp_chd + 
    cpp_stroke + 
    cpp_poststroke_dementia +
    cpp_copd +
    cpp_t2dm +
    cpp_lung_ca +
    cpp_colon_ca +
    cpp_breast_ca, 
  cypp_all_cause = 
    cypp_chd + 
    cypp_stroke + 
    cypp_poststroke_dementia +
    cypp_copd +
    cypp_t2dm +
    cypp_lung_ca +
    cypp_colon_ca +
    cypp_breast_ca)]

# cumsum 
cumsum_dt = tt[, `:=`(pops_cml = round(cumsum(pops)), 
                      all_cause_mrtl_cml = round(cumsum(all_cause_mrtl)), 
                      cpp_cvd_cml = round(cumsum(cpp_cvd)),
                      cpp_chd_cml = round(cumsum(cpp_chd)),
                      cpp_stroke_cml = round(cumsum(cpp_stroke)),
                      cpp_poststroke_dementia_cml = round(cumsum(cpp_poststroke_dementia)),
                      cpp_copd_cml = round(cumsum(cpp_copd)),
                      cpp_t2dm_cml = round(cumsum(cpp_t2dm)),
                      cpp_lung_ca_cml = round(cumsum(cpp_lung_ca)),
                      cpp_colon_ca_cml = round(cumsum(cpp_colon_ca)),
                      cpp_breast_ca_cml = round(cumsum(cpp_breast_ca)),
                      cypp_cvd_cml = round(cumsum(cypp_cvd)),
                      cypp_chd_cml = round(cumsum(cypp_chd)),
                      cypp_stroke_cml = round(cumsum(cypp_stroke)),
                      cypp_poststroke_dementia_cml = round(cumsum(cypp_poststroke_dementia)),
                      # cypp_af_cml = round(cumsum(cypp_af)), 
                      cypp_copd_cml = round(cumsum(cypp_copd)),
                      cypp_t2dm_cml = round(cumsum(cypp_t2dm)),
                      cypp_lung_ca_cml = round(cumsum(cypp_lung_ca)),
                      cypp_colon_ca_cml = round(cumsum(cypp_colon_ca)),
                      cypp_breast_ca_cml = round(cumsum(cypp_breast_ca)),
                      dpp_nonmodelled_cml = round(cumsum(dpp_nonmodelled)),
                      dpp_chd_cml = round(cumsum(dpp_chd)),
                      dpp_stroke_cml = round(cumsum(dpp_stroke)),
                      dpp_copd_cml = round(cumsum(dpp_copd)),
                      dpp_lung_ca_cml = round(cumsum(dpp_lung_ca)),
                      dpp_colon_ca_cml = round(cumsum(dpp_colon_ca)),
                      dpp_breast_ca_cml = round(cumsum(dpp_breast_ca)),
                      dpp_all_cause_cml = round(cumsum(dpp_all_cause)),
                      cpp_all_cause_cml = round(cumsum(cpp_all_cause)),
                      cypp_all_cause_cml = round(cumsum(cypp_all_cause))
),
by = c("friendly_name", "mc", "sex", "agegrp", "qimd", "ethnicity")]

# filter by timepoint
cumsum_dt = cumsum_dt[year %in% c(2030, 2047, 2072)] #starting 2023

# wide to long format
cumsum_dt = melt(cumsum_dt, id.vars = c("friendly_name", "mc","year", "sex", "agegrp", "qimd", "ethnicity"),
                 value.name = "value",
                 variable.name = "variable") # long to wide format

cumsum_dt[, auxv := paste0(variable,'~', friendly_name,'~', year,'~', sex,'~', agegrp,'~', qimd,'~', ethnicity)]
write_fst(cumsum_dt, "cumsum_dt.fst", 100L)

cumsum_dt = read.fst("./cumsum_dt.fst", as.data.table = TRUE)
setkey(cumsum_dt, auxv)
ttt = cumsum_dt[, CKutils::fquantile_byid(value, q = c(0.5, 0.025, 0.975), id = auxv, rounding = TRUE)]
setnames(ttt, new = c("auxv", "median", "lui", 'hui'))
ttt[, c("variable", "friendly_name", "year", "sex", "agegrp", "qimd", "ethnicity") := 
      tstrsplit(auxv, '~', fixed = TRUE)]
setcolorder(ttt, )

write.csv(cumsum_dt, "analysed_result.csv" )
