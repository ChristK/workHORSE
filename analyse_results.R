library(fst)
library(data.table)
library(ggplot2)

# make the smoke prevalence thingy - 1hr
dt = read.fst("/mnt/storage_fast/output/hint/results.fst", as.data.table = TRUE)
smoke_prev = fread("/mnt/storage_fast/output/hint/smk_output.csv", data.table = TRUE)

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
#cumsum_dt = cumsum_dt[year %in% c(2030, 2047, 2072)] #starting 2023
cumsum_dt = cumsum_dt[year %in% c(2030, 2047, 2072, 2045, 2046)] #starting 2023

# wide to long format
cumsum_dt = melt(cumsum_dt, id.vars = c("friendly_name", "mc","year", "sex", "agegrp", "qimd", "ethnicity"),
                 value.name = "value",
                 variable.name = "variable") # long to wide format

# format to different table

ttt = copy(cumsum_dt)
ttt[, auxv := paste0(variable,'~', friendly_name,'~', year,'~', sex,'~', agegrp,'~', qimd,'~', ethnicity)]
setkey(ttt, auxv)
quantile_dt = ttt[, CKutils::fquantile_byid(value, q = c(0.5, 0.025, 0.975), id = auxv, rounding = TRUE)]
setnames(quantile_dt, new = c("auxv", "median", "lui", 'hui'))
quantile_dt[, c("variable", "friendly_name", "year", "sex", "agegrp", "qimd", "ethnicity") := 
              tstrsplit(auxv, '~', fixed = TRUE)]
quantile_dt = quantile_dt[, c("year", "friendly_name", "sex" ,"agegrp" ,"qimd" ,"ethnicity", "variable" ,
                              "median", "lui", "hui" )]
write.csv(quantile_dt, "/mnt/storage_fast/output/hint/output/quantile_by_all_variable.csv" )

remove(ttt)
remove(quantile_dt)

ttt = copy(cumsum_dt)
ttt[, auxv := paste0(variable,'~', friendly_name,'~', year,'~', qimd)]
setkey(ttt, auxv)
ttt = ttt[, value := sum(value), by = .(auxv, variable, mc, friendly_name, year, sex, agegrp, qimd)]
quantile_dt = ttt[, CKutils::fquantile_byid(value, q = c(0.5, 0.025, 0.975), id = auxv, rounding = TRUE)]
setnames(quantile_dt, new = c("auxv", "median", "lui", 'hui'))
quantile_dt[, c("variable", "friendly_name", "year",  "qimd") := 
              tstrsplit(auxv, '~', fixed = TRUE)]
quantile_dt = quantile_dt[, c("year", "friendly_name", "qimd" , "variable" ,
                              "median", "lui", "hui" )]
quantile_dt[year == 2072 & variable == "cpp_all_cause_cml" & friendly_name == "SSS and taxation"]
quantile_dt[year == 2072 & variable == "cpp_all_cause_cml" & friendly_name == "SSS and taxation", signif(median,2)]
quantile_dt[year == 2072 & variable == "cpp_all_cause_cml" & friendly_name == "SSS and taxation", signif(lui,2)]
quantile_dt[year == 2072 & variable == "cpp_all_cause_cml" & friendly_name == "SSS and taxation", signif(hui,2)]
quantile_dt[year == 2072 & variable == "cypp_all_cause_cml" & friendly_name == "taxation"]
write.csv(quantile_dt, "/mnt/storage_fast/output/hint/output/quantile_by_qimd_year.csv" )

remove(ttt)
ttt = copy(cumsum_dt)
ttt[, auxv := paste0(variable,'~', friendly_name,'~', year)]
setkey(ttt, auxv)
ttt = ttt[, value := sum(value), by = .(auxv, mc, variable)]
quantile_dt = ttt[, CKutils::fquantile_byid(value, q = c(0.5, 0.025, 0.975), id = auxv, rounding = TRUE)]
setnames(quantile_dt, new = c("auxv", "median", "lui", 'hui'))
quantile_dt[, c("variable", "friendly_name", "year") := 
              tstrsplit(auxv, '~', fixed = TRUE)]
quantile_dt = quantile_dt[, c("year", "friendly_name", "variable" ,
                              "median", "lui", "hui" )]
quantile_dt[year %in% c(2046, 2045) & friendly_name %in% c("MALA 21", "SSS") & variable == "cpp_all_cause_cml"]
quantile_dt[year == 2045 & friendly_name %in% c("MALA 21", "SSS") & variable == "cpp_all_cause_cml", signif(median,2)]
quantile_dt[year == 2045 & friendly_name %in% c("MALA 21", "SSS") & variable == "cpp_all_cause_cml", signif(lui,2)]
quantile_dt[year == 2045 & friendly_name %in% c("MALA 21", "SSS") & variable == "cpp_all_cause_cml", signif(hui,2)]
write.csv(quantile_dt, "/mnt/storage_fast/output/hint/output/quantile_by_year.csv" )


