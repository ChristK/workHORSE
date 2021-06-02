parameters <- qread("./output/parameters.qs")
input <- qread("./output/input.qs")

out <- function() read_fst("./output/results.fst", as.data.table = TRUE)[]

out_proc <- function() {
    agegroup_filter <-
      c("30-49", "50-69", "70-89")[c("30-49", "50-69", "70-89") %in% input$out_characteristics_select]
    sex_filter <-
      c("men", "women")[c("men", "women") %in% input$out_characteristics_select]
    qimd_filter <-
      c("1 most deprived", "2", "3", "4", "5 least deprived")[c("1 most deprived", "2", "3", "4", "5 least deprived") %in% input$out_characteristics_select]
    ethn_filter <- c(
      "white",
      "indian",
      "pakistani",
      "bangladeshi",
      "other asian",
      "black caribbean",
      "black african",
      "chinese",
      "other"
    )[c(
      "white",
      "indian",
      "pakistani",
      "bangladeshi",
      "other asian",
      "black caribbean",
      "black african",
      "chinese",
      "other"
    ) %in% input$out_characteristics_select]


    dt <- out()[year <= input$inout_year_slider &
        # friendly_name %in% input$inout_scenario_select &
        agegrp %in% agegroup_filter &
        sex %in% sex_filter &
        qimd %in% qimd_filter &
        ethnicity %in% ethn_filter,
    ]

    dt <- sum_dt(dt)

    # Discounting costs and utility
    nam <- grep("_utility$", names(dt), value = TRUE)
    for (j in nam) {
      set(dt, NULL, j, deflate(dt[[j]], input$out_discount_qalys_slider, dt$year, 2019L))
    }
    nam <- grep("_cost$", names(dt), value = TRUE)
    for (j in nam) {
      set(dt, NULL, j, deflate(dt[[j]], input$out_discount_costs_slider, dt$year, 2019L))
    }

    setkey(dt, year, friendly_name, mc)
    dt[,
      `:=` (
        invitation_cost_cml = round(cumsum(invitation_cost)),
        attendees_cost_cml = round(cumsum(attendees_cost)),
        active_days_cost_cml = round(cumsum(active_days_cost)),
        bmi_cost_cml = round(cumsum(bmi_cost)),
        alcohol_cost_cml = round(cumsum(alcohol_cost)),
        smoking_cost_cml = round(cumsum(smoking_cost)),
        nonmodelled_mrtl_cml = round(cumsum(nonmodelled_mrtl)),
        cvd_prvl_cml = round(cumsum(cvd_prvl)),
        stroke_prvl_cml = round(cumsum(stroke_prvl)),
        stroke_dgn_cml = round(cumsum(stroke_dgn)),
        stroke_mrtl_cml = round(cumsum(stroke_mrtl)),
        chd_prvl_cml = round(cumsum(chd_prvl)),
        chd_dgn_cml = round(cumsum(chd_dgn)),
        chd_mrtl_cml = round(cumsum(chd_mrtl)),
        t2dm_prvl_cml = round(cumsum(t2dm_prvl)),
        t2dm_dgn_cml = round(cumsum(t2dm_dgn)),
        af_prvl_cml = round(cumsum(af_prvl)),
        af_dgn_cml = round(cumsum(af_dgn)),
        # htn_prvl_cml = round(cumsum(htn_prvl)),
        # htn_dgn_cml = round(cumsum(htn_dgn)),
        all_cause_mrtl_cml = round(cumsum(all_cause_mrtl)),
        eq5d_cml = cumsum(eq5d),
        healthcare_cost_cml = round(cumsum(healthcare_cost)),
        socialcare_cost_cml = round(cumsum(socialcare_cost)),
        productivity_cost_cml = round(cumsum(productivity_cost)),
        informal_care_cost_cml = round(cumsum(informal_care_cost)),
        cvd_incd_cml = round(cumsum(cvd_incd)),
        stroke_incd_cml = round(cumsum(stroke_incd)),
        chd_incd_cml = round(cumsum(chd_incd)),
        t2dm_incd_cml = round(cumsum(t2dm_incd)),
        af_incd_cml = round(cumsum(af_incd)),
        # htn_incd_cml = round(cumsum(htn_incd)),
        pops_cml = round(cumsum(pops)),
        smkcess_ovrhd_cost_cml = round(cumsum(smkcess_ovrhd_cost)),
        wghtloss_ovrhd_cost_cml = round(cumsum(wghtloss_ovrhd_cost)),
        pa_ovrhd_cost_cml = round(cumsum(pa_ovrhd_cost)),
        alcoholreduc_ovrhd_cost_cml = round(cumsum(alcoholreduc_ovrhd_cost)),
        policy_cost_cml = round(cumsum(policy_cost)),
        net_utility_cml = signif(cumsum(net_utility), 2L),
        net_policy_cost_cml = round(cumsum(net_policy_cost)),
        net_healthcare_cost_cml = round(cumsum(net_healthcare_cost)),
        net_socialcare_cost_cml = round(cumsum(net_socialcare_cost)),
        net_informal_care_cost_cml = round(cumsum(net_informal_care_cost)),
        net_productivity_cost_cml = round(cumsum(net_productivity_cost)),
        cpp_cvd_cml = round(cumsum(cpp_cvd)),
        cpp_chd_cml = round(cumsum(cpp_chd)),
        cpp_stroke_cml = round(cumsum(cpp_stroke)),
        cpp_poststroke_dementia_cml = round(cumsum(cpp_poststroke_dementia)),
        # cpp_htn_cml = round(cumsum(cpp_htn)),
        cpp_copd_cml = round(cumsum(cpp_copd)),
        # cpp_af_cml = round(cumsum(cpp_af)),
        cpp_t2dm_cml = round(cumsum(cpp_t2dm)),
        cpp_lung_ca_cml = round(cumsum(cpp_lung_ca)),
        cpp_colon_ca_cml = round(cumsum(cpp_colon_ca)),
        cpp_breast_ca_cml = round(cumsum(cpp_breast_ca)),

        cypp_cvd_cml = round(cumsum(cypp_cvd)),
        cypp_chd_cml = round(cumsum(cypp_chd)),
        cypp_stroke_cml = round(cumsum(cypp_stroke)),
        cypp_poststroke_dementia_cml = round(cumsum(cypp_poststroke_dementia)),
        # cypp_htn_cml = round(cumsum(cypp_htn)),
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
        total_hcp_cost_cml = round(cumsum(total_hcp_cost)),
        total_hscp_cost_cml = round(cumsum(total_hscp_cost)),
        societal_cost_cml = round(cumsum(societal_cost))
      ),
      by = c("friendly_name", "mc")
    ][, nmb_cml := net_monetary_benefit_cml(.SD, input$health_econ_perspective_checkbox, input$out_wtp_box)]
    # CBR cannot be summed because it is a ratio so it is inappropriate to put here
  }
out_proc()[]

out_proc_qimd <- function () {
    agegroup_filter <-
      c("30-49", "50-69", "70-89")[c("30-49", "50-69", "70-89") %in% input$out_characteristics_select]
    sex_filter <-
      c("men", "women")[c("men", "women") %in% input$out_characteristics_select]
    qimd_filter <-
      c("1 most deprived", "2", "3", "4", "5 least deprived")[c("1 most deprived", "2", "3", "4", "5 least deprived") %in% input$out_characteristics_select]
    ethn_filter <- c(
      "white",
      "indian",
      "pakistani",
      "bangladeshi",
      "other asian",
      "black caribbean",
      "black african",
      "chinese",
      "other"
    )[c(
      "white",
      "indian",
      "pakistani",
      "bangladeshi",
      "other asian",
      "black caribbean",
      "black african",
      "chinese",
      "other"
    ) %in% input$out_characteristics_select]


    dt <- out()[year <= input$inout_year_slider &
        # friendly_name %in% input$inout_scenario_select &
        agegrp %in% agegroup_filter &
        sex %in% sex_filter &
        qimd %in% qimd_filter &
        ethnicity %in% ethn_filter,
    ]

    dt <- sum_dt(dt, c("year", "friendly_name", "mc", "qimd"))


    # Discounting costs and utility
    nam <- grep("_utility$", names(dt), value = TRUE)
    for (j in nam) {
      set(dt, NULL, j, deflate(dt[[j]], input$out_discount_qalys_slider, dt$year, 2019L))
    }
    nam <- grep("_cost$", names(dt), value = TRUE)
    for (j in nam) {
      set(dt, NULL, j, deflate(dt[[j]], input$out_discount_costs_slider, dt$year, 2019L))
    }

    setkey(dt, year, friendly_name, qimd, mc)
    dt[,
      `:=` (
        invitation_cost_cml = round(cumsum(invitation_cost)),
        attendees_cost_cml = round(cumsum(attendees_cost)),
        active_days_cost_cml = round(cumsum(active_days_cost)),
        bmi_cost_cml = round(cumsum(bmi_cost)),
        alcohol_cost_cml = round(cumsum(alcohol_cost)),
        smoking_cost_cml = round(cumsum(smoking_cost)),
        nonmodelled_mrtl_cml = round(cumsum(nonmodelled_mrtl)),
        cvd_prvl_cml = round(cumsum(cvd_prvl)),
        stroke_prvl_cml = round(cumsum(stroke_prvl)),
        stroke_dgn_cml = round(cumsum(stroke_dgn)),
        stroke_mrtl_cml = round(cumsum(stroke_mrtl)),
        chd_prvl_cml = round(cumsum(chd_prvl)),
        chd_dgn_cml = round(cumsum(chd_dgn)),
        chd_mrtl_cml = round(cumsum(chd_mrtl)),
        t2dm_prvl_cml = round(cumsum(t2dm_prvl)),
        t2dm_dgn_cml = round(cumsum(t2dm_dgn)),
        af_prvl_cml = round(cumsum(af_prvl)),
        af_dgn_cml = round(cumsum(af_dgn)),
        # htn_prvl_cml = round(cumsum(htn_prvl)),
        # htn_dgn_cml = round(cumsum(htn_dgn)),
        all_cause_mrtl_cml = round(cumsum(all_cause_mrtl)),
        eq5d_cml = cumsum(eq5d),
        healthcare_cost_cml = round(cumsum(healthcare_cost)),
        socialcare_cost_cml = round(cumsum(socialcare_cost)),
        productivity_cost_cml = round(cumsum(productivity_cost)),
        informal_care_cost_cml = round(cumsum(informal_care_cost)),
        cvd_incd_cml = round(cumsum(cvd_incd)),
        stroke_incd_cml = round(cumsum(stroke_incd)),
        chd_incd_cml = round(cumsum(chd_incd)),
        t2dm_incd_cml = round(cumsum(t2dm_incd)),
        af_incd_cml = round(cumsum(af_incd)),
        # htn_incd_cml = round(cumsum(htn_incd)),
        pops_cml = round(cumsum(pops)),
        smkcess_ovrhd_cost_cml = round(cumsum(smkcess_ovrhd_cost)),
        wghtloss_ovrhd_cost_cml = round(cumsum(wghtloss_ovrhd_cost)),
        pa_ovrhd_cost_cml = round(cumsum(pa_ovrhd_cost)),
        alcoholreduc_ovrhd_cost_cml = round(cumsum(alcoholreduc_ovrhd_cost)),
        policy_cost_cml = round(cumsum(policy_cost)),
        net_utility_cml = signif(cumsum(net_utility), 2L),
        net_policy_cost_cml = round(cumsum(net_policy_cost)),
        net_healthcare_cost_cml = round(cumsum(net_healthcare_cost)),
        net_socialcare_cost_cml = round(cumsum(net_socialcare_cost)),
        net_informal_care_cost_cml = round(cumsum(net_informal_care_cost)),
        net_productivity_cost_cml = round(cumsum(net_productivity_cost)),
        cpp_cvd_cml = round(cumsum(cpp_cvd)),
        cpp_chd_cml = round(cumsum(cpp_chd)),
        cpp_stroke_cml = round(cumsum(cpp_stroke)),
        cpp_poststroke_dementia_cml = round(cumsum(cpp_poststroke_dementia)),
        # cpp_htn_cml = round(cumsum(cpp_htn)),
        cpp_copd_cml = round(cumsum(cpp_copd)),
        # cpp_af_cml = round(cumsum(cpp_af)),
        cpp_t2dm_cml = round(cumsum(cpp_t2dm)),
        cpp_lung_ca_cml = round(cumsum(cpp_lung_ca)),
        cpp_colon_ca_cml = round(cumsum(cpp_colon_ca)),
        cpp_breast_ca_cml = round(cumsum(cpp_breast_ca)),

        cypp_cvd_cml = round(cumsum(cypp_cvd)),
        cypp_chd_cml = round(cumsum(cypp_chd)),
        cypp_stroke_cml = round(cumsum(cypp_stroke)),
        cypp_poststroke_dementia_cml = round(cumsum(cypp_poststroke_dementia)),
        # cypp_htn_cml = round(cumsum(cypp_htn)),
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
        total_hcp_cost_cml = round(cumsum(total_hcp_cost)),
        total_hscp_cost_cml = round(cumsum(total_hscp_cost)),
        societal_cost_cml = round(cumsum(societal_cost))
      ),
      by = c("friendly_name", "mc", "qimd")
    ][, nmb_cml := net_monetary_benefit_cml(.SD, input$health_econ_perspective_checkbox, input$out_wtp_box)
    ]
    # CBR cannot be summed because it is a ratio so it is inappropriate to put here
}

tt <- out_proc_qimd()
tt[, sum(cpp_chd), keyby = qimd]
tt[, sum(cpp_t2dm), keyby = qimd]

tt <- out_proc_qimd()[year == max(year), .(net_utility_cml, nmb_cml, eq5d_cml, pops_cml, mc, friendly_name, qimd)
]
# [, sum_dt(.SD, c("mc", "friendly_name", "qimd"), character(0))]
calc_rigit_scores(tt, c("mc", "friendly_name"))
calc_sei(tt, c("mc", "friendly_name"))
tt <- tt[, .(nmb_cml = sum(nmb_cml), sei = mean(sei)), by = c("mc", "friendly_name")]

