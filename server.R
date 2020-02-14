
server = function(input, output, session) {

  # restrict friendly names to 32 characters at max
  shinyjs::runjs("$('#friendly_name_sc1').attr('maxlength',32)")
  shinyjs::runjs("$('#friendly_name_sc2').attr('maxlength',32)")
  shinyjs::runjs("$('#friendly_name_sc3').attr('maxlength',32)")
  shinyjs::runjs("$('#friendly_name_sc4').attr('maxlength',32)")
  shinyjs::runjs("$('#friendly_name_sc5').attr('maxlength',32)")
  shinyjs::runjs("$('#friendly_name_sc6').attr('maxlength',32)")
  shinyjs::runjs("$('#friendly_name_sc7').attr('maxlength',32)")
  shinyjs::runjs("$('#friendly_name_sc8').attr('maxlength',32)")
  shinyjs::runjs("$('#friendly_name_sc9').attr('maxlength',32)")

# Dynamic scenario tab addition -------------------------------------------
  output$scenario_tabs <- renderUI({
    Tabs <- as.list(rep(0, input$scenarios_number_slider))
    for (i in 0:length(Tabs)) {
      Tabs[i] <- lapply(paste0("Scenario ", i), tabPanel, value = i)
    }

    #Tabs <- lapply(paste("Layer",0:input$subClust,sep=" "), tabPanel)
    do.call(tabsetPanel, c(Tabs, id = "level"))
  })


# Reactive function to choose scenario colours and shapes, and use their value according to their name

  colours <- reactive({
      setnames(transpose(as.data.table(lapply(seq_len(input$scenarios_number_slider), function(i) {
        c(input[[paste0("col_sc", i)]],
          input[[paste0("friendly_name_sc", i)]],
          input[[paste0("symbol_sc", i)]])
      }))), c("colour", "names", "symbol"))
  })



  # firstFrame <- vapplyseq_len(input$scenarios_number_slider), function(i) isTRUE(tr[["frame"]] %in% frameNames[[1]]), logical(1))
  # input$scenarios_number_slider[firstFrame] <- p$x$frames[[1]]$data

  #p$x$data <- c(c(p$x$data, colours), NULL)

  # callModule(module = createPlot, id = "cep1")




  observeEvent(input$scenarios_number_slider, {
# Render uptake/Px tables -------------------------------------------------
    for (i in seq_len(input$scenarios_number_slider)) {
      source(exprs = str2expression(gsub("_sc1",  paste0("_sc", i),
                                      readLines(file.path("server", "render_scenario_tables.R")))), local = TRUE
             )$value
    }

    # Alter ,min/max scenario init year -----------------------------------------
    for (i in seq_len(input$scenarios_number_slider)) {
      source(exprs = str2expression(gsub("_sc1",  paste0("_sc", i),
                                      readLines(file.path("server", "dynamic_scenario_minmax_year.R")))), local = TRUE
             )$value
    }

    # Load/save scenario  -----------------------------------------------------
    for (i in seq_len(input$scenarios_number_slider)) {
      source(exprs = str2expression(gsub("_sc1",  paste0("_sc", i),
                                      readLines(file.path("server", "save_load_scenarios.R")))), local = TRUE
             )$value
    }


    # Conditionally disable inputs -------
    for (i in seq_len(input$scenarios_number_slider)) {
      source(exprs = str2expression(gsub("_sc1",  paste0("_sc", i),
                                       readLines(file.path("server", "show_hide_logic.R")))), local = TRUE
      )$value
    }

    # Prev/Next buttons -------
    for (i in seq_len(input$scenarios_number_slider)) {
      source(exprs = str2expression(gsub("_sc1",  paste0("_sc", i),
                                       readLines(file.path("server", "prev_next_buttons.R")))), local = TRUE
      )$value
    }
  },
  ignoreNULL = FALSE, ignoreInit = FALSE)


  # Check if synthpop exists and generate if not ----------------------------
   observeEvent(input$locality_select_validator, {
     # See https://stackoverflow.com/questions/50165443/async-process-blocking-r-shiny-app
     # and https://github.com/rstudio/promises/issues/23
    init_parameters <- reactiveValuesToList(input)
    future({
       future_lapply(1:20, generate_synthpop,
        locality = init_parameters$locality_select,
        n = design$n, # number of simulants
        synthpop_dir,
        sim_horizon = design$sim_horizon_max,
        init_year = design$init_year_long, # unused but could be a vector
        max_lag = design$maxlag,
        smoking_relapse_limit = design$smoking_relapse_limit,
        ageL = design$ageL,
        ageH = design$ageH,
        jumpiness = design$jumpiness,
        simsmok_calibration = design$simsmok_calibration,
        include_diseases = TRUE,
        design = design,
        n_cpu =  design$n_cpus)
    })
    NULL # This is necessary!! after future when there is no promise to return
  })


# Run simulation ----
  out <- eventReactive(input[[paste0("run_simulation_sc",
                                     input$scenarios_number_slider)]], {

    parameters <- reactiveValuesToList(input)
    qsave(reactiveValuesToList(input), "./parameters.qs") # TODO delete for production

    # progress$inc(1/n, detail = paste("Doing part", i))
    withProgress(message = 'Running workHORSE model.',
                 detail = 'This may take a couple of minutes...', value = 0, {

                   # TODO remove before release
                   if (file.exists("./output/test_results.fst")) {

                     out <- read_fst("./output/test_results.fst", as.data.table = TRUE)
                     out

                   } else {

                     run_simulation(parameters, 1:design$iteration_n)

                   }
                 })
  })


# change to output panel
  observeEvent(input[[paste0("run_simulation_sc", input$scenarios_number_slider)]], {
    check_locality <- length(input$locality_select) > 0
    check_baseline <-
      any(
        input$baseline_sc1,
        input$baseline_sc2,
        input$baseline_sc3,
        input$baseline_sc4,
        input$baseline_sc5,
        input$baseline_sc6,
        input$baseline_sc7,
        input$baseline_sc8,
        input$baseline_sc9
      )
    # TODO What if more than one baseline scenario is selected

    if (!check_locality) {

      createAlert(
        session,
        paste0("alert_locality_sc", input$scenarios_number_slider),
        "alert_locality",
        title = "Oops... area is missing!",
        content = "Please select an area to simulate from the 'Simulation parameters' tab.",
        append = TRUE,
        style = "warning")

      } else closeAlert(session, "alert_locality")

     if (!check_baseline) {

      createAlert(
        session,
        paste0("alert_baseline_sc", input$scenarios_number_slider),
        "alert_baseline",
        title = "Oops...no baseline scenario!",
        content = "Please select a scenario as the baseline one. The baseline button is at top of each scenario page.",
        append =  TRUE,
        style = "warning")

    } else closeAlert(session, "alert_baseline")

      if (check_locality && check_baseline) {

        updateTabsetPanel(session, "inTabset",
          selected = "output_panel")

      }
  })

# Close alerts when error is corrected
  observeEvent(input$locality_select, {
    check_locality <- length(input$locality_select) > 0
    if (check_locality) closeAlert(session, "alert_locality")
  }
    )

  observeEvent({
    input$baseline_sc1
    input$baseline_sc2
    input$baseline_sc3
    input$baseline_sc4
    input$baseline_sc5
    input$baseline_sc6
    input$baseline_sc7
    input$baseline_sc8
    input$baseline_sc9
    1}, {

      check_baseline <-
        any(
          input$baseline_sc1,
          input$baseline_sc2,
          input$baseline_sc3,
          input$baseline_sc4,
          input$baseline_sc5,
          input$baseline_sc6,
          input$baseline_sc7,
          input$baseline_sc8,
          input$baseline_sc9
        )

      if (check_baseline) closeAlert(session, "alert_baseline")
    })



# Output tab --------------------------------------------------------
  output$most_cost_effective_box <- renderInfoBox({
    infoBox(
      "Most Cost-Effective", most_cost_effective(out_proc()), icon = icon("pound-sign"),
      color = "aqua", fill = FALSE
    )
  })

  output$most_equitable_box <- renderInfoBox({
    infoBox(
      "Most Equitable", most_equitable(out_proc_qimd()), icon = icon("balance-scale"),
      color = "purple"
    )
  })

  output$most_effective_box <- renderInfoBox({
    infoBox(
      "Most Effective", most_effective(out_proc()), icon = icon("heart"),
      color = "yellow"
    )
  })


  diff_year <- reactive({
    input$inout_year_slider - input$simulation_period_slider[1]
  })

  # Report button ----
  out_report <- eventReactive(input$produce_report, {
    parameters <- reactiveValuesToList(input)

    # progress$inc(1/n, detail = paste("Doing part", i))
    withProgress(message = 'Running workHORSE model.',
                 detail = 'This may take a couple of hours...',
                 value = 0,
                 {
                   run_simulation(parameters, 21:(design$iteration_n * 5L)) # 100

                 }
  )
})

# Automated graph explanation text  ------------------------------------------------------

  output$automated_text_descr <- renderUI({
    HTML(
      paste0(
        "With time horizon of ",
        diff_year(),
        " years, starting in year ",
        input$simulation_period_slider[1],
        ", the most cost effective scenario was ",
      most_benefit_cost_ratio(
        out_proc(),
        input$health_econ_perspective_checkbox,
        input$out_wtp_box
      ),
      " scenario. ",
      "It had a benefit:cost ratio (including the value of health gains) of £",
      most_benefit_cost_ratio_value(out_proc(),
                                    input$health_econ_perspective_checkbox,
                                    input$out_wtp_box),
      " for every £1 spent."
    , br(), br(),


      "This scenario had a societal incremental cost effectiveness ratio of ",
      round(scn_icer_cml(out_proc(), first(most_benefit_cost_ratio(out_proc(),
                                                            input$health_econ_perspective_checkbox,
                                                            input$out_wtp_box, 2)),
                   input$health_econ_perspective_checkbox) -
      scn_icer_cml(out_proc(), last(most_benefit_cost_ratio(out_proc(),
                                                            input$health_econ_perspective_checkbox,
                                                            input$out_wtp_box, 2)),
                   input$health_econ_perspective_checkbox)),
      " per QALY when compared to the next best scenario. ", br(), br(),


       "The order of the scenarios, starting from the most cost effective was ",
       head(most_benefit_cost_ratio(
         out_proc(),
         input$health_econ_perspective_checkbox,
         input$out_wtp_box, 10L
       ), -1),
      " and finally, ",
      last(most_benefit_cost_ratio(
        out_proc(),
        input$health_econ_perspective_checkbox,
        input$out_wtp_box, 10L
      )),
      "."
      ))
  })


  output$note_cost_eff <- renderUI({
    HTML(
      paste0(
        "With time horizon of ",
        diff_year(),
        " years, starting in year ",
        input$simulation_period_slider[1],
        ": and based on valuing QALYs at ",
        input$out_wtp_box,
        ", the scenario with the greatest societal benefit:cost ratio was ",
        most_benefit_cost_ratio(
          out_proc(),
          input$health_econ_perspective_checkbox,
          input$out_wtp_box
        ),
        " scenario which had a benefit: cost ratio of ",
        most_benefit_cost_ratio_value(
          out_proc(),
          input$health_econ_perspective_checkbox,
          input$out_wtp_box
        ),
        " for every £1 spent and a total net monetary benefit of £",
        most_benefit_cost_ratio_value(
          out_proc(),
          input$health_econ_perspective_checkbox,
          input$out_wtp_box
        ),
        " from a societal perspective.",
        br(),
        br(),
        "This scenario had a societal incremental cost effectiveness ratio of ",
        round(scn_icer_cml(out_proc(), first(most_benefit_cost_ratio(out_proc(),
                                                                     input$health_econ_perspective_checkbox,
                                                                     input$out_wtp_box, 2)),
                           input$health_econ_perspective_checkbox) -
                scn_icer_cml(out_proc(), last(most_benefit_cost_ratio(out_proc(),
                                                                      input$health_econ_perspective_checkbox,
                                                                      input$out_wtp_box, 2)),
                             input$health_econ_perspective_checkbox)),
        " per QALY when compared to the next best scenario ",
        last(most_benefit_cost_ratio(out_proc(),
                                     input$health_econ_perspective_checkbox,
                                     input$out_wtp_box, 2)),
        " scenario.",
        br(),
        br(),
        # "The most cost effective scenario from a healthcare cost perspective was ",
        # most_cost_effective(out_proc()),
        # " scenario which had an incremental cost effectiveness ratio of ",
        # incr_cost_effect_ratio(),
        # " per QALY gained.",
        # br(),
        "Ranking the scenarios from most to least cost effective, the order was:",
        br(),
        br(),
        most_benefit_cost_ratio(
          out_proc(),
          input$health_econ_perspective_checkbox,
          input$out_wtp_box, 10L
        )
        # paste(head(rank_cost_effective(out_proc(
        # )), -1), collapse = ", "),
        # last_rank()
      )
    )
  })

   #TODO: update with new functions for absolute and relative equitable

   output$note_health_ineq <- renderUI({
     HTML(paste0(most_equitable_rel(out_proc_qimd()), " scenario had the biggest impact in terms of reducing relative inequalities, reducing the relative index of inequalities by ", reduce_rel_index_ineq(), ".", br(), br(),
                 last(most_equitable_abs(out_proc_qimd())), " scenario did the least harm to absolute inequalities, increasing the absolute index of inequalities by ", increase_abs_index_ineq(), "."))

   })


   output$note_effvnss <- renderUI({
     HTML(paste0("The most effective scenario in terms of disease cases prevented and postponed was ", most_effective(out_proc()), " scenario which prevented '...' total disease cases,
                 made up of '...' cases of CVD, and '...' cases of other diseases"))
   })


   output$note_social_benef <- renderUI({
     HTML(paste0("As well as healthcare benefits, ", most_equitable(out_proc())," scenario would produce additional social care cost savings of ",
                 social_care_cost_sav(), ", productivity benefits of ", prod_benef(), " (including earnings
and the value of household productivity) and informal care cost savings of", inform_care_cost_sav(), br()))
   })

 output$info_ce_plane <- renderText({
   "Needs to be implemented"
 })

 output$info_equ_plane <- renderText({
   "The future explications for this specific graph"
 })

 output$info_ce1_plane <- renderText({
   "The future explications for this specific graph"
 })

 output$info_cepcs_plane <- renderText({
   "The future explications for this specific graph"
 })

 output$info_cepce_plane <- renderText({
   "The future explications for this specific graph"
 })

 output$info_cep_anim_plane <- renderText({
   "The future explications for this specific graph"
 })

 output$info_equ_rel_plane <- renderText({
   "The future explications for this specific graph"
 })

 output$info_equ1_1_plane <- renderText({
   "The future explications for this specific graph"
 })

 output$info_equ_p_rel_plane <- renderText({
   "The future explications for this specific graph"
 })

 output$info_equ_p_abs_plane <- renderText({
   "The future explications for this specific graph"
 })

output$info_equ_anim_rel_plane <- renderText({
  "The future explications for this specific graph"
})

output$info_equ_anim_abs_plane <- renderText({
  "The future explications for this specific graph"
})

 # user_inputs <- reactive(reactiveValuesToList(input))


  source(file.path("server", "output_inputs.R"),  local = TRUE)$value


# Tooltips ----------------------------------------------------------------
  # source(file.path("server", "tooltips.R"),  local = TRUE)$value


  # outputOptions(output, "nrows", suspendWhenHidden = FALSE)


# last_rank <- function(){
#   lr <- rank_cost_effective(out_proc())[length(rank_cost_effective(out_proc()))]
#   lr <- paste0(" and ", last(rank_cost_effective(out_proc())))
# }

}

