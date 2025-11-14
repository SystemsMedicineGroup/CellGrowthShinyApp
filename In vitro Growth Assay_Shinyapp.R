
######################################################################################################################################
## Latest Shiny App for the cell growth and viability model
######################################################################################################################################
library(shiny)
library(ggplot2)
library(FME)
library(deSolve)
library(dplyr)
library(readr)
library(tidyr)
library(gridExtra)
library(stringr)
library(tibble)

# Load Fortran DLLs only once
setwd("C:/Users/ibanda/OneDrive - The University of Liverpool/Attachments/TRR1 Project/Programming and Data Analysis UoL/Fortran Magic/Trr1 Model - Both models")
if (!"modnt" %in% names(getLoadedDLLs())) dyn.load("modnt.dll")
if (!"modwt" %in% names(getLoadedDLLs())) dyn.load("modwt.dll")


ui <- fluidPage(
  titlePanel("Generic Modelling Shiny App for Trr1 In Vitro Data"),
  
  sidebarLayout(
    sidebarPanel(
      radioButtons("datatype", "Select Data Type", 
                   choices = c("Optical Density (OD)", "Viability (FACS)")),
      
      fileInput("od_file", "Upload OD Data CSV"),
      
      conditionalPanel(
        condition = "input.datatype == 'Optical Density (OD)'",
        selectInput("exp_id", 
                        "Select Experimental ID", choices = NULL),
        selectInput("strain", 
                        "Select Strain", choices = NULL),
        selectInput("OD0", 
                        "Select Initial OD", choices = NULL),
        actionButton("plot_OD", 
                        "Plot OD Data", class = "btn-success")
      ),
      br(),
      fileInput("viability_file", 
                    "Upload Viability Data CSV"),
        conditionalPanel(
        condition = "input.datatype == 'Viability (FACS)'",
        selectInput("dox", 
                        "Select Doxycycline conc. (μg/mL)", choices = NULL),
        actionButton("plot_viability", 
                        "Plot Viability Data", class = "btn-success")
      ),
      br(),
      selectInput("model_choice", "Choose Model", 
                  choices = c("Model 1 (No transition rate)", "Model 2 (With transition rate)")),
      
      uiOutput("param_inputs"),
                  verbatimTextOutput("modelCode"),
      
      sliderInput("num_iter", "Number of MCMC iterations:",
                  min = 500, max = 5000,
                  value = 500,
                  step = 500),
      
      actionButton("run_model", "Run Model", class = "btn-primary"),
      tags$div(
        style = "margin-top: 5px; color: green; font-size: 13px;",
        p("Ensure that both data sets are", 
                        strong("Plotted"), "and the OD type is selected before clicking", 
                        strong("Run Model"))
      ),
      br(), br(),
      downloadButton("download_fitted_params", "Download Fitted Parameters")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("OD Data Visualisation", 
                 conditionalPanel(
                   condition = "input.datatype == 'Optical Density (OD)'",
                   plotOutput("odPlot", height = "600px", width = "100%")
                 )
        ),
        tabPanel("Viability Data Visualisation", 
                 conditionalPanel(
                   condition = "input.datatype == 'Viability (FACS)'",
                   plotOutput("viabilityPlot", height = "600px", width = "100%")
                 )
        ),
        tabPanel("Model Code", 
                    verbatimTextOutput("modelCode")),
        tabPanel("Fit Results", 
                    verbatimTextOutput("fitResults")),
        tabPanel("Growth rate versus DOX", 
                    plotOutput("growthratePlot", height = "600px", width = "100%")),
        tabPanel("Decay rate versus DOX", 
                    plotOutput("decayratePlot", height = "600px", width = "100%")),
        tabPanel("Transition rate versus DOX", 
                 plotOutput("transratePlot", height = "600px", width = "100%")),
        tabPanel("Overlay Plot", 
                    plotOutput("overlayPlot", height = "600px", width = "100%")),
        tabPanel("Predicted Dead Cell Plot", 
                    plotOutput("deadcellPlot", height = "600px", width = "100%")),
        tabPanel("MCMC Results", 
                     selectInput("plot_choice", "Select conc:", choices = NULL), 
                     plotOutput("mcmcPlot", height = "600px", width = "100%")),
        tabPanel("Sensitivity Plot", 
                     selectInput("sensplot_choice", "Select conc:", choices = NULL),
                     plotOutput("sensPlot", height = "600px", width = "100%")),
        tabPanel("Collinearity Index", 
                     selectInput("collinplot_choice", "Select conc:", choices = NULL),
                     plotOutput("collinPlot", height = "600px", width = "100%")),
        tabPanel("VPC (SensRange) Plot", 
                     plotOutput("sensRangePlot", height = "600px", width = "100%"))
      )
    )
  )
)

server <- function(input, output, session) {
  
  # Reactive to read OD data
  od_data <- reactive({
    req(input$od_file)
    read.csv(input$od_file$datapath)
  })
  
  observeEvent(input$od_file, {
    req(input$datatype == "Optical Density (OD)") #
    od_df <- od_data()
    
    # --- Detect experiment ID column ---
    possible_cols <- c("exp_id", "Experiment_ID", "ExperimentalID", "ID")
    exp_col <- possible_cols[possible_cols %in% names(od_df)][1]
    
    # Update exp_id selectInput
    if (!is.null(exp_col)) {
      updateSelectInput(session, "exp_id", choices = sort(unique(od_df[[exp_col]])))
    } else {
      updateSelectInput(session, "exp_id", choices = NULL)
    }
    
    # --- Update strain choices ---
    if ("strain" %in% names(od_df)) {
      updateSelectInput(session, "strain", choices = sort(unique(od_df$strain)))
    } else {
      updateSelectInput(session, "strain", choices = NULL)
    }
    
    # --- Update OD0 choices ---
    if ("OD0" %in% names(od_df)) {
      updateSelectInput(session, "OD0", choices = sort(unique(od_df$OD0)))
    } else {
      updateSelectInput(session, "OD0", choices = NULL)
    }
  })
  
  filtered_od_data <- reactive({
    req(input$datatype == "Optical Density (OD)", input$exp_id, input$strain, input$OD0, od_data())
    od_df <- od_data()
    
    # Detect experiment ID column
    possible_cols <- c("exp_id", "Experiment_ID", "ExperimentalID", "ID")
    exp_col <- possible_cols[possible_cols %in% names(od_df)][1]
    
    validate(
      need(!is.null(exp_col), "No valid experiment ID column found in OD data.")
    )
    
    # Filter by selected values
    od_df <- od_df %>%
      dplyr::filter(
        !!rlang::sym(exp_col) == isolate(input$exp_id),
        strain == isolate(input$strain),
        OD0 == isolate(input$OD0)
      )
    
    return(od_df)
  })
  
  # Reactive expression to prepare data only after action button is clicked
  od_plot_data <- eventReactive(input$plot_OD, {
    od_df <- filtered_od_data()
    
    validate(
      need(nrow(od_df) > 0, "No data available for the selected filters."),
      need(all(c("time", "OD", "dox") %in% names(od_df)), 
           "OD data must contain 'time', 'OD', and 'dox' columns.")
    )
    
    od_df <- od_df %>%
      dplyr::group_by(time, dox) %>%
      dplyr::summarise(
        mean_OD = mean(OD, na.rm = TRUE),
        sd_OD = sd(OD, na.rm = TRUE),
        n = dplyr::n(),
        .groups = "drop"
      ) %>%
      dplyr::rename(OD = mean_OD, sd = sd_OD) %>%
      dplyr::mutate(sd = tidyr::replace_na(sd, 0)) %>%
      dplyr::filter(OD > 0)
    
    return(od_df)
  })
  
  # Output plot only when button is clicked and data is ready
  output$odPlot <- renderPlot({
    req(input$plot_OD > 0) 
    od_df <- od_plot_data()  # Triggers only on button click
    
    ggplot(od_df, aes(x = time, y = OD)) +
      geom_point(size = 1.5) +
      geom_errorbar(aes(ymin = OD - sd, ymax = OD + sd),
                    width = 0.2, alpha = 0.7) +
      labs(
        title = paste("Strain:", input$strain, "OD0:", input$OD0),
        x = "Time (hours)",
        y = expression("OD"[600] ~ "(log"[10]*")")
      ) +
      facet_wrap(~dox, labeller = label_both, ncol = 6) +
      scale_x_continuous(
        limits = c(min(od_df$time), max(od_df$time)),
        breaks = seq(min(od_df$time), max(od_df$time), 3)
        ) +
      scale_y_log10(
        limits = c(0.001, 2),
        minor_breaks = c(seq(1, 9, 1)*0.001, seq(1, 9, 1) * 0.01,
                         seq(1, 9, 1) * 0.1, seq(1, 9) * 1, 2)
      ) +
      theme_minimal() +
      theme(
        title = element_text(size = 18),
        strip.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.text.x = element_text(size = 14)
      )
  })
  
  # Reactive FACS data
  viability_data <- reactive({
    req(input$viability_file)
    read.csv(input$viability_file$datapath)
  })
  
  # Dynamically update dox dropdown based on uploaded file
  observeEvent(input$viability_file, {
    req(input$viability_file)
    df <- viability_data()
    
    if (!"dox" %in% names(df)) {
      showNotification("Column 'dox' not found in uploaded file.", type = "error")
      updateSelectInput(session, "dox", choices = NULL)
      return()
    }
    
    # Convert dox values
    dox_values <- suppressWarnings(as.numeric(as.character(df$dox)))
    dox_values <- dox_values[!is.na(dox_values)]
    
    if (length(dox_values) == 0) {
      showNotification("No valid numeric values found in 'dox' column.", type = "error")
      updateSelectInput(session, "dox", choices = NULL)
      return()
    }
    
    updateSelectInput(session, "dox", choices = sort(unique(dox_values)))
  })
  
  # Triggered data when plotting
  viability_data_triggered <- eventReactive(input$plot_viability, {
    df <- viability_data()
    validate(need(all(c("time", "dox", "dead_cell_prop") %in% names(df)), 
                  "Viability data must have 'time', 'dox' and 'dead_cell_prop' columns."))
    
    # Ensuring that dox is numeric (just in case)
    df$dox <- isolate(as.numeric(as.character(df$dox)))
    
    # Filter by the selected dox value and isolating it
    selected_dox <- isolate(input$dox)
    df <- df %>%
      dplyr::filter(dox == as.numeric(selected_dox)) %>%
      dplyr::rename(dead = dead_cell_prop) %>%
      dplyr::select(time, dead, dox)
    
    df
  })
  
  # Plot
  output$viabilityPlot <- renderPlot({
    df <- viability_data_triggered()
    
    validate(
      need(nrow(df) > 0, paste("No data found for selected dox =", input$dox))
    )
    
    ggplot(df, aes(x = as.factor(time), y = dead)) +
      geom_bar(stat = "identity", fill = "black", width = 0.25) +
      labs(
        title = paste("Viability at DOX =", input$dox),
        x = "Time (hours)",
        y = "Proportion of Cells Dead"
      ) +
      scale_x_discrete(
        labels = c("0" = "0h", "5" = "5h", "10" = "10h", "15" = "15h", "20" = "20h")
      ) +
      scale_y_continuous(
        limits = c(0, 1),
        breaks = seq(0, 1, by = 0.25)
      ) +
      theme_minimal() +
      theme(
        title = element_text(size = 18),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.text.x = element_text(size = 14)
      )
  })
  
  # Available models
  # --- ODE Models ---
  model1 <- function(time, state, pars) {
    with(as.list(c(state, pars)), {
      P <- state["P"] # Proliferating cells
      D <- state["D"] # Dead cells
      
      dP = r * P * (1 - P / K) - kd * P
      dD = kd * P
      
      return(list(c(dP, dD)))
    })
  }
  
  model2 <- function(time, state, pars) {
    with(as.list(c(state, pars)), {
      Q <- state["Q"] # Dormant cells
      P <- state["P"] # Proliferating cells
      D <- state["D"] # Dead cells
      
      dQ = -kt * Q
      dP = kt * Q + r * P * (1 - P / K) - kd * P
      dD = kd * P
      
      return(list(c(dQ, dP, dD)))
    })
  }
  
  ## Show model code in UI
  output$modelCode <- renderText({
    req(input$model_choice)
    
    if (input$model_choice == "Model 1 (No transition rate)") {
      model1_code <- paste(capture.output(model1), collapse = "\n")
      paste("Model 1 function without transition:\n", model1_code)
      
    } else if (input$model_choice == "Model 2 (With transition rate)") {
      model2_code <- paste(capture.output(model2), collapse = "\n")
      paste("Model 2 function with transition:\n", model2_code)
    } else {
      "Select a model to see its function code."
    }
  })
  
  # --- Dynamic parameter UI ---
  output$param_inputs <- renderUI({
    req(input$model_choice)
    
    if (input$model_choice == "Model 1 (No transition rate)") {
      tagList(
        numericInput("r", "Growth rate: r", value = 0.4),
        numericInput("K", "Carrying capacity: K", value = 1),
        numericInput("kd", "Decay rate: kd", value = 0.01)
      )
    } else if (input$model_choice == "Model 2 (With transition rate)") {
      tagList(
        numericInput("r", "Growth rate: r", value = 0.4),
        numericInput("K", "Carrying capacity: K", value = 1),
        numericInput("kd", "Decay rate: kd", value = 0.01),
        numericInput("kt", "Transition rate: kt", value = 0.3)
      )
    }
  })
  
  # --- Number of iterations for MCMC simulations--- 
  niter <- reactive({
    req(input$num_iter)
    niter <- as.numeric(input$num_iter)
      
    })
  
  # --- Simulation function ---
  sim_model <- function(OD0, pars, model_type) {
    times <- seq(0, 24, by = 0.5)
    
    if (model_type == "model1") {
      dllname <- "modnt"
      initcfunc <- "initcal1"
      func_name <- "derivscal1"  
      state <- c(P = OD0, D = 0)
     
    } else if (model_type == "model2") {
      dllname <- "modwt"
      func_name <- "derivscal2"  
      initcfunc <- "initcal2"
      state <- c(Q = OD0, P = 0, D = 0)
    } else {
      stop("Invalid model_type")
    }
    
    do.call(".Fortran", list(initcfunc, as.double(pars)))
    
  # The solution
    out <- tryCatch({
      ode(
        y = state,
        times = times,
        func = func_name,
        parms = NULL,
        dllname = dllname
      )
    }, error = function(e) {
      message("ODE solver failed: ", e$message)
      return(NULL)  
    })
    
    # If ODE failed, return NULL or empty data frame early
    if (is.null(out)) {
      return(NULL)
    }
    
    out <- as.data.frame(out)
          
    out$Q <- if (!"Q" %in% names(out)) 0 else out$Q
    out$OD <- out$Q + out$P + out$D
    out$dead <- out$D / out$OD
    
    return(out)
  }
  
  # --- Model fitting ---
  model_run <- function(od_df, df, model_type, pars, niter = 500) {
    od_mat <- as.matrix(od_df[, c("time", "OD")])
    death_mat <- as.matrix(df[, c("time", "dead")])
    
    OD_0 <- od_df %>%
      arrange(time) %>% 
      slice(1) %>%
      pull(OD)
    
    cost_fn <- function(pars) {
      out <- sim_model(OD0 = OD_0, pars = pars, model_type = model_type)
      cost <- FME::modCost(model = out, obs = od_mat)
      cost <- FME::modCost(model = out, obs = death_mat, cost = cost)
      return(cost)
    }
    
    Sfun <- sensFun(func = cost_fn, parms = pars)
    sens_func <- Sfun
    
   ident <- collin(Sfun)
    
    cost_fn_log <- function(log_pars) {
      pars <- exp(log_pars)
      cost_fn(pars)
    }
    
       # Model Fit
    fit <- FME::modFit(
      f = cost_fn_log,
      p = log(pars*2)
          )
    
    # Fit results
    fit_summary <- summary(fit)
    log_pars <- coef(fit)
    se_log <- sqrt(diag(fit_summary$cov.scaled))
    
    lower_log <- log_pars - 1.96 * se_log
    upper_log <- log_pars + 1.96 * se_log
    
    ci_table <- data.frame(
            parameter = names(log_pars),
            fitted_value = exp(log_pars),
            lower_95_CI = exp(lower_log),
            upper_95_CI = exp(upper_log)
    )
    
    final_model <- sim_model(OD0 = OD_0, pars = exp(log_pars), model_type = model_type)
    
    obs_pred_data <- data.frame(od_mat) %>%
      left_join(final_model, by = "time", suffix = c("", "_pred"))
    
    set.seed(42)
    mcmc <- FME::modMCMC(
            f = cost_fn_log,
            p = log_pars,
            jump = fit_summary$cov.scaled * 2.4^2 / length(log_pars),
            var0 = fit$var_ms_unweighted,
            niter = niter,
            burninlength = 0.1*niter,
            updatecov = 100
    )
      mcmc <- mcmc
      mcmc_samples = as.data.frame(exp(mcmc$pars))
      
    ## Sense Range determination
    sR = sensRange(
      func = function(pars) sim_model(OD0 = OD_0, pars = pars, model_type = model_type),
      parms = exp(fit$par),
      parInput = exp(mcmc$par))
    
    sR <- summary(sR)
    
    list(
      sens_func = sens_func,
      collinearity = ident,
      ci_table = ci_table,
      obs_pred_data = obs_pred_data,
      mcmc = mcmc,
      mcmc_samples = mcmc_samples,
      sR = sR
    )
  }
  
  # --- Reactive for model parameters ---
  model_params <- reactive({
    req(input$r, input$K, input$kd)
    
    if (input$model_choice == "Model 1 (No transition rate)") {
      c(r = input$r, K = input$K, kd = input$kd)
    } else if (input$model_choice == "Model 2 (With transition rate)") {
      req(input$kt)
      c(r = input$r, K = input$K, kd = input$kd, kt = input$kt)
    } else {
      stop("Unknown model choice.")
    }
  })
  
  # --- Store current run result ---
  fit_results <- reactiveVal(NULL)
  
  observeEvent(input$run_model, {
    req(filtered_od_data(), viability_data_triggered())
    
    od_df <- filtered_od_data() %>%
            dplyr::group_by(time, dox) %>%
            dplyr::summarise(
                      mean_OD = mean(OD, na.rm = TRUE),
                      sd_OD = sd(OD, na.rm = TRUE),
                      n = dplyr::n(),
                  .groups = "drop"
      ) %>%
            dplyr::rename(OD = mean_OD, sd = sd_OD) %>%
            dplyr::mutate(sd = tidyr::replace_na(sd, 0)) %>%
            dplyr::filter(OD > 0)
    
    df <- viability_data_triggered() 
    
    dox_values <- sort(unique(od_df$dox))
    results_list <- list()
    
    withProgress(message = "Running models...", value = 0, {
      n <- length(dox_values)
      for (i in seq_along(dox_values)) {
        concn <- dox_values[i]
        data_subset_od <- subset(od_df, dox == concn)
       
        model_type <- if (input$model_choice == "Model 1 (No transition rate)") "model1" else "model2"
        
      result <- tryCatch({
          model_run(
            od_df = data_subset_od,
            df = df,
            model_type = model_type,
            pars = model_params()
          )
        }, error = function(e) {
          showNotification(paste("Model failed at conc", concn, ":", e$message), type = "error")
          return(NULL)
        })
        
        print("OD Data used in the model:")
                print(data_subset_od)
        
        print("Actual Initial OD used:")
        print(data_subset_od %>%
                slice(1) %>%
                pull(OD))
        
        print("Viability Data used in the model:")
        print(df)
        
        print("Running model with initial parameters:")
        print(model_params())
        
        print("Number of MCMC iterations in use:")
        print(as.numeric(input$num_iter))
        
        if (!is.null(result)) {
          results_list[[paste0("dox_", concn)]] <- result
        }
        
        incProgress(1/n, detail = paste("Completed conc:", concn))
      }
    })
    
    fit_results(results_list)
    
    showNotification("Model run complete!", type = "message")
  })
  
  
  ## Parameter estimates after model fitting
  output$fitResults <- renderText({
    results <- fit_results()
    req(results)
    
    ci_all <- do.call(rbind, lapply(names(results), function(name) {
      if (!is.null(results[[name]]) && !is.null(results[[name]]$ci_table)) {
        cbind(concn = name, results[[name]]$ci_table)
      } else {
        data.frame(
          conc = name,
          parameter = NA, 
          fitted_value = NA,
          lower_95_CI = NA, 
          upper_95_CI = NA
        )
      }
    }))
    
    paste(capture.output(print(ci_all)), collapse = "\n")
  })
  
  ## Plot for growth rate (r) against DOX
  output$growthratePlot <- renderPlot({
    results <- fit_results()
    req(results, input$exp_id)
    
    mcmc_all <- do.call(rbind, lapply(names(results), function(name) {
      if (!is.null(results[[name]]) && !is.null(results[[name]]$mcmc_samples)) {
        cbind(concn = name, results[[name]]$mcmc_samples)
      } else {
        NULL  # Skip if missing
      }
    }))
    
    if(input$exp_id == 5) {
    sum_r <- bind_rows(mcmc_all) %>%
      mutate(dox = as.numeric(sub("dox_", "", concn)),
             dox_adj = ifelse(dox == 0, 3e-3, dox)) %>%
      select(dox, r, dox_adj) %>%
      group_by(dox, dox_adj) %>%
       summarise(
        median = median(r),
        q025 = quantile(r, 0.025),
        q975 = quantile(r, 0.975),
        .groups = "drop"
      )
    } else {
      sum_r <- bind_rows(mcmc_all) %>%
        mutate(dox = as.numeric(sub("dox_", "", concn)),
               dox_adj = ifelse(dox == 0, 3e-1, dox)) %>%
        select(dox, r, dox_adj) %>%
        group_by(dox, dox_adj) %>%
        summarise(
          median = median(r),
          q025 = quantile(r, 0.025),
          q975 = quantile(r, 0.975),
          .groups = "drop"
        )
    }
            # Calculate range for breaks
      min_dox <- min(sum_r$dox_adj, na.rm = TRUE)
      max_dox <- max(sum_r$dox_adj, na.rm = TRUE)
      min_par <- round(min(sum_r$q025, na.rm = TRUE), 1)
      max_par <- round(max(sum_r$q975, na.rm = TRUE), 1)
      
      log_min <- floor(log10(min_dox))
      log_max <- ceiling(log10(max_dox))
      
      # Major: powers of 10; Minor: 1–9 × those powers
      major_breaks <- 10^(log_min:log_max)
      minor_breaks <- unlist(lapply(major_breaks, function(b) b * (1:9)))
      minor_breaks <- minor_breaks[minor_breaks >= min_dox & minor_breaks <= max_dox]
      
  ggplot(sum_r, aes(x = dox_adj, y = median)) +
        geom_point(size = 1.5) +
        geom_errorbar(aes(ymin = q025, ymax = q975), 
                      width = 0.025, lwd = 1) +
        scale_x_log10(
          name = "Doxycycline Concentration (µg/mL)",
          limits = c(min_dox*0.7, max_dox*1.5),
          breaks = major_breaks,
          minor_breaks = minor_breaks
        ) +
        scale_y_continuous(
          name = "Estimated Growth Rate (Median & 95% CI)",
          limits = c(min_par*0.8, max_par*1.25),
          breaks = seq(min_par, max_par*1.25, 0.1)
        ) +
        labs(title = paste("Effect of DOX on the growth rate of", 
                           input$OD0,
                           input$strain, "in Experiment", input$exp_id),
             ) +
        theme_minimal() +
        theme(
          axis.title = element_text(size = 16),
          axis.text = element_text(size = 16),
          axis.text.x = element_text(size = 14, vjust = 0.5),
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
        ) +
    geom_vline(
      xintercept = min_dox,
      linetype = "dashed",
      colour = "black",
      lwd = 1.2
    ) +
    annotate(
      "text",
      x = 1.1*min_dox,          
      y = max_par - 0.1,    
      label = "DOX = 0",
      angle = 90,           
      vjust = 0.5,          
      hjust = 1,         
      size = 4
    )
  })
  
  ## Relationship between the decay rate and DOX
    output$decayratePlot <- renderPlot({
    results <- fit_results()
    req(results, input$exp_id)
    
    mcmc_all <- do.call(rbind, lapply(names(results), function(name) {
      if (!is.null(results[[name]]) && !is.null(results[[name]]$mcmc_samples)) {
        cbind(concn = name, results[[name]]$mcmc_samples)
      } else {
        NULL  # Skip if missing
      }
    }))
    
    if(input$exp_id == 5) {
    sum_kd <- bind_rows(mcmc_all) %>%
      mutate(dox = as.numeric(sub("dox_", "", concn)),
             dox_adj = ifelse(dox == 0, 3e-3, dox)) %>%
      select(dox, kd, dox_adj) %>%
      group_by(dox, dox_adj) %>%
      summarise(
        median = median(kd),
        q025 = quantile(kd, 0.025),
        q975 = quantile(kd, 0.975),
      .groups = "drop")
    
    } else {
    sum_kd <- bind_rows(mcmc_all) %>%
          mutate(dox = as.numeric(sub("dox_", "", concn)),
                 dox_adj = ifelse(dox == 0, 6e-1, dox)) %>%
          select(dox, kd, dox_adj) %>%
          group_by(dox, dox_adj) %>%
          summarise(
            median = median(kd),
            q025 = quantile(kd, 0.025),
            q975 = quantile(kd, 0.975),
          .groups = "drop")
      
    }
    
    # Calculate range for breaks
    min_dox <- min(sum_kd$dox_adj, na.rm = TRUE)
    max_dox <- max(sum_kd$dox_adj, na.rm = TRUE)
    min_par <- signif(min(sum_kd$q025, na.rm = TRUE), 1)
    max_par <- signif(max(sum_kd$q975, na.rm = TRUE), 1)
    int_par <- signif(((max_par - min_par)/5), 1)
    
    log_min <- floor(log10(min_dox))
    log_max <- ceiling(log10(max_dox))
    
    # Major: powers of 10; Minor: 1–9 × those powers
    major_breaks <- 10^(log_min:log_max)
    minor_breaks <- unlist(lapply(major_breaks, function(b) b * (1:9)))
    minor_breaks <- minor_breaks[minor_breaks >= min_dox & minor_breaks <= max_dox]
    
    ggplot(sum_kd, aes(x = dox_adj, y = median)) +
      geom_point(size = 1.5) +
      geom_errorbar(aes(ymin = q025, ymax = q975), 
                    width = 0.025, lwd = 1) +
      scale_x_log10(
        name = "Doxycycline Concentration (µg/mL)",
        limits = c(min_dox*0.9, max_dox*1.5),
        breaks = major_breaks,
        minor_breaks = minor_breaks
      ) +
      scale_y_continuous(
        name = "Estimated Decay Rate (Median & 95% CI)",
        limits = c(0, max_par*1.25),
        breaks = seq(0, max_par*1.25, int_par)
      ) +
      labs(title = paste("Effect of DOX on the decay rate of", 
                         input$OD0,
                         input$strain, "in Experiment", input$exp_id),
      ) +
      theme_minimal() +
      theme(
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.text.x = element_text(size = 14, hjust = 0.5),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
      ) +
      geom_vline(
        xintercept = min_dox,
        linetype = "dashed",
        colour = "black",
        lwd = 1.2
      ) +
      annotate(
        "text",
        x = 1.15*min_dox,          
        y = max_par,    
        label = "DOX = 0",
        angle = 90,           
        vjust = 0.5,          
        hjust = 1,         
        size = 4
      )
  })
  
    ## Transition rate verses DOX
    output$transratePlot <- renderPlot({
      results <- fit_results()
      req(results, input$exp_id, input$model_choice)
        
        if (input$model_choice == "Model 2 (With transition rate)") {
          mcmc_all <- do.call(rbind, lapply(names(results), function(name) {
            if (!is.null(results[[name]]) && !is.null(results[[name]]$mcmc_samples)) {
              cbind(concn = name, results[[name]]$mcmc_samples)
            } else {
              NULL
            }
          }))
          
          if(input$exp_id == 5) {
            sum_kt <- bind_rows(mcmc_all) %>%
              mutate(dox = as.numeric(sub("dox_", "", concn)),
                     dox_adj = ifelse(dox == 0, 3e-3, dox)) %>%
              select(dox, kt, dox_adj) %>%
              group_by(dox, dox_adj) %>%
              summarise(
                median = median(kt),
                q025 = quantile(kt, 0.025),
                q975 = quantile(kt, 0.975),
                .groups = "drop")
          } else {
            sum_kt <- bind_rows(mcmc_all) %>%
              mutate(dox = as.numeric(sub("dox_", "", concn)),
                     dox_adj = ifelse(dox == 0, 6e-1, dox)) %>%
              select(dox, kt, dox_adj) %>%
              group_by(dox, dox_adj) %>%
              summarise(
                median = median(kt),
                q025 = quantile(kt, 0.025),
                q975 = quantile(kt, 0.975),
                .groups = "drop") %>%
              filter(q975 < 10)
          }
          
          # Calculate range for breaks
          min_dox <- min(sum_kt$dox_adj, na.rm = TRUE)
          max_dox <- max(sum_kt$dox_adj, na.rm = TRUE)
          min_par <- signif(min(sum_kt$q025, na.rm = TRUE), 1)
          max_par <- signif(max(sum_kt$q975, na.rm = TRUE), 1)
          int_par <- signif(((max_par - min_par)/5), 1)
         
          log_min <- floor(log10(min_dox))
          log_max <- ceiling(log10(max_dox))
          
          major_breaks <- 10^(log_min:log_max)
          minor_breaks <- unlist(lapply(major_breaks, function(b) b * (1:9)))
          minor_breaks <- minor_breaks[minor_breaks >= min_dox & minor_breaks <= max_dox]
          
          ggplot(sum_kt, aes(x = dox_adj, y = median)) +
            geom_point(size = 1.5) +
            geom_errorbar(aes(ymin = q025, ymax = q975), 
                          width = 0.025, lwd = 1) +
            scale_x_log10(
              name = "Doxycycline Concentration (µg/mL)",
              limits = c(min_dox*0.9, max_dox*1.15),
              breaks = major_breaks,
              minor_breaks = minor_breaks
            ) +
            scale_y_continuous(
              name = "Estimated Transition Rate (Median & 95% CI)",
              limits = c(0, max_par*1.5),
              breaks = seq(0, max_par*1.5, int_par)
            ) +
            labs(title = paste("Effect of DOX on the transition rate of", 
                               input$OD0,
                               input$strain, "in Experiment", input$exp_id)
            ) +
            theme_minimal() +
            theme(
              axis.title = element_text(size = 16, vjust = 1),
              axis.text = element_text(size = 16),
              axis.text.x = element_text(size = 14, hjust = 0.5),
              plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
            ) +
            geom_vline(
              xintercept = min_dox,
              linetype = "dashed",
              colour = "black",
              lwd = 1.2
            ) +
            annotate(
              "text",
              x = 1.15*min_dox,          
              y = max_par,    
              label = "DOX = 0",
              angle = 90,           
              vjust = 0.5,          
              hjust = 1,         
              size = 4
            )
          
        } else {
          # For other models, show a message on the plot
          plot.new()
          text(0.5, 0.5, paste0(input$model_choice, " does not estimate transition rate"), cex = 1.5)
        }
      })
    
  ## Overlay of the observed by the predicted
  output$overlayPlot <- renderPlot({
    results <- fit_results()
    req(results)
    
    # Combine all obs_pred_data with the dose/concentration name
    obs_pred_data_all <- do.call(rbind, lapply(names(results), function(name) {
      res_item <- results[[name]]
      if (!is.null(res_item$obs_pred_data)) {
        cbind(concn = name, res_item$obs_pred_data)
      } else {
        NULL  
      }
    }))
    
    # Require valid data
    req(nrow(obs_pred_data_all) > 0)
    
    # Pivot to long format for easier ggplotting
    obs_pred_data <- tidyr::pivot_longer(
      obs_pred_data_all,
      cols = c(OD, OD_pred),
      names_to = "Data",
      values_to = "value"
    )
    obs_pred_data$Data <- factor(obs_pred_data$Data, levels = c("OD", "OD_pred"),
                        labels = c("Observed", "Predicted"))
    
    
    # Plot observed and predicted lines
    ggplot(obs_pred_data, aes(x = time, y = value, colour = Data)) +
      geom_point(data = subset(obs_pred_data, Data == "Observed"), size = 1.5) +
      geom_line(data = subset(obs_pred_data, Data == "Predicted"), lwd = 1.5) +
      labs(
        title = paste0("Observed at ", input$dox, " vs Predicted OD", "in Experiment ", input$exp_id),
        y = expression("OD"[600] ~ "(log"[10]*")"),
        x = "Time (hours)",
        colour = "Data"
      ) +
      facet_wrap(~ concn, labeller = label_both, ncol = 6) +
      scale_color_manual(
        values = c("Observed" = "black", "Predicted" = "blue")
      ) +
      scale_x_continuous(
        limits = c(min(obs_pred_data$time, na.rm = TRUE), max(obs_pred_data$time, na.rm = TRUE)),
        breaks = seq(
          min(obs_pred_data$time, na.rm = TRUE),
          max(obs_pred_data$time, na.rm = TRUE), 
          by = 3
        )
      ) +
      scale_y_log10(
        limits = c(0.001, 2),
        minor_breaks = c(
          seq(1, 9, 1) * 0.001,
          seq(1, 9, 1) * 0.01,
          seq(1, 9, 1) * 0.1,
          seq(1, 9) * 1,
          2
        )
      ) +
      theme_minimal() +
      theme(
        title = element_text(size = 18, hjust = 0.5),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        strip.text = element_text(size = 14),
        legend.position = "bottom",
        legend.text = element_text(size = 14)
      )
  })
  
  ## Predicted proportions of dead cells
  output$deadcellPlot <- renderPlot({
    results <- fit_results()
    req(results)
    
    dead_pred_data_all <- do.call(rbind, lapply(names(results), function(name) {
      res_item <- results[[name]]
      if (!is.null(res_item$obs_pred_data)) {
        cbind(concn = name, res_item$obs_pred_data)
      } else {
        NULL  
      }
    }))
    
    req(nrow(dead_pred_data_all) > 0)
    
    dead_pred_data <- dead_pred_data_all %>%
                      dplyr::rename(Observed = D,
                                    Predicted = dead) %>%
                      tidyr::pivot_longer( cols = c(Observed, Predicted),
                                                names_to = "Data",
                                                values_to = "value")
    
    dead_pred_data$Data <- factor(dead_pred_data$Data, levels = c("Observed", "Predicted"))
    
    ggplot(dead_pred_data, aes(x = time, y = value, colour = Data)) +
      geom_point(data = subset(dead_pred_data, Data == "Observed"), size = 1.5) +
      geom_line(data = subset(dead_pred_data, Data == "Predicted"), lwd = 1.5) +
      labs(
        title = paste0("Observed at ", input$dox, " vs Predicted Viability", "in Experiment ", input$exp_id),
        y = "Proportion of Dead Cells",
        x = "Time (hours)",
        colour = "Data type"
      ) +
     facet_wrap(~ concn, labeller = label_both, ncol = 6) +
      scale_color_manual(
        values = c("Observed" = "black", "Predicted" = "blue")
      ) +
      scale_x_continuous(
          limits = c(min(dead_pred_data$time, na.rm = TRUE), max(dead_pred_data$time, na.rm = TRUE)),
          breaks = seq(
                      min(dead_pred_data$time, na.rm = TRUE),
                      max(dead_pred_data$time, na.rm = TRUE), 
                           by = 3)
                         )+
      scale_y_continuous(
        limits = c(0, 1),
        breaks = seq(0, 1, by = 0.25)
      ) +
      theme_minimal() +
      theme(
        title = element_text(size = 18, hjust = 0.5),
        strip.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.position = "bottom"
      )
  })
  
  ## MCMC plots
  # Update drop down dynamically when results change
  observe({
    results <- fit_results()
    req(results)
    
    updateSelectInput(
      session,
      "plot_choice",
      choices = names(results),
      selected = names(results)[1]
    )
  })
  
  # Render MCMC plot for selected result
  output$mcmcPlot <- renderPlot({
    results <- fit_results()
    req(results, input$plot_choice)
    
    selected_result <- results[[input$plot_choice]]
    req(!is.null(selected_result$mcmc))
    
    # Directly call the plotting function (it plots in base R)
    plot(selected_result$mcmc, Full = TRUE)
  })
  
  ## Sensitivity plot
  ## creating drop down options for all the DOX concentrations
  observe({
    results <- fit_results()
    req(results)
    
    updateSelectInput(
      session,
      "sensplot_choice",
      choices = names(results),
      selected = names(results)[1]
    )
  })
  
  # Render sensitivity plot for selected result
  output$sensPlot <- renderPlot({
    results <- fit_results()
    req(results, input$sensplot_choice)
    
    selected_result <- results[[input$sensplot_choice]]
    req(!is.null(selected_result$sens_func))
    
    # Directly call the plotting function (it plots in base R)
    plot(selected_result$sens_func, which = c("OD", "dead"), lwd = 2)
  })
  
  ## Collinearity Indices
  observe({
    results <- fit_results()
    req(results)
    
    updateSelectInput(
      session,
      "collinplot_choice",
      choices = names(results),
      selected = names(results)[1]
    )
  })
  
  # Render collinearity plot for selected result
  output$collinPlot <- renderPlot({
    results <- fit_results()
    req(results, input$collinplot_choice)
    
    selected_result <- results[[input$collinplot_choice]]
    req(!is.null(selected_result$collinearity))
    
    plot(selected_result$collinearity, log = "y")
  })
  
  ## Sense Range plots for the model predictions with MCMC applications
  # Render senseRange plot for selected result
  output$sensRangePlot <- renderPlot({
    results <- fit_results()
    req(results, filtered_od_data())
    od_df <- filtered_od_data() %>%
      dplyr::group_by(time, dox) %>%
      dplyr::summarise(
        mean_OD = mean(OD, na.rm = TRUE),
        sd_OD = sd(OD, na.rm = TRUE),
        n = dplyr::n(),
        .groups = "drop"
      ) %>%
      dplyr::rename(OD = mean_OD, sd = sd_OD) %>%
      dplyr::filter(OD > 0)
    
    sensR_all <- do.call(rbind, lapply(names(results), function(name) {
      res_item <- results[[name]]
      if (!is.null(res_item$sR)) {
        cbind(concn = name, res_item$sR)
      } else {
        NULL  
      }
    }))
    
  sensR_df <- sensR_all %>%
              tibble::rownames_to_column("label") %>%
              dplyr::mutate(
                state = dplyr::case_when(
                  grepl("^OD", label) ~ "OD",
                  grepl("^dead", label) ~ "dead",
                  TRUE ~ NA_character_),
                  dox = suppressWarnings(as.numeric(sub("dox_", "", concn)))) %>%
              dplyr::filter(state == "OD") %>%
              dplyr::rename(time = x) %>%
              dplyr::select(time, dox, Mean, q05, q50, q95) %>%
              dplyr::left_join(od_df, by = c("time", "dox")) %>%
              tidyr::drop_na() 
              
   ggplot(sensR_df, aes(x = time)) +
      geom_ribbon(aes(ymin = q05, ymax = q95, fill = "Prediction Interval (90%)"), alpha = 0.5) +
      geom_line(aes(y = q50, colour = "Predicted (median)")) +
      geom_point(aes(y = OD, colour = "Observed"), size = 1.5) +
      scale_color_manual(values = c("Predicted (median)" = "blue", "Observed" = "black"))+
      scale_fill_manual(values = c("Prediction Interval (90%)" = "gray80")) +
         guides(
           colour = guide_legend(order = 1),
       fill = guide_legend(order = 1)) +
      facet_wrap(~ dox, labeller = label_both, ncol = 6) +
      labs(
        title = paste0("VPC plot for ", 
                       input$strain, " at OD0 of ", input$OD0,
                       " in Experiment ", input$exp_id),
        y = expression("OD"[600] ~ "(log"[10]*")"),
        x = "Time (hours)",
        colour = "Data type",
        fill = ""
      ) +
      scale_x_continuous(
            limits = c(0, max(sensR_df$time, na.rm = TRUE)),
            breaks = seq(0, max(sensR_df$time, na.rm = TRUE), 
              by = 3)
      )+
      scale_y_log10(
        limits = c(0.001, 2),
        minor_breaks = c(seq(1, 9, 1)*0.001, seq(1, 9, 1) * 0.01,
                         seq(1, 9, 1) * 0.1, seq(1, 9) * 1, 2)
      ) +
      theme_minimal() +
      theme(
        title = element_text(size = 18),
        strip.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.position = "bottom"
      )
  })
    
  # Download handler for fitted parameters
  fitted_parms <- reactive({
    results <- fit_results()
    req(results)  
    
    do.call(rbind, lapply(names(results), function(name) {
      if (!is.null(results[[name]]) && !is.null(results[[name]]$ci_table)) {
        cbind(concn = name, results[[name]]$ci_table)
      } else {
        data.frame(
          conc = name,
          parameter = NA, 
          fitted_value = NA,
          lower_95_CI = NA, 
          upper_95_CI = NA
        )
      }
    }))
  })
  
  output$download_fitted_params <- downloadHandler(
    filename = function() {
      paste("fitted_parameters_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(fitted_parms(),
                file, row.names = FALSE)
    }
  )
}

shinyApp(ui, server)


###########################################################################################################################################