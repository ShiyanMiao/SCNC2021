# app.R — Shiny  # File name and framework used (Shiny app in R).
# model page: single-marker (y ~ X_j) and multiple-marker (y ~ all X_*) linear models.  # Describes two model types provided on the model tab.
# power page: NCP/λ method based on λ_j ≈ n · (β_j^2/σ^2) · 1/VIF_j with noncentral-t power;  # Explains power page uses lambda/NCP formula and noncentral-t power.
#             solver uses "halve-then-bisect" (start from big n, halve to bracket, then bisection).  # Notes sample-size solver strategy: halve then bisection.
# structure inspired by "Mastering Shiny" (Ch.1–Data/Ch.2–Inputs/Outputs/Ch.6–Layout/Ch.9–File upload).  # Credits app structure to Mastering Shiny chapters.

suppressPackageStartupMessages({  # Start a block to load packages quietly (no startup messages).
  library(shiny)            # Load Shiny core framework.
  library(DT)               # Load DT for interactive tables.
  library(shinydashboard)   # Load shinydashboard for infoBox components.
})  # End quiet package loading.

# -----------------------------  # Section divider: data generation.
# Data generation — split X and y  # We generate predictors X and response y separately for clarity.
# -----------------------------

# Generate clean X (genotypes; no missing values)  # Function to create an n-by-m matrix of markers.
sim_X_clean <- function(n = 200,      # individuals (rows)  # Argument: number of rows (sample size).
                        m = 50,       # markers (columns)   # Argument: number of columns (predictors).
                        maf = 0.3,    # minor allele frequency (HWE)  # Probability used in binomial draws.
                        seed = 1) {   # RNG seed for X     # Argument: seed to reproduce X.
  set.seed(seed) # Ensure reproducibility for X  # Fix random seed for consistent X.
  G <- matrix(rbinom(n * m, 2, maf), nrow = n, ncol = m)  # Fill matrix by binomial(2, maf) draws then reshape to n×m.
  colnames(G) <- paste0("X_", seq_len(m))  # Name columns X_1, X_2, ..., X_m.
  G  # Return the X matrix.
}  # End sim_X_clean.

# Generate y given X (linear model) — flexible effect size & error SD  # Function to simulate y from linear model.
sim_y_given_X <- function(X,  # Full predictor matrix used to form the linear predictor.
                          n_causal = 3,  # number of causal markers (first n_causal columns used)  # Count of nonzero betas (prefix).
                          beta = 0.2,    # effect size per causal marker (no upper bound)  # Coefficient magnitude for causal X.
                          error_sd = 1,  # error SD (σ), now user-controlled  # Standard deviation of noise.
                          seed = 1) {    # RNG seed for y  # Seed to reproduce y draws.
  stopifnot(is.matrix(X))  # Validate X is a matrix to avoid recycling or class issues.
  set.seed(seed) # Ensure reproducibility for y  # Fix random seed for noise.
  m <- ncol(X)  # Number of predictors (used to size coefficient vector).
  b <- rep(0, m)                      # start with no effects  # Initialize all betas to zero.
  if (n_causal > 0) {  # Only set effects when at least one causal variable requested.
    b[seq_len(min(n_causal, m))] <- beta  # set first n_causal to given effect size  # Assign beta to first causal positions.
  }
  eta <- as.numeric(X %*% b)          # linear predictor  # Compute Xβ as numeric vector.
  y <- eta + rnorm(nrow(X), 0, error_sd)  # add noise with chosen SD (σ)  # Add N(0, σ^2) noise to get y.
  y  # Return simulated response vector.
}  # End sim_y_given_X.

# Single-marker scanning（y ~ X_j）  # Routine to fit y on each marker separately.
single_marker_scan <- function(dat) {  # Accept a data.frame: first column y, others X_* predictors.
  # Fit y ~ X_j for each marker independently and collect beta_hat, t, p.  # Describe the approach per marker.
  x_cols <- names(dat)[startsWith(names(dat), "X_")]  # Pick predictor columns whose names begin with "X_".
  out <- lapply(x_cols, function(v) {  # Loop over each marker column name.
    fit <- lm(dat$y ~ dat[[v]])  # Fit simple linear regression y ~ X_j for this marker.
    co <- summary(fit)$coefficients  # Extract coefficient table (estimate, SE, t, p).
    c(beta_hat = unname(co[2, 1]),  # Store slope estimate for X_j (remove dimnames).
      t        = unname(co[2, 3]),  # Store t-statistic for the X_j coefficient.
      p        = unname(co[2, 4]))  # Store two-sided p-value for the X_j coefficient.
  })  # End lapply over markers.
  as.data.frame(do.call(rbind, out), row.names = x_cols)  # Stack rows and set rownames to marker IDs.
}  # End single_marker_scan.

# -----------------------------  # Section divider: UI.
# UI  # Define the user interface with a navbar layout.
# -----------------------------
ui <- navbarPage(  # Create a multi-tab navigation bar UI.
  title = "QTL Design Sandbox",  # App title displayed in the navbar.
  
  # ---- data ----  # Data tab: generate/upload and preview data.
  tabPanel("data", sidebarLayout(  # Add "data" tab with a sidebar and main area.
    sidebarPanel(  # Left panel with input controls for data management.
      h4("Active data source"),  # Small heading for source selection.
      radioButtons("active_source", "Use data from",  # Radio to switch between generated vs uploaded data.
                   choices = c("Generated" = "generated", "Uploaded" = "uploaded"),  # Two data sources to choose.
                   selected = "generated"),  # Default to generated data.
      tags$hr(),  # Visual separator line.
      
      h4("Generate clean X"),  # Heading for X-generation parameters.
      numericInput("n_ind",  "Number of individuals", 200, min = 20, step = 20),  # Control sample size n.
      numericInput("n_mark", "Number of markers",      50,  min = 5,  step = 5),  # Control number of markers m.
      sliderInput("maf", "Minor allele frequency", min = 0.05, max = 0.5, value = 0.3, step = 0.01),  # Control MAF.
      numericInput("seed_x", "Seed for X", 1, step = 1),  # Seed for reproducible X.
      actionButton("btn_gen_X", "Generate X"),  # Button to trigger X generation.
      tags$hr(),  # Separator.
      
      h4("Generate y given X"),  # Heading for y-generation parameters.
      numericInput("n_causal", "Number of causal markers", 3, min = 0, step = 1),  # Control number of causal markers.
      numericInput("beta", "Effect size per causal marker (β)", 0.2, step = 0.05),  # Control beta magnitude.
      numericInput("error_sd", "Error SD (σ)", 1, min = 1e-6, step = 0.05),  # Control noise SD.
      numericInput("seed_y", "Seed for y", 1, step = 1),  # Seed for reproducible y.
      actionButton("btn_gen_Y", "Generate y"),  # Button to simulate y based on current X.
      tags$hr(),  # Separator.
      
      h4("Or upload CSV (clean)"),  # Heading for CSV upload option.
      fileInput("upload", "Upload CSV", accept = ".csv"),  # File chooser for CSV with y and X_* columns.
      helpText("Columns: y, X_1 ... X_m (no missing values). Uploaded and generated data both remain available; the active source drives displays and fits.")  # Guidance on CSV format and source behavior.
    ),  # End sidebarPanel.
    mainPanel(  # Right panel to display tables and plots.
      h4("Preview X (active source)"),  # Heading for X preview table.
      DTOutput("x_head"),  # DataTable placeholder for first rows of X.
      tags$hr(),  # Separator.
      h4("Preview data (y + X) — active source"),  # Heading for combined data preview.
      DTOutput("data_head"),  # DataTable placeholder for y and X_*.
      tags$hr(),  # Separator.
      h4("y histogram (active source)"),  # Heading for histogram of y.
      plotOutput("y_hist", height = 240),  # Plot area to display histogram of y.
      tags$hr(),  # Separator.
      verbatimTextOutput("data_info")  # Text area summarizing dataset structure.
    )  # End mainPanel.
  )),  # End tabPanel "data".
  
  # ---- model ----  # Model tab: run SIM and MMM fits.
  tabPanel("model", sidebarLayout(  # Add "model" tab with sidebar and main.
    sidebarPanel(  # Left panel with model controls.
      radioButtons(  # Let user choose modeling mode.
        "model_kind", "Model type",  # Input id and label for model selection.
        c(
          "Single-marker models (y ~ X_j)"      = "single",  # Option for per-marker regressions.
          "Multiple-marker model (y ~ all X_*)" = "multi"    # Option for full OLS with all markers.
        ),
        selected = "single"  # Default mode is single-marker scan.
      ),
      actionButton("btn_fit", "Fit model")  # Button to trigger model fitting.
    ),  # End model sidebar.
    mainPanel(  # Right panel for model outputs.
      conditionalPanel(  # Show this block only when single-marker mode is selected.
        condition = "input.model_kind == 'single'",  # Client-side condition on selection.
        h4("Single-marker scan results"),  # Heading for SIM results table.
        DTOutput("single_res")  # DataTable placeholder for SIM statistics.
      ),
      conditionalPanel(  # Show this block only when multiple-marker mode is selected.
        condition = "input.model_kind == 'multi'",  # Client-side condition on selection.
        h4("Multiple-marker model summary: y ~ all X_*"),  # Heading for MMM summary.
        verbatimTextOutput("multi_summary")  # Text output for summary(lm(...)).
      ),
      tags$hr(),  # Separator between results and hints.
      uiOutput("model_hint"),  # Dynamic UI with scale checks and notes.
      uiOutput("model_status")  # shows if results are stale after data refresh  # Warns when data changed after last fit.
    )  # End model main panel.
  )),  # End tabPanel "model".
  
  # ---- power (λ method only) ----  # Power tab using lambda/NCP approach only.
  tabPanel("power", sidebarLayout(  # Add "power" tab with sidebar and main.
    sidebarPanel(  # Left panel for power inputs.
      helpText("NCP/λ method (single coefficient): λ_j ≈ n · (β_j^2/σ^2) · 1/VIF_j; power via noncentral-t."),  # Brief method note.
      numericInput("n_fixed", "n to check (given n)", 1000, min = 10, step = 10),  # Input n to evaluate power at.
      numericInput("alpha_lambda", "Significance level (α)", 0.05, min = 1e-8, step = 0.01),  # Per-test alpha for power.
      numericInput("beta_j", "Assumed β_j", 0.2, step = 0.01),  # Assumed effect size for a single coefficient.
      numericInput("sigma_j", "Assumed σ (error SD)", 1, min = 1e-8, step = 0.01),  # Assumed residual std deviation.
      checkboxInput("manual_vif", "Enter VIF manually", FALSE),  # Toggle to input VIF by hand.
      conditionalPanel(  # Show numeric VIF input only when manual mode is checked.
        condition = "input.manual_vif == true",
        numericInput("vif_manual", "VIF_j (manual)", 1, min = 1, step = 0.1)  # Manual VIF value (≥1).
      ),
      conditionalPanel(  # Otherwise, allow choosing a marker to compute VIF from data.
        condition = "input.manual_vif == false",
        uiOutput("marker_choice_ui"),  # Dynamic selectInput listing X_* columns.
        helpText("If data available, VIF_j is computed from X_j ~ X_-j; otherwise switch to manual VIF.")  # Note on VIF estimation.
      ),
      numericInput("target_power_lambda", "Target power (default 0.90)", 0.90, min = 0.5, max = 0.999, step = 0.01),  # Desired power level.
      actionButton("btn_lambda", "Compute power / solve n (λ method)")  # Trigger computing power and minimal n.
    ),  # End power sidebar.
    mainPanel(  # Right panel for power results.
      h4("Required sample size (λ method)"),  # Heading for main sample-size result.
      infoBoxOutput("ibox_n", width = 12),  # headline using infoBox  # Info box to display n* and quick status.
      tags$hr(),  # Separator.
      h4("NCP/λ results"),  # Heading for detailed lambda results.
      uiOutput("lambda_text")  # Dynamic block with powers, interval, and narrative.
    )  # End power main panel.
  ))  # End tabPanel "power".
)  # End UI navbarPage.

# -----------------------------  # Section divider: server.
# SERVER  # Define server logic handling reactivity, models, and power.
# -----------------------------
server <- function(input, output, session) {  # Main server function with Shiny inputs/outputs/session.
  
  # Store generated and uploaded data SEPARATELY and choose an active source.  # Keep both sources and switch by radio.
  X_gen   <- reactiveVal(NULL)  # generated X (matrix)  # Holds last generated X.
  y_gen   <- reactiveVal(NULL)  # generated y (numeric)  # Holds last generated y.
  df_gen  <- reactiveVal(NULL)  # generated full data.frame(y, X_*)  # Holds combined data for generated source.
  df_upl  <- reactiveVal(NULL)  # uploaded full data.frame(y, X_*)  # Holds uploaded dataset.
  
  stale_model <- reactiveVal(TRUE)  # Flag whether model outputs are stale (data changed after last fit).
  stale_power <- reactiveVal(TRUE)  # Flag whether power outputs are stale (inputs/data changed after last compute).
  
  # Helper: rule-of-thumb total n for regression (very rough)  # n ≳ 10p heuristic as quick scale check.
  rule_of_thumb_n <- function(p) ceiling(10 * max(1, p))  # Return ceil(10*p) with p ≥ 1.
  
  # ---- safe default infoBox at startup ----  # Initialize power info box to a neutral prompt.
  output$ibox_n <- renderInfoBox({  # Render the infoBox for power result.
    shinydashboard::infoBox(  # Create a colored info box widget.
      title = "Sample size",  # Box title label.
      value = "—",  # Placeholder value (dash).
      subtitle = "Click 'Compute power / solve n' to compute.",  # Prompt text for user action.
      icon = icon("calculator"), color = "light-blue", fill = TRUE  # Style with icon/color/fill.
    )
  })  # End default infoBox render.
  
  # ---- Source switching helpers (do not overwrite model outputs) ----  # Functions to switch active source and mark stale.
  use_generated <- function() {  # Switch to generated source.
    updateRadioButtons(session, "active_source", selected = "generated")  # Update radio input to “generated”.
    stale_model(TRUE); stale_power(TRUE)  # Mark downstream outputs as stale.
  }  # End use_generated.
  use_uploaded <- function() {  # Switch to uploaded source.
    updateRadioButtons(session, "active_source", selected = "uploaded")  # Update radio input to “uploaded”.
    stale_model(TRUE); stale_power(TRUE)  # Mark downstream outputs as stale.
  }  # End use_uploaded.
  observeEvent(input$active_source, { stale_model(TRUE); stale_power(TRUE) }, ignoreInit = TRUE)  # Whenever source changes, flag outputs stale.
  
  # ---- Upload CSV (kept separate; does NOT erase generated data) ----  # Read and store uploaded data.
  observeEvent(input$upload, {  # React when a CSV is uploaded.
    file <- input$upload; req(file)  # Ensure a file is present.
    dat <- read.csv(file$datapath, check.names = FALSE)  # Read CSV without altering column names.
    validate(need("y" %in% names(dat), "Uploaded CSV must contain column 'y'."))  # Require response column y.
    validate(need(any(startsWith(names(dat), "X_")), "Uploaded CSV must contain X_* columns."))  # Require predictor columns.
    df_upl(dat)            # store uploaded  # Save to uploaded data reactive.
    use_uploaded()         # switch active source to uploaded  # Make uploaded the active source.
  }, ignoreInit = TRUE)  # Do not run at app start.
  
  # ---- Generate X ----  # Build X matrix using current controls.
  observeEvent(input$btn_gen_X, {  # React when Generate X is clicked.
    G <- sim_X_clean(n = input$n_ind, m = input$n_mark, maf = input$maf, seed = input$seed_x)  # Create X with parameters.
    X_gen(G)     # store generated X  # Save X to reactive.
    y_gen(NULL)  # y not yet generated  # Clear y because X changed.
    df_gen(NULL) # full data not ready  # Clear combined data until y is generated.
    use_generated()  # Switch active source to generated and mark stale.
  }, ignoreInit = TRUE)  # Ignore on startup.
  
  # ---- Generate y given X ----  # Build y using current X and y parameters.
  observeEvent(input$btn_gen_Y, {  # React when Generate y is clicked.
    req(X_gen())  # Require X exists before generating y.
    y <- sim_y_given_X(X = X_gen(),  # Pass the current generated X.
                       n_causal = input$n_causal,  # Use current number of causal markers.
                       beta     = input$beta,      # Use current effect size.
                       error_sd = input$error_sd,  # Use current noise SD.
                       seed     = input$seed_y)    # Use seed for y.
    y_gen(y)  # Save generated y.
    df_gen(data.frame(y = y, X_gen(), check.names = FALSE))  # Build combined data frame y + X_*.
    use_generated()  # Keep generated as active and mark stale.
  }, ignoreInit = TRUE)  # Ignore on startup.
  
  # ---- Current active data helpers ----  # Helpers to provide active full data and X-only.
  current_data <- reactive({  # Return active full dataset (y + X) or NULL.
    if (input$active_source == "uploaded" && !is.null(df_upl())) return(df_upl())  # Prefer uploaded when active and present.
    if (input$active_source == "generated" && !is.null(df_gen())) return(df_gen())  # Otherwise use generated when present.
    NULL  # No full data yet.
  })  # End current_data.
  
  # For power/VIF computing: allow X-only (no y) when using generated X  # Return X matrix even when y missing.
  current_X_only <- reactive({  # Provide X-only for VIF or power tab.
    if (input$active_source == "generated" && !is.null(X_gen())) return(X_gen())  # Use generated X when active.
    if (input$active_source == "uploaded" && !is.null(df_upl())) {  # If uploaded active, extract X_* columns.
      x_cols <- names(df_upl())[startsWith(names(df_upl()), "X_")]  # Identify X_* columns.
      if (length(x_cols) > 0) return(as.matrix(df_upl()[, x_cols, drop = FALSE]))  # Return as matrix if any present.
    }
    NULL  # No X available.
  })  # End current_X_only.
  
  # ---------------- Data tab outputs ----------------  # Render tables/plots and info for data tab.
  
  output$x_head <- renderDT({  # Show first rows of X depending on active source.
    if (input$active_source == "generated") {  # When generated source selected.
      if (is.null(X_gen()))  # If no generated X yet.
        return(datatable(data.frame(Info = "No generated X. Click 'Generate X'."), options = list(dom = 't')))  # Show message table.
      datatable(head(as.data.frame(X_gen()), 10),  # Otherwise show first 10 rows of generated X.
                options = list(scrollX = TRUE, pageLength = 10))  # Enable scrolling and page length.
    } else {  # Uploaded source branch.
      if (is.null(df_upl()))  # No uploaded dataset yet.
        return(datatable(data.frame(Info = "No uploaded data. Upload a CSV."), options = list(dom = 't')))  # Show message table.
      x_cols <- names(df_upl())[startsWith(names(df_upl()), "X_")]  # Identify X_* columns in uploaded data.
      if (!length(x_cols))  # If no predictors present.
        return(datatable(data.frame(Info = "Uploaded data has no X_* columns."), options = list(dom = 't')))  # Show warning table.
      datatable(head(df_upl()[, x_cols, drop = FALSE], 10),  # Show first 10 rows of X_* from uploaded data.
                options = list(scrollX = TRUE, pageLength = 10))  # Enable scrolling and page length.
    }
  })  # End output$x_head.
  
  output$data_head <- renderDT({  # Show first rows of full dataset (y + X_*).
    dat <- current_data()  # Get active full dataset.
    if (is.null(dat))  # If not available yet.
      return(datatable(data.frame(Info = "No full data (y + X). Generate y or upload CSV."), options = list(dom = 't')))  # Show message.
    datatable(head(dat, 10), options = list(scrollX = TRUE, pageLength = 10))  # Otherwise show first 10 rows.
  })  # End output$data_head.
  
  output$y_hist <- renderPlot({  # Plot histogram of y for the active source.
    dat <- current_data()  # Get active full dataset.
    if (is.null(dat) || !"y" %in% names(dat)) {  # If y not present yet.
      plot.new(); title("No y yet (generate y or upload CSV including y).")  # Show placeholder plot with message.
      return(invisible())  # Exit silently after placeholder.
    }
    hist(dat$y, breaks = 30, main = "Histogram of y (active source)", xlab = "y")  # Draw histogram with 30 bins.
  })  # End output$y_hist.
  
  output$data_info <- renderPrint({  # Print dataset info in text form.
    src <- input$active_source  # Read which source is active.
    if (src == "generated") {  # Generated branch.
      if (is.null(X_gen()) && is.null(df_gen())) {  # If neither X nor y generated.
        cat("Active source: Generated — No X or y yet.\n"); return(invisible())  # Print note and exit.
      }
      if (!is.null(df_gen())) {  # If full generated dataset exists.
        dat <- df_gen()  # Get generated full data.
        p <- sum(startsWith(names(dat), "X_"))  # Count markers.
        cat("Active source: Generated\nRows:", nrow(dat), "  Markers:", p, "\nColumns:", paste(names(dat), collapse = ", "), "\n")  # Print stats.
      } else {  # Only X exists (y not generated).
        cat("Active source: Generated (X only). Rows:", nrow(X_gen()), "  Markers:", ncol(X_gen()), "\n")  # Print X-only stats.
      }
    } else {  # Uploaded branch.
      if (is.null(df_upl())) { cat("Active source: Uploaded — none.\n"); return(invisible()) }  # No uploaded data present.
      dat <- df_upl()  # Get uploaded data.
      p <- sum(startsWith(names(dat), "X_"))  # Count markers.
      cat("Active source: Uploaded\nRows:", nrow(dat), "  Markers:", p, "\nColumns:", paste(names(dat), collapse = ", "), "\n")  # Print stats.
    }
  })  # End output$data_info.
  
  # ---------------- Model tab ----------------  # SIM and MMM logic and renders.
  
  # Single-marker scan (event-based; stale after data changes)  # Compute SIM results when user clicks Fit.
  single_res <- eventReactive(input$btn_fit, {  # Recompute when Fit button is pressed.
    dat <- current_data()  # Get active full dataset.
    validate(need(!is.null(dat), "No full data. Generate y or switch to an uploaded dataset with y."))  # Require dataset presence.
    validate(need("y" %in% names(dat), "Missing column 'y'."))  # Require y column.
    validate(need(any(startsWith(names(dat), "X_")), "No X_* columns found."))  # Require predictors.
    stale_model(FALSE)  # results are now fresh  # Mark model outputs as up to date.
    single_marker_scan(dat)  # Run SIM and return results table.
  }, ignoreInit = TRUE)  # Do not run at startup.
  
  output$single_res <- renderDT({  # Render the SIM results as an interactive table.
    res <- single_res()  # Get SIM results.
    req(res)  # Ensure results exist.
    datatable(  # Create DataTable widget.
      transform(  # Post-process columns for nicer display.
        cbind(marker = rownames(res), res),  # Add marker name as a first column.
        p        = signif(p, 4),  # Format p-values to 4 significant digits.
        beta_hat = signif(beta_hat, 4),  # Format estimates to 4 significant digits.
        t        = round(t, 3)  # Round t-statistics to 3 decimals.
      ),
      options = list(pageLength = 10)  # Show 10 rows per page.
    )
  })  # End output$single_res.
  
  output$multi_summary <- renderPrint({  # Render summary of multiple-marker OLS fit.
    req(input$btn_fit > 0)  # Only after Fit button has been clicked at least once.
    dat <- current_data()  # Get active full dataset.
    validate(need(!is.null(dat), "No full data. Generate y or switch to uploaded data."))  # Must have data.
    x_cols <- names(dat)[startsWith(names(dat), "X_")]  # Identify predictor columns.
    validate(need(length(x_cols) > 0, "No X_* columns to fit."))  # Must have predictors.
    validate(need(nrow(dat) > length(x_cols) + 1,  # Ensure degrees of freedom > 0.
                  sprintf("Multiple-marker OLS requires n > p + 1. Current n=%d, p=%d. Increase n or reduce p.",
                          nrow(dat), length(x_cols))))  # Give a helpful message if underdetermined.
    fit <- lm(reformulate(x_cols, "y"), data = dat)   # MMM = OLS with all markers  # Fit y ~ X_1 + ... + X_p.
    stale_model(FALSE)  # results are now fresh  # Mark model outputs as up to date.
    summary(fit)  # Print standard lm summary (coefficients, R^2, etc.).
  })  # End output$multi_summary.
  
  output$model_hint <- renderUI({  # Render a small guidance block about scale and methods.
    if (is.null(current_data()))
      return(tagList(h4("Scale check / intuition"), p("No data yet.")))  # Show placeholder if no data.
    dat <- current_data()  # Get dataset.
    p <- sum(startsWith(names(dat), "X_"))  # Count predictors.
    n <- nrow(dat)  # Read sample size.
    tagList(  # Return a group of UI elements.
      h4("Scale check / intuition"),  # Section title.
      p(sprintf("Current: n = %d, p = %d.", n, p)),  # Show current n and p.
      p(sprintf("Rule-of-thumb: n \u2265 10 \u00D7 p \u2248 %d (for linear regression).",  # Heuristic size guidance.
                rule_of_thumb_n(p))),
      p("SIM: fit y ~ X_j for each marker (separate OLS)."),  # Quick description for SIM.
      p("MMM: fit y ~ X_1 + ... + X_p by OLS; requires n > p + 1 (else consider penalized models, not covered here).")  # MMM caveat.
    )
  })  # End output$model_hint.
  
  output$model_status <- renderUI({  # Render stale status badge for model results.
    if (isTRUE(stale_model())) {  # If stale, ask user to refresh.
      tags$em("Model outputs are outdated due to data changes. Click 'Fit model' to refresh.")
    } else {
      NULL  # Otherwise show nothing.
    }
  })  # End output$model_status.
  
  # ---------------- Power tab (λ method only) ----------------  # Lambda/NCP power utilities and UI bindings.
  
  # Compute VIF_j from current X (if possible)  # Estimate VIF by regressing X_j on all other X.
  compute_vif_j <- function(X, j) {  # X is the design matrix; j is the column index for the focal predictor.
    p <- ncol(X)  # Number of predictors.
    if (p <= 1) return(1)  # only one predictor -> VIF = 1  # No multicollinearity with single predictor.
    others <- setdiff(seq_len(p), j)  # Indices of all predictors except j.
    # Guard: need n > length(others) + 1 to fit OLS for R^2  # Ensure enough rows for regression.
    if (nrow(X) <= length(others) + 1) return(NA_real_)  # Return NA when insufficient data.
    df_tmp <- data.frame(xj = X[, j], X[, others, drop = FALSE])  # Build data frame with response xj and other predictors.
    fit <- try(lm(xj ~ ., data = df_tmp), silent = TRUE)  # Regress xj on remaining columns.
    if (inherits(fit, "try-error")) return(NA_real_)  # If lm failed, return NA VIF.
    r2 <- max(0, min(1, summary(fit)$r.squared))  # Clamp R^2 to [0,1] to avoid numerical issues.
    vif <- 1 / max(1e-12, (1 - r2))  # Compute VIF = 1/(1 - R^2) with small lower bound on denominator.
    vif  # Return VIF value.
  }  # End compute_vif_j.
  
  # Given n, beta, sigma, VIF, p, alpha -> power for single β_j by noncentral t (two-sided)  # Core power calculator.
  lambda_power_for_n <- function(n, beta, sigma, vif, p, alpha) {  # Compute power using NCP formula and t critical value.
    n <- as.integer(n)  # Force integer sample size.
    p <- as.integer(p)  # Force integer predictor count.
    if (n <= p + 1) return(0)  # not identifiable  # Return 0 power if model not estimable.
    nu <- n - p - 1  # residual df for OLS  # Degrees of freedom for residuals.
    # quoted formula (knowledge point): λ_j ≈ n · (β_j^2 / σ^2) · 1/VIF_j  # Define lambda as per quoted approximation.
    lambda <- n * (beta^2 / (sigma^2)) * (1 / vif)  # Compute lambda (noncentrality parameter squared).
    delta  <- sqrt(lambda)   # noncentrality for t  # Noncentrality parameter for t distribution.
    tcrit  <- qt(1 - alpha/2, df = nu)  # Two-sided critical t value at given alpha.
    # two-sided power = P(|T_{nu,delta}| > tcrit)  # Power is tail probability beyond ±tcrit for noncentral t.
    pow <- pt(-tcrit, df = nu, ncp = delta) + (1 - pt(tcrit, df = nu, ncp = delta))  # Sum of both tail probabilities.
    as.numeric(max(0, min(1, pow)))  # Clamp to [0,1] and return numeric scalar.
  }  # End lambda_power_for_n.
  
  # Halve-then-bisect solver: start from a big n (e.g., 1000), halve to bracket, then bisection to target power  # n* solver.
  solve_n_lambda_halving_then_bisect <- function(target_power, beta, sigma, vif, p, alpha,  # Find minimal n for given target power.
                                                 n_init = 1000, tol = 1, max_iter = 60) {  # Controls: starting n, tolerance, max iterations.
    # Ensure an initial upper bound hi (identifiable)  # Make sure hi is at least p+2 for identifiability.
    hi <- max(as.integer(n_init), p + 2)  # Set initial upper bound n.
    
    # If starting power is too small, grow hi until it reaches/overcomes target  # Expand hi until power(hi) ≥ target.
    pow_hi <- lambda_power_for_n(hi, beta, sigma, vif, p, alpha)  # Compute power at hi.
    grow_guard <- 0  # Safety counter for loop.
    while (pow_hi < target_power && hi < 1e6 && grow_guard < 30) {  # Increase hi if still below target power.
      hi <- as.integer(ceiling(hi * 1.5))  # Enlarge hi by 1.5× to accelerate search.
      pow_hi <- lambda_power_for_n(hi, beta, sigma, vif, p, alpha)  # Recompute power at new hi.
      grow_guard <- grow_guard + 1  # Increment guard.
    }
    if (pow_hi < target_power) {  # If even after growth target is not reachable.
      return(list(n_star = NA_integer_, lo = NA_integer_, hi = NA_integer_,  # Return NA to signal failure.
                  pow_lo = NA_real_, pow_hi = NA_real_))  # Include unset powers.
    }
    
    # Halving: shrink hi by halves to find lo with power(lo) < target <= power(hi)  # Bracket target by halving hi.
    lo <- max(p + 2, floor(hi / 2))  # Set tentative lo as half of hi (respect identifiability).
    pow_lo <- lambda_power_for_n(lo, beta, sigma, vif, p, alpha)  # Compute power at lo.
    
    shrink_guard <- 0  # Safety counter for halving loop.
    while (pow_lo >= target_power && lo > (p + 2) && shrink_guard < 30) {  # If lo still meets target, keep halving.
      hi <- lo  # Move hi down to lo (tighten upper bound).
      pow_hi <- pow_lo  # Update upper power.
      lo <- max(p + 2, floor(lo / 2))  # Halve lo again.
      pow_lo <- lambda_power_for_n(lo, beta, sigma, vif, p, alpha)  # Recompute power at new lo.
      shrink_guard <- shrink_guard + 1  # Increment guard.
    }
    # Now expect pow_lo < target_power <= pow_hi  # Bracket should contain target.
    if (!(pow_lo < target_power && pow_hi >= target_power)) {  # If bracketing failed (numerical oddities).
      # Fallback: return current hi as best  # Return hi as conservative n*.
      return(list(n_star = hi, lo = lo, hi = hi, pow_lo = pow_lo, pow_hi = pow_hi))
    }
    
    # Bisection within [lo, hi] until interval length <= tol  # Standard binary search for minimal n.
    it <- 0  # Iteration counter.
    while ((hi - lo) > tol && it < max_iter) {  # Continue until interval small enough or max iterations hit.
      it  <- it + 1  # Increment iteration.
      mid <- as.integer(floor((lo + hi) / 2))  # Middle n between lo and hi.
      pow_mid <- lambda_power_for_n(mid, beta, sigma, vif, p, alpha)  # Power at mid.
      if (pow_mid >= target_power) {  # If mid achieves target,
        hi <- mid  # move upper bound down to mid.
        pow_hi <- pow_mid  # Update power at hi.
      } else {  # Otherwise,
        lo <- mid  # move lower bound up to mid.
        pow_lo <- pow_mid  # Update power at lo.
      }
    }
    
    list(n_star = hi, lo = lo, hi = hi, pow_lo = pow_lo, pow_hi = pow_hi)  # Return minimal n and final bracket/powers.
  }  # End solver.
  
  # UI for choosing marker j when auto-computing VIF  # Build a selectInput listing markers.
  output$marker_choice_ui <- renderUI({  # Dynamic UI based on current X availability.
    Xonly <- current_X_only()  # Fetch X matrix if available.
    if (is.null(Xonly)) return(tags$em("No X available; use manual VIF."))  # If not, ask to use manual VIF.
    xnames <- colnames(Xonly)  # Get marker names.
    selectInput("marker_j", "Choose marker j", choices = xnames, selected = xnames[1])  # Build dropdown with markers.
  })  # End marker_choice_ui.
  
  # Headline infoBox updates after clicking λ button  # Compute headline result and render info box.
  observeEvent(input$btn_lambda, {  # React to Compute power / solve n click.
    isolate({  # Avoid unintended reactive dependencies while reading inputs.
      Xonly <- current_X_only()  # Try to fetch X to infer p and compute VIF.
      if (is.null(Xonly) && !isTRUE(input$manual_vif)) {  # If no X and no manual VIF,
        output$ibox_n <- renderInfoBox({  # Render a warning info box.
          shinydashboard::infoBox(
            title = "Sample size (λ)", value = "—",  # Empty value.
            subtitle = "No X available to compute VIF. Use manual VIF.",  # Guidance text.
            icon = icon("triangle-exclamation"), color = "red", fill = TRUE  # Warning styling.
          )
        })
        return()  # Exit early.
      }
      if (!is.null(Xonly)) {  # If X available,
        p <- ncol(Xonly)  # Use its number of columns as p.
      } else {  # Otherwise still missing p.
        output$ibox_n <- renderInfoBox({  # Render warning that p cannot be determined.
          shinydashboard::infoBox(
            title = "Sample size (λ)", value = "—",
            subtitle = "Please provide X to determine p (markers).",
            icon = icon("triangle-exclamation"), color = "red", fill = TRUE
          )
        })
        return()  # Exit early.
      }
      # VIF  # Decide VIF source (manual vs computed).
      if (isTRUE(input$manual_vif)) {  # Manual VIF path.
        vif_j <- max(1, as.numeric(input$vif_manual))  # Take nondecreasing VIF ≥ 1.
        marker_lab <- "(manual VIF)"  # Label for display.
      } else {  # Auto-computed VIF path.
        j <- match(input$marker_j, colnames(Xonly))  # Map selected marker name to column index.
        vif_j <- compute_vif_j(Xonly, j)  # Compute VIF via auxiliary regression.
        marker_lab <- if (is.finite(vif_j)) {  # Build label with VIF when finite.
          paste0(input$marker_j, sprintf(" (VIF = %.3f)", vif_j))
        } else "(VIF unavailable)"  # Otherwise show unavailable tag.
      }
      if (!is.finite(vif_j)) {  # If VIF not computable,
        output$ibox_n <- renderInfoBox({  # Show error info box prompting manual VIF.
          shinydashboard::infoBox(
            title = "Sample size (λ)", value = "—",
            subtitle = "Could not compute VIF; try manual VIF.",
            icon = icon("triangle-exclamation"), color = "red", fill = TRUE
          )
        })
        return()  # Exit early.
      }
      
      n_given <- as.integer(input$n_fixed)  # Read n to check power at.
      alpha   <- input$alpha_lambda  # Read significance level alpha.
      beta    <- input$beta_j  # Read assumed effect size β_j.
      sigma   <- input$sigma_j  # Read assumed residual SD σ.
      target  <- input$target_power_lambda  # Read target power.
      
      # Compute power at given n; solve minimal n via "halve-then-bisect"  # Evaluate both power(n_given) and n*.
      pow_n <- lambda_power_for_n(n_given, beta, sigma, vif_j, p, alpha)  # Power at user-specified n.
      res   <- solve_n_lambda_halving_then_bisect(  # Find minimal n to meet target power.
        target_power = target,
        beta = beta, sigma = sigma, vif = vif_j, p = p, alpha = alpha,
        n_init = n_given, tol = 1, max_iter = 60
      )
      n_min  <- res$n_star  # Extract minimal n.
      lo_fin <- res$lo  # Extract final lower bracket.
      hi_fin <- res$hi  # Extract final upper bracket (also n*).
      pow_lo <- res$pow_lo  # Power at lower bracket.
      pow_hi <- res$pow_hi  # Power at upper bracket.
      
      if (is.na(n_min)) {  # If solver could not achieve target,
        output$ibox_n <- renderInfoBox({  # Show unattainable status.
          shinydashboard::infoBox(
            title = "Sample size (λ)",
            value = "—",
            subtitle = "Target power may be unattainable under current settings.",
            icon = icon("triangle-exclamation"), color = "red", fill = TRUE
          )
        })
      } else {  # Otherwise display the minimal n* and status at n_given.
        meets <- if (n_given >= n_min) "meets target" else "below target"  # Compare n_given to n*.
        output$ibox_n <- renderInfoBox({  # Render info box with n*.
          shinydashboard::infoBox(
            title = "Sample size (λ)",
            value = paste0("n* ≈ ", n_min),
            subtitle = sprintf("Given n=%d → power≈%.3f (%s)", n_given, pow_n, meets),
            icon = icon("calculator"), color = "blue", fill = TRUE
          )
        })
      }
      stale_power(FALSE)  # Mark power results as fresh.
    })
  })  # End observeEvent for λ button.
  
  # Detailed λ results text (also shows final bracket and endpoint powers)  # Render narrative power block.
  output$lambda_text <- renderUI({  # Build detailed UI after button press.
    input$btn_lambda  # Depend on the compute button.
    isolate({  # Avoid extra reactive triggers while assembling text.
      Xonly <- current_X_only()  # Try to get X to infer p and compute VIF if needed.
      if (is.null(Xonly) && !isTRUE(input$manual_vif)) {  # If neither X nor manual VIF available,
        return(tags$em("No X available to compute VIF. Check 'Enter VIF manually'."))  # Prompt for manual VIF.
      }
      if (!is.null(Xonly)) {  # If X is available,
        p <- ncol(Xonly)  # Set p from its columns.
      } else {  # Otherwise we cannot compute power without p.
        return(tags$em("Please provide X to determine p (number of markers)."))  # Ask for X.
      }
      
      if (isTRUE(input$manual_vif)) {  # Manual VIF path for text block.
        vif_j <- max(1, as.numeric(input$vif_manual))  # Read VIF ensuring ≥1.
        marker_lab <- "(manual VIF)"  # Label as manual.
      } else {  # Auto VIF path.
        j <- match(input$marker_j, colnames(Xonly))  # Map selection to index.
        if (is.na(j)) return(tags$em("Invalid marker selection."))  # Guard invalid selection.
        vif_j <- compute_vif_j(Xonly, j)  # Compute VIF by regression.
        if (!is.finite(vif_j))
          return(tags$em("Could not compute VIF from X; dimension may be too high. Use manual VIF."))  # Guard computation failure.
        marker_lab <- paste0(input$marker_j, sprintf(" (VIF = %.3f)", vif_j))  # Include VIF in label.
      }
      
      n_given <- as.integer(input$n_fixed)  # Read evaluation n.
      alpha   <- input$alpha_lambda  # Read alpha.
      beta    <- input$beta_j  # Read beta.
      sigma   <- input$sigma_j  # Read sigma.
      target  <- input$target_power_lambda  # Read target power.
      
      pow_n <- lambda_power_for_n(n_given, beta, sigma, vif_j, p, alpha)  # Power at n_given.
      res   <- solve_n_lambda_halving_then_bisect(  # Find n* by halving plus bisection.
        target_power = target,
        beta = beta, sigma = sigma, vif = vif_j, p = p, alpha = alpha,
        n_init = n_given, tol = 1, max_iter = 60
      )
      n_min  <- res$n_star  # Minimal n meeting target power.
      lo_fin <- res$lo  # Final lower bracket.
      hi_fin <- res$hi  # Final upper bracket (=n*).
      pow_lo <- res$pow_lo  # Power at lo bracket.
      pow_hi <- res$pow_hi  # Power at hi bracket.
      
      tagList(  # Build a rich text block with details.
        p(HTML(sprintf("<b>Marker:</b> %s; <b>p</b> = %d; <b>α</b> = %.3g; <b>β<sub>j</sub></b> = %.3f; <b>σ</b> = %.3f.",
                       marker_lab, p, alpha, beta, sigma))),  # Summarize key inputs.
        p(HTML("Using quoted formula: λ<sub>j</sub> ≈ n · (β<sub>j</sub><sup>2</sup>/σ<sup>2</sup>) · 1/VIF<sub>j</sub>.")),  # Cite lambda formula.
        h4(sprintf("Given n = %d → power ≈ %.3f", n_given, pow_n)),  # Show power at given n.
        if (is.na(n_min)) {  # If not attainable,
          tags$em("Target power may be unattainable under current settings.")  # Explain limitation.
        } else {  # Otherwise report n* and bracket.
          tagList(
            h4(sprintf("Minimal n to reach target power (%.2f): n* ≈ %d", target, n_min)),  # Headline n*.
            p(sprintf("Final bracket: [%d, %d] with power(lo)=%.3f, power(hi)=%.3f",
                      lo_fin, hi_fin, pow_lo, pow_hi))  # Show the final interval and endpoint powers.
          )
        },
        p(if (!is.na(n_min)) {  # Conclude with actionable comparison to current n.
          if (n_given >= n_min) "Current n meets the target power." else "Current n is below the target; increase n."
        } else { "" })  # If unattainable, leave empty.
      )
    })
  })  # End output$lambda_text.
}  # End server function.

shinyApp(ui, server)  # Launch the Shiny app by pairing UI and server.
