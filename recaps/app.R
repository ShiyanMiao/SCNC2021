# app.R — Shiny application for QTL design sandbox (English-only UI)
# Every line includes an English comment for clarity and maintainability.

suppressPackageStartupMessages({                      # Suppress package startup messages for a clean console
  library(shiny)                                      # Core Shiny package for web apps
  library(DT)                                         # DataTables output for interactive tables
  library(shinydashboard)                             # Dashboard components such as infoBox
  library(ggplot2)                                    # Grammar of graphics for plotting
  library(htmltools)                                  # HTML construction utilities
  library(base64enc)                                  # Base64 encoding for embedding images in HTML
})

# -----------------------------
# Data generation — split X and y
# -----------------------------

sim_X_clean <- function(n = 200,                      # Function to simulate a clean genotype matrix X
                        m = 50,                       # Number of markers (columns)
                        maf = 0.3,                    # Minor allele frequency for binomial draws
                        seed = 1) {                   # Random seed for reproducibility
  set.seed(seed)                                      # Set seed so results are reproducible
  G <- matrix(                                        # Create a matrix of genotypes
    rbinom(n * m, 2, maf),                            # Draw genotypes assuming 0..2 copies under HWE
    nrow = n, ncol = m                                # Shape matrix to n rows by m columns
  )
  colnames(G) <- paste0("X_", seq_len(m))             # Name columns as X_1, X_2, ... X_m
  G                                                   # Return generated matrix
}

sim_y_given_X <- function(X,                           # Function to simulate phenotype y given X
                          n_causal = 3,                # Number of causal markers (first n_causal in order)
                          beta = 0.2,                  # Effect size per causal marker
                          error_sd = 1,                # Residual standard deviation
                          seed = 1,                    # Random seed
                          pve_target = NULL,           # Optional target proportion of variance explained
                          pve_mode = c("scale_beta",   # Mode for target PVE: scale betas
                                       "set_sigma")) { # Or adjust sigma
  stopifnot(is.matrix(X))                              # Ensure X is a matrix
  set.seed(seed)                                       # Set seed for reproducibility
  m <- ncol(X)                                         # Number of markers
  b <- rep(0, m)                                       # Initialize beta vector with zeros
  if (n_causal > 0)                                    # If any causal markers are requested
    b[seq_len(min(n_causal, m))] <- beta               # Set first n_causal betas to provided effect size
  eta0 <- as.numeric(X %*% b)                          # Compute linear predictor with initial betas
  var_eta0 <- stats::var(eta0)                         # Variance of linear predictor
  
  if (!is.null(pve_target) &&                          # Apply PVE control only if target provided
      is.finite(pve_target) && pve_target > 0 && pve_target < 1 &&
      n_causal > 0 && var_eta0 > 0) {                  # Ensure conditions allow scaling
    pve_mode <- match.arg(pve_mode)                    # Normalize PVE mode argument
    if (pve_mode == "scale_beta") {                    # If scaling beta to hit PVE with fixed sigma
      c_scale <- sqrt((pve_target / (1 - pve_target)) *# Compute scaling constant for betas
                        (error_sd^2 / var_eta0))        # Based on sigma and current eta variance
      b <- b * c_scale                                 # Scale betas
    } else {                                           # Otherwise adjust sigma to hit PVE with fixed betas
      error_sd <- sqrt(var_eta0 * (1 - pve_target) /   # Recompute residual sigma from target PVE
                         pve_target)                   # Using variance decomposition
    }
  }
  
  eta <- as.numeric(X %*% b)                           # Recompute linear predictor with final betas
  y <- eta + stats::rnorm(nrow(X), 0, error_sd)        # Simulate y as eta plus Gaussian noise
  realized_pve <- if ((stats::var(y) > 0))             # Compute realized sample PVE if var(y) > 0
    stats::var(eta) / stats::var(y)                    # PVE = Var(eta) / Var(y)
  else
    NA_real_                                           # Otherwise set as NA
  attr(y, "beta_vec") <- b                             # Attach beta vector as attribute
  attr(y, "sigma") <- error_sd                         # Attach sigma used as attribute
  attr(y, "pve") <- realized_pve                       # Attach realized PVE as attribute
  y                                                    # Return simulated y
}

# -----------------------------
# SIM & MMM statistics helpers
# -----------------------------

compute_sim_stats <- function(dat) {                   # Compute per-marker single-marker regression stats
  stopifnot("y" %in% names(dat))                       # Require column y
  x_cols <- names(dat)[startsWith(names(dat), "X_")]   # Identify X_* columns
  n <- nrow(dat)                                       # Number of rows (individuals)
  if (!length(x_cols))                                 # If no marker columns
    return(data.frame())                               # Return empty data frame
  
  fit0 <- lm(y ~ 1, data = dat)                        # Null model: intercept only
  rss0 <- sum(residuals(fit0)^2)                       # Residual sum of squares for null model
  
  out <- lapply(x_cols, function(v) {                  # For each marker column v
    fit1 <- lm(reformulate(v, "y"), data = dat)        # Full model: y ~ X_v
    co <- summary(fit1)$coefficients                   # Coefficient table
    rss1 <- sum(residuals(fit1)^2)                     # RSS for full model
    tval <- unname(co[2, 3])                           # t-statistic for marker coefficient
    pval_t <- unname(co[2, 4])                         # p-value for marker coefficient
    LRT  <- n * log(rss0 / rss1)                       # Likelihood ratio test statistic (approx, Gaussian)
    LOD  <- (n / 2) * log10(rss0 / rss1)               # LOD score
    c(                                                  # Return named numeric vector of stats
      beta_hat = unname(co[2, 1]),                     # Estimated beta
      t = tval,                                        # t-statistic
      p = pval_t,                                      # p-value
      RSS0 = rss0,                                     # Null RSS
      RSS1 = rss1,                                     # Full RSS
      LRT = LRT,                                       # LRT statistic
      LOD = LOD                                        # LOD score
    )
  })
  as.data.frame(do.call(rbind, out),                   # Bind results for all markers into a data frame
                row.names = x_cols,                    # Use marker names as rownames
                optional = TRUE)                       # Avoid converting strings to factors
}

compute_mmm_stats <- function(dat) {                   # Compute per-marker stats from multiple-marker model
  stopifnot("y" %in% names(dat))                       # Require y column
  x_cols <- names(dat)[startsWith(names(dat), "X_")]   # Identify X_* columns
  n <- nrow(dat)                                       # Number of observations
  p <- length(x_cols)                                  # Number of markers
  if (p == 0)                                          # If no predictors
    return(data.frame())                               # Return empty data frame
  fit_full <- lm(reformulate(x_cols, "y"), data = dat) # Full model including all markers
  rss_full <- sum(residuals(fit_full)^2)               # RSS for the full model
  sm <- summary(fit_full)$coefficients                 # Coefficient summary for full model
  
  out <- lapply(seq_along(x_cols), function(j) {       # Loop across markers to compute reduced-vs-full LRT
    v <- x_cols[j]                                     # Current marker name
    red <- if (p == 1)                                 # Build reduced model formula (drop v)
      "1"                                              # If only one predictor, reduced is intercept-only
    else
      paste(setdiff(x_cols, v), collapse = " + ")      # Otherwise include all other markers
    fit_red <- lm(as.formula(paste("y ~", red)), data = dat) # Fit reduced model without v
    rss_red <- sum(residuals(fit_red)^2)               # RSS for reduced model
    if (v %in% rownames(sm)) {                         # If coefficient row exists for v in full model
      beta_hat <- sm[v, 1]                             # Extract beta for v
      tval <- sm[v, 3]                                 # Extract t-statistic for v
      pval <- sm[v, 4]                                 # Extract p-value for v
    } else {                                           # If not present (e.g., collinearity)
      beta_hat <- NA                                   # Set beta to NA
      tval <- NA                                       # Set t to NA
      pval <- NA                                       # Set p to NA
    }
    LRT  <- n * log(rss_red / rss_full)                # LRT comparing reduced vs full
    LOD  <- (n / 2) * log10(rss_red / rss_full)        # Corresponding LOD score
    c(                                                  # Return named vector of results
      beta_hat = beta_hat,                             # Beta for v in full model
      t = tval,                                        # t-statistic
      p = pval,                                        # p-value
      RSS_reduced = rss_red,                           # RSS reduced
      RSS_full = rss_full,                             # RSS full
      LRT = LRT,                                       # LRT
      LOD = LOD                                        # LOD
    )
  })
  as.data.frame(do.call(rbind, out),                   # Bind results to a data frame
                row.names = x_cols,                    # Use marker names as rownames
                optional = TRUE)                       # Prevent factor conversion
}

# -----------------------------
# Small helpers for output report
# -----------------------------

`%||%` <- function(a, b) {                              # Null-or-empty coalescing operator
  if (is.null(a) ||                                     # Check for NULL
      (is.character(a) && length(a) == 1 && nchar(a) == 0)) { # Or empty single-length character
    b                                                   # Return fallback b
  } else {                                              # Otherwise
    a                                                   # Return a
  }
}

num_fmt <- function(x, digits = 3) {                    # Helper to format numeric values with significant digits
  if (is.numeric(x))                                    # Only format numerics
    return(signif(x, digits))                           # Return formatted numeric
  x                                                     # Return unchanged for non-numeric
}

df_to_table <- function(df,                             # Convert a small data frame to a basic HTML table
                        max_rows = 50,                  # Maximum rows to display
                        digits = 3) {                   # Digits for numeric formatting
  if (is.null(df) || !nrow(df))                         # If empty or NULL
    return(tags$em("No rows."))                         # Return an emphasis message
  df2 <- head(df, max_rows)                             # Take head to limit rows
  thead <- tags$thead(tags$tr(lapply(colnames(df2), tags$th))) # Build header row
  tbody <- tags$tbody(lapply(seq_len(nrow(df2)), function(i) { # Build table body
    tags$tr(lapply(df2[i, , drop = TRUE], function(val)        # For each cell in row i
      tags$td(as.character(num_fmt(val, digits)))))            # Create TD with formatted value
  }))
  tags$table(class = "table", thead, tbody)             # Return full table tag
}

plot_to_base64 <- function(plt_expr,                    # Render a ggplot object to a PNG and return as data URI
                           width = 900,                 # Image width in pixels
                           height = 420,                # Image height in pixels
                           res = 110) {                 # Resolution in DPI
  tf <- tempfile(fileext = ".png")                      # Create a temporary file path for the image
  png(tf, width = width, height = height, res = res)    # Open PNG device
  print(plt_expr)                                       # Print the ggplot to the device
  dev.off()                                             # Close the device to flush image
  on.exit(try(unlink(tf), silent = TRUE))               # Ensure the temp file is removed on exit
  base64enc::dataURI(file = tf, mime = "image/png")     # Return the file as a data URI string
}

# -----------------------------
# UI
# -----------------------------

ui <- navbarPage(                                       # Top-level UI: navbar with multiple tabs
  title = "QTL Design Sandbox",                         # App title shown in the navbar
  
  # ---- readme (first tab) ----
  tabPanel("readme", fluidPage(                         # README tab with usage instructions
    h2("How to use this app"),                          # Section title
    tags$ol(                                            # Ordered list of steps
      tags$li(                                          # Step for the data tab
        tags$b("data"), ": generate or upload data.",   # Short description
        tags$ul(                                        # Detailed bullets
          tags$li("Generate: configure n (sample size), m (markers), MAF, and seeds; click ", tags$code("Generate X"), "."), # X generation hint
          tags$li("Then set number of causal markers, effect size beta, error sigma, and optional PVE control; click ", tags$code("Generate y"), "."), # y generation hint
          tags$li("Upload: alternatively upload a CSV with columns ", tags$code("y"), ", ", tags$code("X_1..X_m"), ".") # Upload instructions
        )
      ),
      tags$li(                                          # Step for the model tab
        tags$b("model"), ": fit SIM (y ~ X_j) and MMM (y ~ all X_*).", # Model overview
        tags$ul(                                        # Details
          tags$li("SIM: per-marker regression returns coefficients, t, p, LRT, LOD."), # SIM details
          tags$li("MMM: full model with all markers; per-marker statistics from reduced-vs-full comparison."), # MMM details
          tags$li("Visualization thresholds: choose on p / LOD / LRT (p supports Bonferroni adjustment).") # Threshold info
        )
      ),
      tags$li(                                          # Step for the power tab
        tags$b("power"), ": plan sample size.",         # Power overview
        tags$ul(                                        # Details
          tags$li("Closure (lambda/NCP): specify alpha, beta_j, sigma, and VIF; compute power for given n and solve minimal n* for target power."), # Closure method
          tags$li("Bootstrap: choose model and marker; set n-grid and B to estimate rejection curve and read off threshold n.") # Bootstrap method
        )
      ),
      tags$li(                                          # Step for the output tab
        tags$b("output"), ": export an HTML report.",   # Output overview
        tags$ul(                                        # Details
          tags$li("Pick which sections to include (data / SIM / MMM / power / session)."), # Section selection
          tags$li("Click ", tags$code("Build report"), " to preview; then ", tags$code("Download HTML"), " for a self-contained file.") # Buttons usage
        )
      )
    ),
    tags$hr(),                                          # Horizontal rule separator
    tags$p(tags$small("Tip: compute results in the 'model' and 'power' tabs first; the 'output' tab pulls the latest results when building the report.")) # Helpful tip
  )),
  
  # ---- data ----
  tabPanel("data", sidebarLayout(                        # Data tab with controls and outputs
    sidebarPanel(                                       # Left sidebar for inputs
      h4("Active data source"),                         # Section heading
      radioButtons("active_source", "Use data from",    # Choice of generated vs uploaded dataset
                   choices = c("Generated" = "generated", "Uploaded" = "uploaded"),
                   selected = "generated"),
      tags$hr(),                                        # Divider
      
      h4("Generate clean X"),                           # Section for generating X
      numericInput("n_ind",  "Number of individuals", 200, min = 20, step = 20), # Number of rows
      numericInput("n_mark", "Number of markers",      50,  min = 5,  step = 5),  # Number of columns
      sliderInput("maf", "Minor allele frequency", min = 0.05, max = 0.5, value = 0.3, step = 0.01), # MAF slider
      numericInput("seed_x", "Seed for X", 1, step = 1),# Seed for X simulation
      actionButton("btn_gen_X", "Generate X"),          # Button to generate X
      tags$hr(),                                        # Divider
      
      h4("Generate y given X"),                         # Section for simulating y
      numericInput("n_causal", "Number of causal markers", 3, min = 0, step = 1), # Number of causal markers
      numericInput("beta", "Effect size per causal marker (beta)", 0.2, step = 0.05), # Beta per causal marker
      numericInput("error_sd", "Error SD (sigma)", 1, min = 1e-6, step = 0.05),  # Residual sigma
      numericInput("seed_y", "Seed for y", 1, step = 1),# Seed for y simulation
      
      selectInput("pve_mode", "PVE control mode",       # How to achieve target PVE
                  choices = c("No control" = "none",
                              "Fix sigma and scale beta (scale_beta)" = "scale_beta",
                              "Fix beta and solve sigma (set_sigma)"   = "set_sigma"),
                  selected = "none"),
      sliderInput("pve_target", "Target PVE (variance explained)", min = 0.01, max = 0.99, value = 0.5, step = 0.01), # Target PVE slider
      
      actionButton("btn_gen_Y", "Generate y"),          # Button to generate y
      tags$hr(),                                        # Divider
      
      h4("Or upload CSV (clean)"),                      # Upload section
      fileInput("upload", "Upload CSV", accept = ".csv"), # File input control
      helpText("Columns: y, X_1 ... X_m (no missing values). Uploaded and generated data remain available; the active source controls displays and fits.") # Help message
    ),
    mainPanel(                                          # Right main panel for outputs
      h4("Preview X (active source)"),                  # Heading for X preview
      DTOutput("x_head"),                               # Table: head of X
      tags$hr(),                                        # Divider
      h4("Preview data (y + X) — active source"),       # Heading for full data preview
      DTOutput("data_head"),                            # Table: head of full data
      tags$hr(),                                        # Divider
      h4("y histogram (active source)"),                # Heading for histogram
      plotOutput("y_hist", height = 240),               # Plot of y histogram
      tags$hr(),                                        # Divider
      h4("Realized PVE from simulation"),               # Heading for realized PVE
      uiOutput("pve_text"),                             # UI text for PVE
      tags$hr(),                                        # Divider
      verbatimTextOutput("data_info")                   # Text block with dataset info
    )
  )),
  
  # ---- model ----
  tabPanel("model", sidebarLayout(                       # Model tab containing SIM and MMM
    sidebarPanel(                                       # Left panel for model settings
      radioButtons("model_kind", "Model type",          # Choose SIM or MMM
                   c("Single-marker models (SIM: y ~ X_j)" = "single",
                     "Multiple-marker model (MMM: y ~ all X_*)" = "multi"),
                   selected = "single"),
      tags$hr(),                                        # Divider
      h4("Selection threshold"),                        # Threshold controls heading
      selectInput("thr_stat", "Statistic",              # Select which statistic to threshold on
                  choices = c("p-value" = "p", "LOD score" = "LOD", "LRT statistic" = "LRT"),
                  selected = "p"),
      numericInput("thr_value", "Threshold value", 0.05, min = 0, step = 0.01), # Numeric threshold value
      checkboxInput("thr_bonf", "Bonferroni across markers (for p-value)", FALSE), # Bonferroni toggle
      helpText("For p-value: selected if p <= threshold (after Bonferroni if enabled). For LOD/LRT: selected if statistic >= threshold."), # Help
      tags$hr(),                                        # Divider
      actionButton("btn_fit", "Fit model")              # Button to trigger model fitting
    ),
    mainPanel(                                          # Main panel with results
      conditionalPanel(                                 # Show SIM section when SIM is selected
        condition = "input.model_kind == 'single'",     # Condition expression
        h4("SIM results (per-marker)"),                 # SIM results heading
        DTOutput("single_res"),                         # Table of SIM stats
        uiOutput("sel_summary_single"),                 # Summary of selected markers
        tags$hr(),                                      # Divider
        h4("SIM visualization"),                        # Plot heading
        plotOutput("sim_vis", height = 260)             # Bar plot of chosen stat across markers
      ),
      conditionalPanel(                                 # Show MMM section when MMM is selected
        condition = "input.model_kind == 'multi'",      # Condition expression
        h4("MMM results (per-marker; reduced-vs-full LRT)"), # MMM results heading
        DTOutput("multi_res"),                          # Table of MMM stats
        uiOutput("sel_summary_multi"),                  # Summary of selected markers
        tags$hr(),                                      # Divider
        h4("MMM visualization"),                        # Plot heading
        plotOutput("mmm_vis", height = 260)             # Bar plot for MMM
      ),
      tags$hr(),                                        # Divider
      uiOutput("model_hint"),                           # Helpful notes depending on n and p
      uiOutput("model_status")                          # Staleness indicator for model outputs
    )
  )),
  
  # ---- power ----
  tabPanel("power", sidebarLayout(                       # Power tab containing closure and bootstrap methods
    sidebarPanel(                                       # Left panel for power settings
      h4("Choose method"),                              # Heading
      radioButtons("power_method", NULL,                # Choose closure or bootstrap
                   choices = c("Closure (lambda + bisection)" = "closure",
                               "Bootstrap (resampling)" = "bootstrap"),
                   selected = "closure"),
      tags$hr(),                                        # Divider
      conditionalPanel(                                 # Settings for closure method
        condition = "input.power_method == 'closure'",  # Show only for closure
        helpText("NCP/lambda method (single coefficient): lambda_j approx n * (beta_j^2/sigma^2) * 1/VIF_j; power via noncentral-t."), # Help
        numericInput("n_fixed", "n to check (given n)", 1000, min = 10, step = 10), # n for power evaluation
        numericInput("alpha_lambda", "Significance level (alpha)", 0.05, min = 1e-8, step = 0.01), # alpha
        numericInput("beta_j", "Assumed beta_j", 0.2, step = 0.01), # beta assumption
        numericInput("sigma_j", "Assumed sigma (error SD)", 1, min = 1e-8, step = 0.01), # sigma assumption
        checkboxInput("manual_vif", "Enter VIF manually", FALSE), # Manual VIF checkbox
        conditionalPanel(                                 # Manual VIF input block
          condition = "input.manual_vif == true",        # Show when manual VIF is checked
          numericInput("vif_manual", "VIF_j (manual)", 1, min = 1, step = 0.1) # Manual VIF value
        ),
        conditionalPanel(                                 # Auto VIF block
          condition = "input.manual_vif == false",       # Show when manual VIF is unchecked
          uiOutput("marker_choice_ui"),                  # Marker selector to compute VIF
          helpText("If data is available, VIF_j is computed from regressing X_j on X_-j; otherwise switch to manual VIF.") # Help
        ),
        numericInput("target_power_lambda", "Target power (default 0.90)", 0.90, min = 0.5, max = 0.999, step = 0.01), # Target power
        actionButton("btn_lambda", "Compute power / solve n (lambda method)") # Button to compute closure results
      ),
      conditionalPanel(                                 # Settings for bootstrap method
        condition = "input.power_method == 'bootstrap'",# Show only for bootstrap
        helpText("Bootstrap resampling: estimate rejection proportion across an n-grid."), # Help text
        selectInput("boot_model", "Model", choices = c("single","multi"), selected = "single"), # Model choice
        numericInput("alpha_boot", "Significance level (alpha)", 0.05, min = 1e-8, step = 0.01), # Alpha for bootstrap
        checkboxInput("bonf_boot", "Bonferroni across markers (multiple testing)", FALSE), # Bonf toggle
        uiOutput("marker_boot_ui"),                      # Marker selector for bootstrap
        numericInput("n_min_boot", "n grid: min", 100, min = 10, step = 10), # Grid min
        numericInput("n_max_boot", "n grid: max", 1000, min = 10, step = 10), # Grid max
        numericInput("n_step_boot", "n grid: step", 100, min = 1, step = 1),  # Grid step
        numericInput("B_boot", "Bootstrap repeats B", 200, min = 10, step = 10), # Number of bootstrap draws
        sliderInput("target_power_boot", "Target power", min = 0.5, max = 0.99, value = 0.9, step = 0.01), # Target power slider
        actionButton("btn_boot", "Run bootstrap power") # Execute bootstrap
      )
    ),
    mainPanel(                                          # Main panel for power results
      conditionalPanel(                                 # Closure results panel
        condition = "input.power_method == 'closure'",
        h4("Required sample size (lambda method)"),     # Heading for closure results
        infoBoxOutput("ibox_n", width = 12),            # Info box summarizing minimal n
        tags$hr(),                                      # Divider
        h4("NCP/lambda results"),                       # Heading for detailed text
        uiOutput("lambda_text")                         # Textual details for closure method
      ),
      conditionalPanel(                                 # Bootstrap results panel
        condition = "input.power_method == 'bootstrap'",
        h4("Bootstrap power results"),                  # Heading
        infoBoxOutput("ibox_boot", width = 12),         # Info box for bootstrap minimal n
        tags$hr(),                                      # Divider
        plotOutput("boot_plot", height = 300),          # Power curve plot
        DTOutput("boot_table")                          # Table of n vs power
      )
    )
  )),
  
  # ---- output (report) ----
  tabPanel("output", sidebarLayout(                      # Output tab: build and download HTML report
    sidebarPanel(                                       # Left panel for report options
      h4("Export HTML report"),                         # Heading
      textInput("report_title", "Title", "QTL Design Report"), # Report title input
      checkboxGroupInput("report_sections", "Include sections", # Sections to include
                         choices = c("data", "SIM", "MMM", "power", "session"),
                         selected = c("data", "SIM", "MMM", "power", "session")),
      numericInput("report_max_rows", "Max rows per table", 50, min = 5, step = 5), # Table row cap
      actionButton("btn_build_report", "Build report"), # Button to build preview
      tags$hr(),                                        # Divider
      downloadButton("btn_download_html", "Download HTML") # Button to download the report
    ),
    mainPanel(                                          # Right panel for preview
      uiOutput("report_preview_hint"),                  # Hint above iframe preview
      uiOutput("report_preview")                        # Iframe holding the report
    )
  ))
)

# -----------------------------
# SERVER
# -----------------------------

server <- function(input, output, session) {            # Server function defines all reactivity and outputs
  report_dir <- file.path(tempdir(), "qtl_reports")     # Directory under tempdir to serve generated reports
  dir.create(report_dir, showWarnings = FALSE, recursive = TRUE) # Ensure directory exists
  addResourcePath("qtl_reports", report_dir)            # Expose directory to the Shiny resource path
  
  X_gen   <- reactiveVal(NULL)                          # Reactive storage for generated X
  y_gen   <- reactiveVal(NULL)                          # Reactive storage for generated y
  df_gen  <- reactiveVal(NULL)                          # Reactive storage for generated full data (y + X)
  df_upl  <- reactiveVal(NULL)                          # Reactive storage for uploaded data
  stale_model <- reactiveVal(TRUE)                      # Flag indicating model outputs are stale
  stale_power <- reactiveVal(TRUE)                      # Flag indicating power outputs are stale
  rule_of_thumb_n <- function(p) ceiling(10 * max(1, p))# Quick heuristic: n >= 10 * p for regression
  
  lambda_cache <- reactiveVal(NULL)                     # Cache closure (lambda) method results
  boot_cache   <- reactiveVal(NULL)                     # Cache bootstrap results
  sim_cache    <- reactiveVal(NULL)                     # Cache latest SIM table
  mmm_cache    <- reactiveVal(NULL)                     # Cache latest MMM table
  
  output$ibox_n <- renderInfoBox({                      # Default info box for closure results before computing
    shinydashboard::infoBox("Sample size", "—",         # Title and placeholder value
                            "Click 'Compute power / solve n' to compute.", # Subtitle text
                            icon = icon("calculator"),  # Calculator icon
                            color = "light-blue", fill = TRUE) # Styling
  })
  output$ibox_boot <- renderInfoBox({                   # Default info box for bootstrap before computing
    shinydashboard::infoBox("Bootstrap power", "—",     # Title and placeholder
                            "Click 'Run bootstrap power' to compute.", # Subtitle
                            icon = icon("chart-line"),  # Chart icon
                            color = "light-blue", fill = TRUE) # Styling
  })
  
  use_generated <- function(){                          # Helper to switch active source to generated
    updateRadioButtons(session,"active_source",selected="generated") # Update radio selection
    stale_model(TRUE); stale_power(TRUE)                # Mark outputs as stale
  }
  use_uploaded  <- function(){                          # Helper to switch active source to uploaded
    updateRadioButtons(session,"active_source",selected="uploaded")  # Update radio selection
    stale_model(TRUE); stale_power(TRUE)                # Mark outputs as stale
  }
  observeEvent(input$active_source, {                   # When user switches source manually
    stale_model(TRUE); stale_power(TRUE)                # Mark outputs as stale
  }, ignoreInit = TRUE)                                 # Ignore initial trigger
  
  observeEvent(input$upload, {                          # When a CSV is uploaded
    file <- input$upload; req(file)                     # Require a file to be present
    dat <- read.csv(file$datapath, check.names = FALSE) # Read CSV with original column names
    validate(need("y" %in% names(dat), "Uploaded CSV must contain column 'y'.")) # Must include y
    validate(need(any(startsWith(names(dat), "X_")),    # Must include X_* columns
                  "Uploaded CSV must contain X_* columns."))
    df_upl(dat); use_uploaded()                         # Store uploaded data and switch active source
  }, ignoreInit = TRUE)                                 # Do not run on startup
  
  observeEvent(input$btn_gen_X, {                       # When user clicks Generate X
    G <- sim_X_clean(n = input$n_ind,                   # Call simulator with UI parameters
                     m = input$n_mark,
                     maf = input$maf,
                     seed = input$seed_x)
    X_gen(G); y_gen(NULL); df_gen(NULL); use_generated()# Store X, clear y/data, mark generated as active
  }, ignoreInit = TRUE)                                 # Ignore initial
  
  observeEvent(input$btn_gen_Y, {                       # When user clicks Generate y
    req(X_gen())                                        # Require that X exists
    pve_target <- if (identical(input$pve_mode,"none")) NULL else input$pve_target # Determine target PVE
    pve_mode   <- if (identical(input$pve_mode,"none")) "scale_beta" else input$pve_mode # Choose mode (unused if NULL)
    y <- sim_y_given_X(X_gen(),                         # Simulate y with current settings
                       n_causal = input$n_causal,
                       beta = input$beta,
                       error_sd = input$error_sd,
                       seed = input$seed_y,
                       pve_target = pve_target,
                       pve_mode = pve_mode)
    y_gen(y)                                            # Store simulated y
    df_gen(data.frame(y = as.numeric(y), X_gen(), check.names = FALSE)) # Build and store full data frame
    use_generated()                                     # Switch active source to generated
  }, ignoreInit = TRUE)                                 # Ignore initial
  
  current_data <- reactive({                            # Reactive getter for the active full dataset
    if (input$active_source == "uploaded"  && !is.null(df_upl())) return(df_upl()) # Return uploaded if active
    if (input$active_source == "generated" && !is.null(df_gen())) return(df_gen()) # Return generated if active
    NULL                                                # Otherwise return NULL
  })
  current_X_only <- reactive({                          # Reactive getter for X-only matrix for active source
    if (input$active_source == "generated" && !is.null(X_gen())) return(X_gen())   # Generated X
    if (input$active_source == "uploaded" && !is.null(df_upl())) {                 # Uploaded case
      x_cols <- names(df_upl())[startsWith(names(df_upl()), "X_")]                 # Identify X_* columns
      if (length(x_cols) > 0) return(as.matrix(df_upl()[, x_cols, drop = FALSE]))  # Return as matrix
    }
    NULL                                                # Otherwise NULL
  })
  
  output$x_head <- renderDT({                           # Render preview of X matrix
    if (input$active_source == "generated") {           # If using generated source
      if (is.null(X_gen()))                             # If X not yet generated
        return(datatable(data.frame(Info = "No generated X. Click 'Generate X'."), options = list(dom = 't'))) # Show message
      datatable(head(as.data.frame(X_gen()), 10), options = list(scrollX = TRUE, pageLength = 10)) # Show first 10 rows
    } else {                                            # If using uploaded source
      if (is.null(df_upl()))                            # If no upload yet
        return(datatable(data.frame(Info = "No uploaded data. Upload a CSV."), options = list(dom = 't'))) # Message
      x_cols <- names(df_upl())[startsWith(names(df_upl()), "X_")]                 # Identify X_* columns
      if (!length(x_cols))                              # If none present
        return(datatable(data.frame(Info = "Uploaded data has no X_* columns."), options = list(dom = 't'))) # Message
      datatable(head(df_upl()[, x_cols, drop = FALSE], 10), options = list(scrollX = TRUE, pageLength = 10)) # Show head
    }
  })
  output$data_head <- renderDT({                        # Render preview of full dataset
    dat <- current_data()                               # Get active data
    if (is.null(dat))                                   # If none
      return(datatable(data.frame(Info = "No full data (y + X). Generate y or upload CSV."), options = list(dom = 't'))) # Message
    datatable(head(dat, 10), options = list(scrollX = TRUE, pageLength = 10)) # Show head
  })
  output$y_hist <- renderPlot({                         # Render histogram of y
    dat <- current_data()                               # Get data
    if (is.null(dat) || !"y" %in% names(dat)) {         # If y not available
      plot.new(); title("No y yet (generate y or upload CSV including y).") # Show placeholder plot
      return(invisible())                               # Nothing else to do
    }
    hist(dat$y, breaks = 30, main = "Histogram of y (active source)", xlab = "y") # Draw histogram of y
  })
  output$pve_text <- renderUI({                         # Render realized or estimated PVE text
    dat <- current_data()                               # Get data
    if (is.null(dat) || !"y" %in% names(dat))           # If no y
      return(tags$em("No y."))                          # Show message
    pve_attr <- attr(y_gen(), "pve")                    # Realized PVE attribute from simulated y
    sigma_attr <- attr(y_gen(), "sigma")                # Sigma attribute from simulated y
    if (!is.null(pve_attr) && is.finite(pve_attr)) {    # If realized PVE available
      HTML(sprintf("Realized sample PVE approx <b>%.3f</b> (sigma approx %.3f).",
                   pve_attr, ifelse(is.null(sigma_attr), NA_real_, sigma_attr))) # Show realized PVE text
    } else {                                            # Otherwise estimate a quick R^2 using first few X
      x_cols <- names(dat)[startsWith(names(dat), "X_")]# Identify X_* columns
      if (!length(x_cols)) return(tags$em("No X_* to estimate PVE.")) # Show message if none
      q <- min(3L, length(x_cols))                      # Use first up to 3 markers
      qnames <- paste0("X_", seq_len(q))                # Build marker names
      fit <- try(lm(reformulate(qnames, "y"), data = dat), silent = TRUE) # Fit quick model
      if (inherits(fit, "try-error")) return(tags$em("Could not estimate PVE.")) # If failed, message
      r2 <- summary(fit)$r.squared                      # Extract R^2 as proxy
      HTML(sprintf("Estimated PVE (via R^2 of y ~ X_1..X_%d) approx <b>%.3f</b>.", q, r2)) # Display estimate
    }
  })
  output$data_info <- renderPrint({                     # Render dataset info text
    src <- input$active_source                          # Read active source selection
    if (src == "generated") {                           # If generated
      if (is.null(X_gen()) && is.null(df_gen())) {      # If neither X nor full data yet
        cat("Active source: Generated - No X or y yet.\n"); return(invisible()) # Print message and exit
      }
      if (!is.null(df_gen())) {                         # If full data available
        dat <- df_gen()                                 # Get generated data
        p <- sum(startsWith(names(dat), "X_"))          # Count markers
        cat("Active source: Generated\nRows:", nrow(dat), "  Markers:", p,
            "\nColumns:", paste(names(dat), collapse = ", "), "\n") # Print info
      } else {                                          # If only X is available
        cat("Active source: Generated (X only). Rows:", nrow(X_gen()), "  Markers:", ncol(X_gen()), "\n") # Print X-only info
      }
    } else {                                            # If uploaded
      if (is.null(df_upl())) {                          # If upload is absent
        cat("Active source: Uploaded - none.\n"); return(invisible()) # Print and exit
      }
      dat <- df_upl()                                   # Get uploaded data
      p <- sum(startsWith(names(dat), "X_"))            # Count markers
      cat("Active source: Uploaded\nRows:", nrow(dat), "  Markers:", p,
          "\nColumns:", paste(names(dat), collapse = ", "), "\n") # Print info
    }
  })
  
  sim_stats  <- eventReactive(input$btn_fit, {          # Compute SIM stats when Fit model is clicked
    dat <- current_data()                               # Get active data
    validate(need(!is.null(dat), "No full data. Generate y or switch to an uploaded dataset with y.")) # Ensure data exists
    validate(need("y" %in% names(dat), "Missing column 'y'.")) # Ensure y exists
    validate(need(any(startsWith(names(dat), "X_")), "No X_* columns found.")) # Ensure predictors exist
    stale_model(FALSE)                                  # Mark model outputs as fresh
    res <- compute_sim_stats(dat)                       # Compute SIM stats
    sim_cache(res)                                      # Cache for report
    res                                                 # Return result
  }, ignoreInit = TRUE)                                 # Ignore initial trigger
  
  mmm_stats <- eventReactive(input$btn_fit, {           # Compute MMM stats when Fit model is clicked
    dat <- current_data()                               # Get data
    validate(need(!is.null(dat), "No full data. Generate y or switch to uploaded data.")) # Ensure data exists
    validate(need(any(startsWith(names(dat), "X_")), "No X_* columns to fit.")) # Ensure X exists
    validate(need(nrow(dat) > sum(startsWith(names(dat), "X_")) + 1, "MMM requires n > p + 1. Increase n or reduce p.")) # Ensure n > p + 1
    stale_model(FALSE)                                  # Mark as fresh
    res <- compute_mmm_stats(dat)                       # Compute MMM stats
    mmm_cache(res)                                      # Cache for report
    res                                                 # Return result
  }, ignoreInit = TRUE)                                 # Ignore initial trigger
  
  apply_threshold <- function(df,                       # Helper to mark selected markers based on threshold
                              stat = c("p","LOD","LRT"),
                              thr = 0.05,
                              bonf = FALSE) {
    stat <- match.arg(stat)                             # Normalize stat argument
    if (nrow(df) == 0) return(df)                       # If empty, return as-is
    p <- nrow(df)                                       # Number of markers
    if (stat == "p") {                                  # p-value thresholding
      a <- if (bonf) thr / max(1L, p) else thr          # Effective alpha under Bonferroni if requested
      df$selected <- is.finite(df$p) & (df$p <= a)      # Selected if p <= a
      attr(df, "alpha_eff") <- a                        # Attach effective alpha
    } else if (stat == "LOD") {                         # LOD thresholding
      df$selected <- is.finite(df$LOD) & (df$LOD >= thr)# Selected if LOD >= thr
    } else {                                            # LRT thresholding
      df$selected <- is.finite(df$LRT) & (df$LRT >= thr)# Selected if LRT >= thr
    }
    df                                                 # Return modified data frame
  }
  
  output$single_res <- renderDT({                       # Render SIM results table
    df <- sim_stats(); req(df)                          # Ensure SIM results exist
    df2 <- df                                          # Copy for formatting
    if ("p" %in% names(df2))        df2$p        <- signif(df2$p, 4)   # Format p-values
    if ("beta_hat" %in% names(df2)) df2$beta_hat <- signif(df2$beta_hat, 4) # Format betas
    if ("t" %in% names(df2))        df2$t        <- round(df2$t, 3)    # Format t-stat
    if ("LRT" %in% names(df2))      df2$LRT      <- round(df2$LRT, 3)  # Format LRT
    if ("LOD" %in% names(df2))      df2$LOD      <- round(df2$LOD, 3)  # Format LOD
    df2 <- apply_threshold(df2, stat = input$thr_stat, thr = input$thr_value, bonf = isTRUE(input$thr_bonf)) # Apply threshold
    df2 <- cbind(marker = rownames(df2), df2)          # Add marker name as a column
    datatable(df2, options = list(pageLength = 10))    # Render as DataTable
  })
  output$sel_summary_single <- renderUI({               # Render summary text for SIM selections
    df <- sim_stats(); req(df)                          # Ensure data exists
    df2 <- apply_threshold(df, stat = input$thr_stat, thr = input$thr_value, bonf = isTRUE(input$thr_bonf)) # Apply threshold
    k  <- sum(isTRUE(df2$selected), na.rm = TRUE)       # Count selected markers
    a_eff <- attr(df2, "alpha_eff")                     # Extract effective alpha if applicable
    if (input$thr_stat == "p" && is.finite(a_eff)) {    # If p-value threshold
      HTML(sprintf("Selected markers (SIM) = <b>%d</b> / %d (effective alpha = %.3g).", k, nrow(df2), a_eff)) # Message with alpha
    } else {                                            # Otherwise
      HTML(sprintf("Selected markers (SIM) = <b>%d</b> / %d.", k, nrow(df2))) # Plain message
    }
  })
  
  output$multi_res <- renderDT({                        # Render MMM results table
    df <- mmm_stats(); req(df)                          # Ensure MMM exists
    df2 <- df                                          # Copy
    if ("p" %in% names(df2))        df2$p        <- signif(df2$p, 4)   # Format p
    if ("beta_hat" %in% names(df2)) df2$beta_hat <- signif(df2$beta_hat, 4) # Format beta
    if ("t" %in% names(df2))        df2$t        <- round(df2$t, 3)    # Format t
    if ("LRT" %in% names(df2))      df2$LRT      <- round(df2$LRT, 3)  # Format LRT
    if ("LOD" %in% names(df2))      df2$LOD      <- round(df2$LOD, 3)  # Format LOD
    df2 <- apply_threshold(df2, stat = input$thr_stat, thr = input$thr_value, bonf = isTRUE(input$thr_bonf)) # Apply threshold
    df2 <- cbind(marker = rownames(df2), df2)          # Add marker column
    datatable(df2, options = list(pageLength = 10))    # Render
  })
  output$sel_summary_multi <- renderUI({                # Render summary text for MMM selections
    df <- mmm_stats(); req(df)                          # Ensure data exists
    df2 <- apply_threshold(df, stat = input$thr_stat, thr = input$thr_value, bonf = isTRUE(input$thr_bonf)) # Apply threshold
    k  <- sum(isTRUE(df2$selected), na.rm = TRUE)       # Count selected
    a_eff <- attr(df2, "alpha_eff")                     # Effective alpha
    if (input$thr_stat == "p" && is.finite(a_eff)) {    # If p-statistic
      HTML(sprintf("Selected markers (MMM) = <b>%d</b> / %d (effective alpha = %.3g).", k, nrow(df2), a_eff)) # Message with alpha
    } else {                                            # Otherwise
      HTML(sprintf("Selected markers (MMM) = <b>%d</b> / %d.", k, nrow(df2))) # Plain message
    }
  })
  
  plot_stat <- function(df,                             # Generic plotting helper for per-marker statistic
                        stat,                           # Which statistic to plot
                        thr,                            # Threshold value to draw as dashed line
                        bonf = FALSE,                   # Whether to use Bonferroni for p-values
                        title_prefix = "") {            # Prefix for the plot title
    df$marker <- factor(rownames(df), levels = rownames(df)) # Keep marker order as-is
    if (stat == "p") {                                  # If plotting p-values
      a_eff <- if (bonf) thr / max(1L, nrow(df)) else thr # Effective alpha
      df$yval <- -log10(df$p)                           # Use -log10(p) scale
      ylab <- "-log10(p)"                               # Y-axis label
      thr_line <- -log10(a_eff)                         # Threshold line in transformed scale
      sub <- sprintf("Dashed line at -log10(alpha%s)=%.3f",
                     if (bonf) " (Bonferroni)" else "", thr_line) # Subtitle text
    } else if (stat == "LOD") {                         # If plotting LOD
      df$yval <- df$LOD                                 # Use LOD values
      ylab <- "LOD"                                     # Label
      thr_line <- thr                                   # Threshold line
      sub <- sprintf("Dashed line at LOD = %.3f", thr)  # Subtitle
    } else {                                            # If plotting LRT
      df$yval <- df$LRT                                 # Use LRT values
      ylab <- "LRT"                                     # Label
      thr_line <- thr                                   # Threshold line
      sub <- sprintf("Dashed line at LRT = %.3f", thr)  # Subtitle
    }
    df <- df[is.finite(df$yval), , drop = FALSE]        # Keep only finite values
    if (!nrow(df)) { plot.new(); title("No finite values to plot"); return(invisible()) } # Guard for empty
    ggplot(df, aes(x = marker, y = yval)) +             # Build bar chart
      geom_col() +                                      # Columns for each marker
      geom_hline(yintercept = thr_line, linetype = "dashed") + # Horizontal dashed threshold line
      labs(x = "Marker", y = ylab,                      # Axis labels
           title = sprintf("%s %s", title_prefix, switch(stat, p = "p-value", LOD = "LOD", LRT = "LRT")), # Title
           subtitle = sub) +                            # Subtitle showing threshold
      theme_minimal() +                                 # Minimal theme
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) # Rotate x labels for readability
  }
  
  output$sim_vis <- renderPlot({                        # Render SIM plot
    df <- sim_stats(); req(df)                          # Ensure SIM results exist
    plot_stat(df, stat = input$thr_stat, thr = input$thr_value, bonf = isTRUE(input$thr_bonf), title_prefix = "SIM ·") # Plot
  })
  
  output$mmm_vis <- renderPlot({                        # Render MMM plot
    df <- mmm_stats(); req(df)                          # Ensure MMM results exist
    plot_stat(df, stat = input$thr_stat, thr = input$thr_value, bonf = isTRUE(input$thr_bonf), title_prefix = "MMM ·") # Plot
  })
  
  output$model_hint <- renderUI({                       # Render helpful hints about scales and formulas
    if (is.null(current_data()))                        # If no data yet
      return(tagList(h4("Scale check / intuition"), p("No data yet."))) # Show message
    dat <- current_data()                               # Get data
    p <- sum(startsWith(names(dat), "X_"))             # Count markers
    n <- nrow(dat)                                      # Count individuals
    tagList(                                            # Return a list of UI elements
      h4("Scale check / intuition"),                    # Heading
      p(sprintf("Current: n = %d, p = %d.", n, p)),     # Show current n and p
      p(sprintf("Rule-of-thumb: n >= 10 * p approx %d (for linear regression).", rule_of_thumb_n(p))), # Heuristic
      p("SIM: per-marker null (y ~ 1) vs full (y ~ X_j): report t, p, LRT = n * log(RSS0/RSS1), LOD = (n/2) * log10(RSS0/RSS1)."), # SIM note
      p("MMM: full (y ~ all X) vs reduced (drop X_j): report coef t/p plus LRT_j, LOD_j.") # MMM note
    )
  })
  output$model_status <- renderUI({                     # Show staleness warning for model outputs
    if (isTRUE(stale_model()))
      tags$em("Model outputs are outdated due to data changes. Click 'Fit model' to refresh.")
    else
      NULL
  })
  
  compute_vif_j <- function(X, j){                      # Compute VIF for marker j by regressing X_j on other X
    p <- ncol(X); if (p <= 1) return(1)                 # If only one predictor, VIF is 1
    others <- setdiff(seq_len(p), j)                    # Indices of other predictors
    if (nrow(X) <= length(others) + 1) return(NA_real_) # Not enough rows to regress
    df_tmp <- data.frame(xj = X[, j], X[, others, drop = FALSE]) # Build temporary data frame
    fit <- try(lm(xj ~ ., data = df_tmp), silent = TRUE) # Fit regression
    if (inherits(fit, "try-error")) return(NA_real_)    # If failed, return NA
    r2 <- max(0, min(1, summary(fit)$r.squared))        # Clamp R^2 to [0,1]
    1 / max(1e-12, 1 - r2)                              # VIF = 1 / (1 - R^2), guard against zero
  }
  lambda_power_for_n <- function(n, beta, sigma, vif, p, alpha){ # Compute power under noncentral-t for given n
    n <- as.integer(n); p <- as.integer(p)              # Cast to integers
    if (n <= p + 1) return(0)                           # No degrees of freedom if n too small
    nu <- n - p - 1                                     # Degrees of freedom for t-test
    lambda <- n * (beta^2 / (sigma^2)) * (1 / vif)      # Noncentrality parameter squared
    delta  <- sqrt(lambda)                              # Noncentrality parameter
    tcrit  <- qt(1 - alpha/2, df = nu)                  # Two-sided critical t value
    pow <- pt(-tcrit, df = nu, ncp = delta) + (1 - pt(tcrit, df = nu, ncp = delta)) # Two-sided power
    as.numeric(max(0, min(1, pow)))                     # Clamp to [0,1] and return numeric
  }
  solve_n_lambda_halving_then_bisect <- function(target_power, # Solve minimal n achieving target power
                                                 beta, sigma, vif, p, alpha,
                                                 n_init = 1000, tol = 1, max_iter = 60){
    hi <- max(as.integer(n_init), p + 2)                # Start upper bound at provided n or minimal valid n
    pow_hi <- lambda_power_for_n(hi, beta, sigma, vif, p, alpha) # Power at hi
    grow_guard <- 0                                     # Loop guard counter for growth
    while (pow_hi < target_power && hi < 1e6 && grow_guard < 30) { # Expand upper bound until reaching target
      hi <- as.integer(ceiling(hi * 1.5))               # Increase hi by factor 1.5
      pow_hi <- lambda_power_for_n(hi, beta, sigma, vif, p, alpha) # Recompute power
      grow_guard <- grow_guard + 1                      # Increment guard
    }
    if (pow_hi < target_power)                          # If target unreachable within bounds
      return(list(n_star = NA_integer_, lo = NA_integer_, hi = NA_integer_, pow_lo = NA_real_, pow_hi = NA_real_)) # Return NA result
    lo <- max(p + 2, floor(hi / 2))                     # Initialize lower bound
    pow_lo <- lambda_power_for_n(lo, beta, sigma, vif, p, alpha) # Power at lo
    shrink_guard <- 0                                   # Guard for shrinking phase
    while (pow_lo >= target_power && lo > (p + 2) && shrink_guard < 30) { # Try shrinking hi downward
      hi <- lo; pow_hi <- pow_lo                        # Move hi to lo
      lo <- max(p + 2, floor(lo / 2))                   # Halve lo
      pow_lo <- lambda_power_for_n(lo, beta, sigma, vif, p, alpha) # Recompute power at new lo
      shrink_guard <- shrink_guard + 1                  # Increment guard
    }
    if (!(pow_lo < target_power && pow_hi >= target_power)) # If bracket not valid
      return(list(n_star = hi, lo = lo, hi = hi, pow_lo = pow_lo, pow_hi = pow_hi)) # Fallback return with hi
    it <- 0                                             # Iteration counter for bisection
    while ((hi - lo) > tol && it < max_iter) {          # Bisection until tolerance met
      it  <- it + 1                                     # Increment iteration
      mid <- as.integer(floor((lo + hi) / 2))           # Midpoint
      pow_mid <- lambda_power_for_n(mid, beta, sigma, vif, p, alpha) # Power at mid
      if (pow_mid >= target_power) { hi <- mid; pow_hi <- pow_mid } else { lo <- mid; pow_lo <- pow_mid } # Update bounds
    }
    list(n_star = hi, lo = lo, hi = hi, pow_lo = pow_lo, pow_hi = pow_hi) # Return result with minimal hi
  }
  
  output$marker_choice_ui <- renderUI({                 # UI for choosing a marker j for VIF
    Xonly <- current_X_only()                           # Get X matrix
    if (is.null(Xonly)) return(tags$em("No X available; use manual VIF.")) # If none, show message
    xnames <- colnames(Xonly)                           # Marker names
    selectInput("marker_j", "Choose marker j", choices = xnames, selected = xnames[1]) # Dropdown of markers
  })
  
  observeEvent(input$btn_lambda, {                      # When user triggers closure computation
    isolate({                                           # Use isolate to read current values once
      Xonly <- current_X_only()                         # Get X matrix
      if (is.null(Xonly) && !isTRUE(input$manual_vif)) {# Require X unless manual VIF is chosen
        output$ibox_n <- renderInfoBox({                # Update info box with error
          shinydashboard::infoBox("Sample size (lambda)", "—",
                                  "No X available to compute VIF. Use manual VIF.",
                                  icon = icon("triangle-exclamation"), color = "red", fill = TRUE)
        })
        return()                                        # Exit early
      }
      if (!is.null(Xonly)) {                            # If X is present
        p <- ncol(Xonly)                                # Determine number of markers
      } else {                                          # If X is absent
        output$ibox_n <- renderInfoBox({                # Inform user that p is unknown
          shinydashboard::infoBox("Sample size (lambda)", "—",
                                  "Please provide X to determine p (number of markers).",
                                  icon = icon("triangle-exclamation"), color = "red", fill = TRUE)
        })
        return()                                        # Exit early
      }
      if (isTRUE(input$manual_vif)) {                   # If manual VIF selected
        vif_j <- max(1, as.numeric(input$vif_manual))   # Use provided VIF, lower-bounded by 1
        marker_lab <- "(manual VIF)"                    # Label for report
      } else {                                          # Otherwise compute VIF from data
        j <- match(input$marker_j, colnames(Xonly))     # Index of chosen marker
        vif_j <- compute_vif_j(Xonly, j)                # Compute VIF
        marker_lab <- if (is.finite(vif_j)) paste0(input$marker_j, sprintf(" (VIF = %.3f)", vif_j)) else "(VIF unavailable)" # Label text
      }
      if (!is.finite(vif_j)) {                          # If VIF not available
        output$ibox_n <- renderInfoBox({                # Show error in info box
          shinydashboard::infoBox("Sample size (lambda)", "—",
                                  "Could not compute VIF; try manual VIF.",
                                  icon = icon("triangle-exclamation"), color = "red", fill = TRUE)
        })
        return()                                        # Exit early
      }
      n_given <- as.integer(input$n_fixed)              # Read given n
      alpha <- input$alpha_lambda                       # Read alpha
      beta  <- input$beta_j                             # Read beta
      sigma <- input$sigma_j                            # Read sigma
      target <- input$target_power_lambda               # Read target power
      pow_n <- lambda_power_for_n(n_given, beta, sigma, vif_j, p, alpha) # Compute power at given n
      res   <- solve_n_lambda_halving_then_bisect(      # Solve for minimal n
        target, beta, sigma, vif_j, p, alpha, n_init = n_given, tol = 1, max_iter = 60
      )
      n_min <- res$n_star                               # Extract minimal n
      lo_fin <- res$lo                                  # Lower bracket
      hi_fin <- res$hi                                  # Upper bracket
      pow_lo <- res$pow_lo                              # Power at lower bracket
      pow_hi <- res$pow_hi                              # Power at upper bracket
      if (is.na(n_min)) {                               # If unattainable
        output$ibox_n <- renderInfoBox({                # Show error state
          shinydashboard::infoBox("Sample size (lambda)", "—",
                                  "Target power may be unattainable under current settings.",
                                  icon = icon("triangle-exclamation"), color = "red", fill = TRUE)
        })
      } else {                                          # If solvable
        meets <- if (n_given >= n_min) "meets target" else "below target" # Compare given n to minimal n
        output$ibox_n <- renderInfoBox({                # Report minimal n and power at given n
          shinydashboard::infoBox("Sample size (lambda)", paste0("n* approx ", n_min),
                                  sprintf("Given n = %d -> power approx %.3f (%s)", n_given, pow_n, meets),
                                  icon = icon("calculator"), color = "blue", fill = TRUE)
        })
      }
      stale_power(FALSE)                                # Mark power outputs as fresh
      lambda_cache(list(                                # Cache closure results for the report
        marker_label = marker_lab,
        p = p,
        alpha = alpha,
        beta = beta,
        sigma = sigma,
        vif = vif_j,
        n_given = n_given,
        power_given = pow_n,
        target = target,
        n_star = n_min,
        bracket = c(lo_fin, hi_fin),
        pow_lo = pow_lo,
        pow_hi = pow_hi
      ))
    })
  })
  output$lambda_text <- renderUI({                      # Render detailed lambda method text
    input$btn_lambda                                    # Depend on button to update afterwards
    isolate({                                           # Use isolate to avoid reactive churn
      Xonly <- current_X_only()                         # Get X matrix
      if (is.null(Xonly) && !isTRUE(input$manual_vif))  # If we need X but do not have it
        return(tags$em("No X available to compute VIF. Check 'Enter VIF manually'.")) # Message
      if (!is.null(Xonly))                              # If X exists
        p <- ncol(Xonly)                                # Number of markers
      else                                              # If X absent
        return(tags$em("Please provide X to determine p (number of markers).")) # Message
      if (isTRUE(input$manual_vif)) {                   # Manual VIF branch
        vif_j <- max(1, as.numeric(input$vif_manual))   # Use manual VIF
        marker_lab <- "(manual VIF)"                    # Label
      } else {                                          # Computed VIF branch
        j <- match(input$marker_j, colnames(Xonly))     # Index of marker
        if (is.na(j)) return(tags$em("Invalid marker selection.")) # Guard on index
        vif_j <- compute_vif_j(Xonly, j)                # Compute VIF
        if (!is.finite(vif_j)) return(tags$em("Could not compute VIF; use manual VIF.")) # Guard on VIF
        marker_lab <- paste0(input$marker_j, sprintf(" (VIF = %.3f)", vif_j)) # Label text
      }
      n_given <- as.integer(input$n_fixed)              # Read given n
      alpha <- input$alpha_lambda                       # Read alpha
      beta  <- input$beta_j                             # Read beta
      sigma <- input$sigma_j                            # Read sigma
      target <- input$target_power_lambda               # Read target power
      pow_n <- lambda_power_for_n(n_given, beta, sigma, vif_j, p, alpha) # Power at given n
      res   <- solve_n_lambda_halving_then_bisect(      # Solve for minimal n
        target, beta, sigma, vif_j, p, alpha, n_init = n_given, tol = 1, max_iter = 60
      )
      n_min <- res$n_star                               # Minimal n
      if (is.null(n_min) || is.na(n_min)) {             # If unattainable
        tagList(                                        # Compose UI elements
          p(HTML(sprintf("<b>Marker:</b> %s; <b>p</b> = %d; <b>alpha</b> = %.3g; <b>beta_j</b> = %.3f; <b>sigma</b> = %.3f.",
                         marker_lab, p, alpha, beta, sigma))), # Summary line
          p(HTML("Using formula: lambda_j approx n * (beta_j^2/sigma^2) * 1/VIF_j.")), # Formula note
          tags$em("Target power may be unattainable under current settings.") # Warning
        )
      } else {                                          # If attainable
        tagList(                                        # Compose details
          p(HTML(sprintf("<b>Marker:</b> %s; <b>p</b> = %d; <b>alpha</b> = %.3g; <b>beta_j</b> = %.3f; <b>sigma</b> = %.3f.",
                         marker_lab, p, alpha, beta, sigma))), # Summary line
          p(HTML("Using formula: lambda_j approx n * (beta_j^2/sigma^2) * 1/VIF_j.")), # Formula note
          h4(sprintf("Minimal n to reach target power (%.2f): n* approx %d", target, n_min)), # Minimal n result
          p(sprintf("Final bracket: [%d, %d] with power(lo) = %.3f, power(hi) = %.3f",
                    res$lo, res$hi, res$pow_lo, res$pow_hi)) # Bracket details
        )
      }
    })
  })
  
  output$marker_boot_ui <- renderUI({                   # Build marker selector for bootstrap
    dat <- current_data()                               # Get active data
    if (is.null(dat)) return(tags$em("No data yet (need y + X_*)")) # Message if none
    x_cols <- names(dat)[startsWith(names(dat), "X_")]  # Identify X_* columns
    if (!length(x_cols)) return(tags$em("No X_* columns in current data.")) # Message if none
    selectInput("marker_boot", "Marker to test", choices = x_cols, selected = x_cols[1]) # Dropdown with markers
  })
  bootstrap_power <- function(dat,                      # Bootstrap power function
                              model = c("single", "multi"), # Model choice
                              alpha = 0.05,            # Significance level
                              n_grid,                  # Grid of n values
                              B = 200,                 # Number of bootstrap resamples
                              marker,                  # Marker name to test
                              bonferroni = FALSE) {    # Whether to use Bonferroni across markers
    model <- match.arg(model)                           # Normalize model argument
    x_cols <- names(dat)[startsWith(names(dat), "X_")]  # Identify X_* columns
    if (length(x_cols) == 0) stop("No X_* columns.", call. = FALSE) # Guard on predictors
    if (!("y" %in% names(dat))) stop("Missing column 'y'.", call. = FALSE) # Guard on response
    if (is.null(marker) || !(marker %in% x_cols)) stop("Invalid marker selection.", call. = FALSE) # Guard on marker
    p <- length(x_cols)                                 # Number of markers
    alpha_eff <- if (bonferroni) alpha / max(1L, p) else alpha # Effective alpha
    pow <- numeric(length(n_grid))                      # Allocate power vector
    for (i in seq_along(n_grid)) {                      # Loop over grid of n
      n_i <- n_grid[i]                                  # Current n
      rej <- 0L                                         # Rejection counter
      for (b in seq_len(B)) {                           # Bootstrap replications
        idx <- sample.int(nrow(dat), size = n_i, replace = TRUE) # Sample with replacement
        d2  <- dat[idx, , drop = FALSE]                 # Bootstrap sample
        if (model == "single") {                        # Single-marker model
          fit <- try(lm(reformulate(marker, "y"), data = d2), silent = TRUE) # Fit y ~ X_marker
          if (!inherits(fit, "try-error")) {            # If fit succeeded
            co <- summary(fit)$coefficients             # Coefficients
            pval <- if (nrow(co) >= 2) co[2, 4] else NA_real_ # p-value for marker
            if (is.finite(pval) && pval < alpha_eff) rej <- rej + 1L # Count rejection
          }
        } else {                                        # Multiple-marker model
          fit <- try(lm(reformulate(x_cols, "y"), data = d2), silent = TRUE) # Fit y ~ all X
          if (!inherits(fit, "try-error")) {            # If fit succeeded
            sm <- summary(fit)$coefficients             # Coefficients summary
            if (marker %in% rownames(sm)) {             # If marker present
              pval <- sm[marker, 4]                     # Extract p-value
              if (is.finite(pval) && pval < alpha_eff) rej <- rej + 1L # Count rejection
            }
          }
        }
      }
      pow[i] <- rej / B                                 # Empirical power at n_i
    }
    data.frame(n = n_grid, power = pow)                 # Return data frame of results
  }
  observeEvent(input$btn_boot, {                        # When user triggers bootstrap power
    dat <- current_data()                               # Get active data
    if (is.null(dat)) {                                 # If no data
      output$ibox_boot <- renderInfoBox({               # Show error
        shinydashboard::infoBox("Bootstrap power", "—",
                                "No data. Generate y or upload CSV.",
                                icon = icon("triangle-exclamation"), color = "red", fill = TRUE)
      })
      return()                                          # Exit
    }
    x_cols <- names(dat)[startsWith(names(dat), "X_")]  # Identify X_* columns
    if (length(x_cols) == 0) {                          # If none
      output$ibox_boot <- renderInfoBox({               # Show error
        shinydashboard::infoBox("Bootstrap power", "—",
                                "No X_* columns in current data.",
                                icon = icon("triangle-exclamation"), color = "red", fill = TRUE)
      })
      return()                                          # Exit
    }
    if (!("y" %in% names(dat))) {                       # If y missing
      output$ibox_boot <- renderInfoBox({               # Show error
        shinydashboard::infoBox("Bootstrap power", "—",
                                "Current data has no 'y' column.",
                                icon = icon("triangle-exclamation"), color = "red", fill = TRUE)
      })
      return()                                          # Exit
    }
    req(input$marker_boot)                              # Require marker selection
    if (!(input$marker_boot %in% x_cols)) {             # If marker invalid
      output$ibox_boot <- renderInfoBox({               # Show error
        shinydashboard::infoBox("Bootstrap power", "—",
                                "Selected marker not found in X_* columns.",
                                icon = icon("triangle-exclamation"), color = "red", fill = TRUE)
      })
      return()                                          # Exit
    }
    n_min <- as.integer(input$n_min_boot)               # Read n grid min
    n_max <- as.integer(input$n_max_boot)               # Read n grid max
    n_step <- as.integer(input$n_step_boot)             # Read n grid step
    if (is.na(n_min) || is.na(n_max) || is.na(n_step) || n_max < n_min || n_step <= 0) { # Validate grid
      output$ibox_boot <- renderInfoBox({               # Show error
        shinydashboard::infoBox("Bootstrap power", "—",
                                "Invalid n-grid settings.",
                                icon = icon("triangle-exclamation"), color = "red", fill = TRUE)
      })
      return()                                          # Exit
    }
    n_grid <- seq(n_min, n_max, by = n_step)            # Build grid
    dfp <- try(bootstrap_power(                         # Try to compute bootstrap power
      dat,
      model = input$boot_model,
      alpha = input$alpha_boot,
      n_grid = n_grid,
      B = as.integer(input$B_boot),
      marker = input$marker_boot,
      bonferroni = isTRUE(input$bonf_boot)
    ), silent = TRUE)
    if (inherits(dfp, "try-error")) {                   # If computation failed
      msg <- conditionMessage(attr(dfp, "condition"))   # Extract error message
      output$ibox_boot <- renderInfoBox({               # Show error
        shinydashboard::infoBox("Bootstrap power", "—",
                                paste("Error:", msg),
                                icon = icon("triangle-exclamation"), color = "red", fill = TRUE)
      })
      return()                                          # Exit
    }
    target <- input$target_power_boot                   # Read target power
    n_star <- NA_integer_                               # Initialize minimal n
    for (i in seq_len(nrow(dfp)))                       # Scan results
      if (is.finite(dfp$power[i]) && dfp$power[i] >= target) { n_star <- dfp$n[i]; break } # First n achieving target
    if (is.na(n_star)) {                                # If not reached
      output$ibox_boot <- renderInfoBox({               # Show warning
        shinydashboard::infoBox("Bootstrap power", "—",
                                "Target power not reached within the n-grid.",
                                icon = icon("triangle-exclamation"), color = "red", fill = TRUE)
      })
    } else {                                            # If reached
      output$ibox_boot <- renderInfoBox({               # Show minimal n
        shinydashboard::infoBox("Bootstrap power", paste0("n* approx ", n_star),
                                sprintf("Target power %.2f reached within grid.", target),
                                icon = icon("chart-line"), color = "blue", fill = TRUE)
      })
    }
    output$boot_plot <- renderPlot({                    # Render power curve plot
      ggplot(dfp, aes(n, power)) + geom_line() + geom_point() +
        geom_hline(yintercept = target, linetype = "dashed") +
        labs(y = "Power (reject proportion)", x = "n",
             title = sprintf("Bootstrap power (%s; marker = %s; B = %d)",
                             input$boot_model, input$marker_boot, input$B_boot)) +
        theme_minimal()
    })
    output$boot_table <- renderDT({                     # Render power table
      datatable(transform(dfp, power = round(power, 3)), options = list(pageLength = 10))
    })
    boot_cache(list(                                    # Cache bootstrap results for the report
      data = dfp,
      target = target,
      model = input$boot_model,
      marker = input$marker_boot,
      B = as.integer(input$B_boot)
    ))
  })
  
  report_path <- reactiveVal(NULL)                      # Reactive to hold the last generated HTML file path
  
  observeEvent(input$btn_build_report, {                # Build report when button is clicked
    dat <- current_data()                               # Get active data
    max_rows <- input$report_max_rows %||% 50           # Rows per table
    sections <- input$report_sections %||% character(0) # Sections to include
    title <- input$report_title %||% "QTL Design Report"# Report title
    now <- Sys.time()                                   # Timestamp
    
    css <- HTML(                                        # CSS for the report
      "body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif; margin: 24px; }
       h1,h2,h3 { margin-top: 1.2em; }
       .muted { color: #666; }
       .section { margin-bottom: 28px; }
       table.table { border-collapse: collapse; width: 100%; }
       table.table th, table.table td { border: 1px solid #ddd; padding: 6px 8px; }
       table.table th { background: #f6f8fa; text-align: left; }
       .code { background: #f6f8fa; padding: 12px; border-radius: 6px; white-space: pre-wrap; font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, 'Liberation Mono', 'Courier New', monospace; }
       .pill { display:inline-block; padding:2px 8px; border-radius:999px; background:#eef; margin-left:8px; font-size:12px; }"
    )
    
    blocks <- list(                                     # Initialize blocks with title and timestamp
      tags$h1(title),
      tags$p(class = "muted", sprintf("Generated: %s", format(now, "%Y-%m-%d %H:%M:%S")))
    )
    
    if ("data" %in% sections) {                         # Add data section if selected
      blk <- list(tags$div(
        class = "section",
        tags$h2("Data"),
        if (is.null(dat))                               # If no data
          tags$em("No data available. Generate or upload on the 'data' tab.") # Message
        else                                            # If data exists
          tagList(
            tags$p(sprintf("Rows: %d, Markers (p): %d", nrow(dat), sum(startsWith(names(dat), "X_")))), # Dimensions
            tags$h3("Preview (head)"),                  # Subheading
            df_to_table(dat, max_rows = min(max_rows, 10)), # Show small preview
            {                                           # Optional histogram
              if ("y" %in% names(dat)) {               # If y exists
                plt <- ggplot(dat, aes(x = y)) + geom_histogram(bins = 30) + theme_minimal() +
                  labs(title = "Histogram of y (active)", x = "y")           # Build histogram plot
                img <- plot_to_base64(plt)             # Encode as data URI
                tags$div(tags$img(src = img, style = "max-width:100%;height:auto;border:1px solid #eee")) # Embed image
              } else tags$em("No 'y' column to plot.") # Message if no y
            },
            {                                           # Realized PVE if available
              pve_attr <- attr(y_gen(), "pve")
              sigma_attr <- attr(y_gen(), "sigma")
              if (!is.null(pve_attr) && is.finite(pve_attr))
                tags$p(HTML(sprintf("Realized sample PVE approx <b>%.3f</b> (sigma approx %.3f).",
                                    pve_attr, ifelse(is.null(sigma_attr), NA_real_, sigma_attr))))
            }
          )
      ))
      blocks <- c(blocks, list(blk))                    # Append block
    }
    
    if ("SIM" %in% sections) {                          # Add SIM section if selected
      df_sim <- sim_cache()                             # Try cached SIM
      if (is.null(df_sim) && !is.null(dat) && "y" %in% names(dat) && any(startsWith(names(dat), "X_"))) {
        df_sim <- try(compute_sim_stats(dat), silent = TRUE) # Compute on the fly if needed
        if (inherits(df_sim, "try-error")) df_sim <- NULL    # Fallback to NULL on error
      }
      blk <- list(tags$div(
        class = "section",
        tags$h2("Single-marker (SIM) results"),
        if (is.null(df_sim) || !nrow(df_sim))           # If unavailable
          tags$em("Not computed yet. Please click 'Fit model' on the 'model' tab.") # Message
        else                                            # If available
          tagList(
            tags$p(HTML(sprintf("Threshold on <b>%s</b> = %.3g%s",
                                input$thr_stat, input$thr_value,
                                if (identical(input$thr_stat, "p") && isTRUE(input$thr_bonf)) " (Bonferroni)" else ""))), # Threshold summary
            df_to_table(cbind(marker = rownames(df_sim), df_sim), max_rows = max_rows), # Table
            {                                           # Plot
              plt <- plot_stat(df_sim, stat = input$thr_stat, thr = input$thr_value, bonf = isTRUE(input$thr_bonf), title_prefix = "SIM ·")
              img <- plot_to_base64(plt)
              tags$div(tags$img(src = img, style = "max-width:100%;height:auto;border:1px solid #eee"))
            }
          )
      ))
      blocks <- c(blocks, list(blk))                    # Append block
    }
    
    if ("MMM" %in% sections) {                          # Add MMM section if selected
      df_mmm <- mmm_cache()                             # Try cached MMM
      if (is.null(df_mmm) && !is.null(dat) && any(startsWith(names(dat), "X_")) &&
          nrow(dat) > sum(startsWith(names(dat), "X_")) + 1) {
        df_mmm <- try(compute_mmm_stats(dat), silent = TRUE) # Compute on the fly if possible
        if (inherits(df_mmm, "try-error")) df_mmm <- NULL    # Fallback to NULL on error
      }
      blk <- list(tags$div(
        class = "section",
        tags$h2("Multiple-marker (MMM) results"),
        if (is.null(df_mmm) || !nrow(df_mmm))           # If unavailable
          tags$em("Not computed or insufficient n > p + 1.") # Message
        else                                            # If available
          tagList(
            tags$p(HTML(sprintf("Threshold on <b>%s</b> = %.3g%s",
                                input$thr_stat, input$thr_value,
                                if (identical(input$thr_stat, "p") && isTRUE(input$thr_bonf)) " (Bonferroni)" else ""))), # Threshold summary
            df_to_table(cbind(marker = rownames(df_mmm), df_mmm), max_rows = max_rows), # Table
            {                                           # Plot
              plt <- plot_stat(df_mmm, stat = input$thr_stat, thr = input$thr_value, bonf = isTRUE(input$thr_bonf), title_prefix = "MMM ·")
              img <- plot_to_base64(plt)
              tags$div(tags$img(src = img, style = "max-width:100%;height:auto;border:1px solid #eee"))
            }
          )
      ))
      blocks <- c(blocks, list(blk))                    # Append block
    }
    
    if ("power" %in% sections) {                        # Add power section if selected
      lam <- lambda_cache()                              # Fetch lambda cache
      boot <- boot_cache()                               # Fetch bootstrap cache
      blk <- list(
        tags$div(
          class = "section",
          tags$h2("Power"),
          tags$h3("Lambda / NCP method"),               # Lambda subsection
          if (is.null(lam))                              # If no lambda results
            tags$em("Not computed yet. Use the 'closure' method on the 'power' tab.") # Message
          else                                           # Otherwise print details
            tagList(
              tags$p(HTML(sprintf("Marker: %s <span class='pill'>p = %d</span>", lam$marker_label, lam$p))), # Marker info
              tags$p(HTML(sprintf("alpha = %.3g; beta_j = %.3f; sigma = %.3f; VIF = %.3f", lam$alpha, lam$beta, lam$sigma, lam$vif))), # Params
              tags$p(sprintf("Given n = %d -> power approx %.3f", lam$n_given, lam$power_given)), # Power at given n
              if (is.na(lam$n_star)) tags$em("Target power unattainable under settings.") else
                tags$p(HTML(sprintf("Target power = %.2f -> minimal n*: <b>%d</b> (bracket [%d, %d])",
                                    lam$target, lam$n_star, lam$bracket[1], lam$bracket[2]))) # Minimal n
            ),
          tags$h3("Bootstrap method"),                   # Bootstrap subsection
          if (is.null(boot))                             # If no bootstrap results
            tags$em("Not computed yet. Run bootstrap on the 'power' tab.") # Message
          else                                           # Otherwise include plot and table
            tagList(
              tags$p(sprintf("Model = %s; Marker = %s; B = %d; target = %.2f",
                             boot$model, boot$marker, boot$B, boot$target)), # Settings
              df_to_table(transform(boot$data, power = round(power, 3)), max_rows = max_rows), # Table
              {                                         # Plot curve
                plt <- ggplot(boot$data, aes(n, power)) + geom_line() + geom_point() +
                  geom_hline(yintercept = boot$target, linetype = "dashed") +
                  labs(y = "Power", x = "n", title = sprintf("Bootstrap power (%s; %s)", boot$model, boot$marker)) +
                  theme_minimal()
                img <- plot_to_base64(plt)
                tags$div(tags$img(src = img, style = "max-width:100%;height:auto;border:1px solid #eee"))
              }
            )
        )
      )
      blocks <- c(blocks, list(blk))                    # Append block
    }
    
    if ("session" %in% sections) {                      # Add session info section if selected
      ses <- paste(capture.output(utils::sessionInfo()), collapse = "\n") # Capture sessionInfo text
      blk <- list(tags$div(
        class = "section",
        tags$h2("Session info"),
        tags$div(class = "code", ses)                   # Show as monospaced block
      ))
      blocks <- c(blocks, list(blk))                    # Append block
    }
    
    html <- tags$html(                                  # Construct full HTML document
      tags$head(
        tags$meta(charset = "utf-8"),                   # Charset meta
        tags$title(title),                              # Title tag
        tags$style(css)                                 # Inline CSS
      ),
      tags$body(do.call(tagList, blocks))               # Body with all blocks
    )
    
    tf <- tempfile(fileext = ".html")                   # Temp file path for the report
    htmltools::save_html(html, file = tf, background = "white") # Save HTML to file
    report_path(tf)                                     # Store temp path in reactive
    file.copy(tf, file.path(report_dir, basename(tf)), overwrite = TRUE) # Copy into served directory
    output$report_preview_hint <- renderUI(tags$p(class = "muted", "Preview below shows the generated HTML.")) # Hint text
    output$report_preview <- renderUI({                 # Build iframe preview
      tags$iframe(src = file.path("qtl_reports", basename(tf)), style = "width:100%;height:800px;border:1px solid #ddd;") # Iframe pointing to saved HTML
    })
  })
  
  output$btn_download_html <- downloadHandler(          # Download handler for the report
    filename = function() {                             # Provide download filename
      paste0((input$report_title %||% "QTL-Design-Report"), "_", format(Sys.time(), "%Y%m%d-%H%M%S"), ".html") # Title_timestamp.html
    },
    content = function(file) {                          # Content function to copy built file
      pth <- report_path()                              # Fetch generated report path
      validate(need(!is.null(pth) && file.exists(pth), "Please click 'Build report' first.")) # Ensure it exists
      file.copy(pth, file, overwrite = TRUE)            # Copy to download path
    }
  )
}

shinyApp(ui, server)                                    # Launch the Shiny application
