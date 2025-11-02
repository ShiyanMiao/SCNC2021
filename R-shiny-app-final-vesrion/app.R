# app.R — Shiny application for QTL design sandbox 

suppressPackageStartupMessages({
  library(shiny)               # Core Shiny
  library(DT)                  # DataTables
  library(shinydashboard)      # infoBox etc.
  library(ggplot2)             # plots
  library(htmltools)           # HTML utilities
  library(base64enc)           # base64 for images
  library(shinythemes)         # nice theme
  library(stats)               # explicit stats funcs
})

# -----------------------------
# Readme (with one-shot install/load summary)
# -----------------------------
readme_html <- htmltools::HTML('
<h1>QTL Design Sandbox</h1>
<h3>A Teaching Guide for Quantitative Trait Loci (QTL) Power Analysis in R Shiny</h3>
<hr />

<h2>Project Overview</h2>
<p>
QTL Design Sandbox is an interactive R Shiny application that lets you simulate quantitative traits,
fit marker-selection models, and evaluate statistical power and false-discovery rate (FDR) under different experimental designs.
It is developed as part of SCNC2021 – Research Project in Science at The Australian National University by Shiyan Miao (u8027892),
supervised by Dr Emi Tanaka and Dr Yidi Deng.
</p>

<h3>Background</h3>
<p>Marker-assisted selection (MAS) improves breeding efficiency by identifying genetic markers associated with target traits.
However, power to detect markers depends on factors such as sample size, effect sizes, population structure, and model choice.
This app provides a simulation-based teaching framework to explore these factors and support experimental design.</p>
<hr />

<h2>Learning Goals</h2>
<ul>
  <li>Understand how sample size, heritability (PVE), and model choice influence power in QTL studies.</li>
  <li>Practise simulation and statistical modelling using R/Shiny.</li>
  <li>Develop skills in reproducible research and visual analysis.</li>
</ul>
<hr />

<h2>App Structure & Required R Packages</h2>
<ol>
  <li><b>Readme Page</b> – Project info and session details.
    <pre class="code"><code># Required
install.packages(c("htmltools", "shinythemes"))
library(htmltools); library(shinythemes)</code></pre>
  </li>

  <li><b>Data Page</b> – Generate clean X, simulate y, or upload CSV; preview & QC.
    <pre class="code"><code># Required
install.packages(c("shiny", "DT"))
library(shiny); library(DT)</code></pre>
  </li>

  <li><b>Model Page</b> – Fit SIM (y ~ X_j) or MMM (y ~ all X_*); show t, p, LRT, LOD; plot with thresholds.
    <pre class="code"><code># Required
install.packages(c("stats", "ggplot2"))
library(stats); library(ggplot2)</code></pre>
  </li>

  <li><b>Power Page</b> – Closure/NCP (bisection) and Bootstrap power; FDR curve (BH).
    <pre class="code"><code># Required
install.packages(c("ggplot2"))
library(ggplot2)  # plotting power/FDR curves</code></pre>
  </li>

  <li><b>Output Page</b> – Build an HTML report that embeds tables/plots.
    <pre class="code"><code># Required
install.packages(c("htmltools", "base64enc"))
library(htmltools); library(base64enc)</code></pre>
  </li>

  <li><b>UI Elements</b> – Info boxes and theme.
    <pre class="code"><code># Required
install.packages(c("shinydashboard", "shinythemes"))
library(shinydashboard); library(shinythemes)</code></pre>
  </li>
</ol>

<h3>One-shot install & load</h3>
<p>Install once, then load each time you run the app.</p>
<pre class="code"><code># Install all (run once)
install.packages(c(
  "shiny",        # core app framework
  "shinythemes",  # themes
  "shinydashboard", # infoBox etc.
  "DT",           # tables
  "ggplot2",      # plots
  "htmltools",    # HTML helpers
  "base64enc"     # embed plots in report
  # stats is base R (no install needed)
))

# Load for current session
library(shiny)
library(shinythemes)
library(shinydashboard)
library(DT)
library(ggplot2)
library(htmltools)
library(base64enc)
library(stats)  # base R, for clarity
</code></pre>
<hr />

<h2>Step-by-Step Teaching Guide</h2>

<h3>1) Data Page</h3>
<p><i>Purpose:</i> Provide genotype data (X) and optionally simulate phenotype (y).</p>
<ul>
  <li><b>Upload CSV</b>: Must include <code>y</code> and <code>X_*</code> columns without missing values.</li>
  <li><b>Generate X</b>: Set <code>n</code>, <code>m</code>, <code>MAF</code>, and seed.</li>
  <li><b>Simulate y</b>: Choose number of causal markers, effect size, error SD; optionally target PVE.</li>
</ul>

<h3>2) Model Page</h3>
<p><i>Purpose:</i> Fit SIM (per-marker) and MMM (all markers) models; view estimates and selection by threshold.</p>
<ul>
  <li><b>Thresholding</b>: p/LOD/LRT with optional Bonferroni for p-values.</li>
  <li><b>Visualisation</b>: Bar plot with dashed threshold lines.</li>
</ul>

<h3>3) Power Page</h3>
<p><i>Purpose:</i> Estimate power vs. sample size using two complementary methods and evaluate FDR (BH).</p>

<h4>3.1 Closure/NCP (“bisection”) method — How it works</h4>
<ol>
  <li>For a target coefficient <code>beta_j</code>, error SD <code>sigma</code>, and collinearity via <code>VIF_j</code>, compute the noncentrality: <code>lambda_j ≈ n * (beta_j^2 / sigma^2) * (1 / VIF_j)</code>.</li>
  <li>With degrees of freedom <code>ν = n - p - 1</code>, approximate power from the noncentral t distribution at level <code>alpha</code>.</li>
  <li>Use a bracket-halving + bisection scheme to find the smallest <code>n*</code> where power ≥ target (e.g., 0.90).</li>
</ol>

<h4>3.2 Bootstrap method — How it works</h4>
<ol>
  <li>For each grid value of <code>n</code>, repeatedly resample (with replacement) rows from the current dataset.</li>
  <li>Refit the chosen model (SIM/MMM) each time; record whether the chosen marker passes the selection rule (e.g., p ≤ alpha with/without Bonferroni).</li>
  <li>Estimated power at that <code>n</code> is the rejection proportion across <code>B</code> resamples.</li>
</ol>

<h4>3.3 FDR curve (BH) — How it works</h4>
<ol>
  <li>For each resample and each marker, compute p-values and apply Benjamini–Hochberg (BH) to get q-values.</li>
  <li>Count selected markers at q ≤ <code>q_threshold</code>. Estimate false discoveries either by (a) known truth in simulation or (b) permutation of <code>y</code> (recommended when truth is unknown).</li>
  <li>Report E[FDP] vs. <code>n</code> with mean ± SD over <code>B</code> resamples.</li>
</ol>

<h4>3.4 Why bootstrap power can be smaller than bisection power</h4>
<ul>
  <li><b>Sampling variability:</b> Bootstrap reflects finite-sample noise and instability in refits; theory assumes ideal conditions.</li>
  <li><b>Model misspecification / collinearity:</b> Real X may inflate variance (VIF), lowering empirical detection compared to the NCP formula.</li>
  <li><b>Discrete thresholding:</b> With small-to-moderate n, estimates fluctuate around the cut-off; fewer resamples pass.</li>
  <li><b>Multiple testing adjustments:</b> Bonferroni/BH are applied inside bootstrap; the NCP calculation targets a single-coefficient test.</li>
</ul>
<p><i>In short:</i> Bootstrap estimates real-world power; the NCP/bisection method estimates theoretical power. A small gap is normal.</p>
<hr />

<h3>4) Output Page</h3>
<p><i>Purpose:</i> Export a self-contained HTML summary of data, models, power, and FDR results.</p>
<ul>
  <li>Tables are rendered from in-app objects; plots are embedded as base64 PNG.</li>
  <li>Choose sections to include and download the final report.</li>
</ul>

<hr />
<h2>Session Info</h2>
<p>The current R session details are shown below for reproducibility.</p>
')



# -----------------------------                                   # Section divider for readability
# Data generation — split X and y                                 # Title: functions to simulate X and y
# -----------------------------                                   # Section divider

sim_X_clean <- function(n = 200,                                  # Function: simulate a clean genotype matrix X
                        m = 50,                                    #   m = number of markers (columns)
                        maf = 0.3,                                 #   maf = minor allele frequency for Binomial(2, maf)
                        seed = 1) {                                #   seed = RNG seed for reproducibility
  set.seed(seed)                                                   # Fix random seed so results are reproducible
  G <- matrix(rbinom(n * m, 2, maf), nrow = n, ncol = m)          # Draw n*m genotypes and reshape to n×m matrix
  colnames(G) <- paste0("X_", seq_len(m))                          # Name columns X_1, X_2, …, X_m
  G                                                                 # Return the simulated matrix
}

sim_y_given_X <- function(X,                                       # Function: simulate phenotype y conditional on X
                          n_causal = 3,                             #   Number of causal markers (first n_causal columns)
                          beta = 0.2,                               #   Effect size assigned to each causal marker
                          error_sd = 1,                             #   Residual standard deviation (sigma)
                          seed = 1,                                 #   RNG seed for reproducibility
                          pve_target = NULL,                        #   Optional target PVE (proportion of variance explained)
                          pve_mode = c("scale_beta", "set_sigma")) {#   Strategy to hit PVE: scale betas or set sigma
  stopifnot(is.matrix(X))                                           # Guard: X must be a numeric matrix
  set.seed(seed)                                                    # Fix random seed
  m <- ncol(X)                                                      # Number of markers (columns) in X
  b <- rep(0, m)                                                    # Initialize effect vector with zeros
  if (n_causal > 0)                                                 # If at least one causal marker is requested
    b[seq_len(min(n_causal, m))] <- beta                            #   Assign beta to the first n_causal markers
  eta0 <- as.numeric(X %*% b)                                       # Linear predictor before any PVE adjustments
  var_eta0 <- stats::var(eta0)                                      # Variance of the linear predictor
  
  if (!is.null(pve_target) &&                                       # If a PVE target is provided and valid
      is.finite(pve_target) && pve_target > 0 && pve_target < 1 &&  #   PVE must be in (0,1)
      n_causal > 0 && var_eta0 > 0) {                               #   Need signal variance to be positive
    pve_mode <- match.arg(pve_mode)                                 # Normalize pve_mode argument
    if (pve_mode == "scale_beta") {                                 # Option 1: hold sigma, scale betas to hit PVE
      c_scale <- sqrt((pve_target / (1 - pve_target)) * (error_sd^2 / var_eta0)) # Scaling factor for betas
      b <- b * c_scale                                              # Apply scaling to effect vector
    } else {                                                        # Option 2: hold betas, solve sigma to hit PVE
      error_sd <- sqrt(var_eta0 * (1 - pve_target) / pve_target)    # Recompute residual SD from variance identity
    }
  }
  
  eta <- as.numeric(X %*% b)                                        # Final linear predictor with adjusted betas/sigma
  y <- eta + stats::rnorm(nrow(X), 0, error_sd)                     # Simulate phenotype: y = eta + Gaussian noise
  realized_pve <- if ((stats::var(y) > 0))                          # Compute realized sample PVE if Var(y) > 0
    stats::var(eta) / stats::var(y)                                 #   PVE = Var(eta) / Var(y)
  else                                                              # If Var(y) is zero (degenerate)
    NA_real_                                                         #   Mark PVE as NA
  attr(y, "beta_vec") <- b                                          # Attach used beta vector as attribute
  attr(y, "sigma") <- error_sd                                      # Attach used sigma as attribute
  attr(y, "pve") <- realized_pve                                    # Attach realized PVE as attribute
  y                                                                  # Return the simulated phenotype vector
}

# -----------------------------                                   # Section divider
# SIM & MMM statistics helpers                                     # Title: helpers to compute per-marker statistics
# -----------------------------                                   # Section divider

compute_sim_stats <- function(dat) {                                # Compute single-marker (SIM) stats for each X_j
  stopifnot("y" %in% names(dat))                                    # Guard: response column y must exist
  x_cols <- names(dat)[startsWith(names(dat), "X_")]                # Identify predictor columns named X_*
  n <- nrow(dat)                                                    # Number of observations (rows)
  if (!length(x_cols))                                              # If no X_* columns are present
    return(data.frame())                                            #   Return empty data frame
  
  fit0 <- lm(y ~ 1, data = dat)                                     # Null model: intercept only
  rss0 <- sum(residuals(fit0)^2)                                    # Null residual sum of squares
  
  out <- lapply(x_cols, function(v) {                                # For each marker column v
    fit1 <- lm(reformulate(v, "y"), data = dat)                     #   Fit full model: y ~ X_v
    co <- summary(fit1)$coefficients                                #   Extract coefficient table
    rss1 <- sum(residuals(fit1)^2)                                  #   Full-model RSS
    tval <- unname(co[2, 3])                                        #   t-statistic for X_v
    pval_t <- unname(co[2, 4])                                      #   p-value for X_v
    LRT  <- n * log(rss0 / rss1)                                    #   Likelihood ratio stat (Gaussian approx.)
    LOD  <- (n / 2) * log10(rss0 / rss1)                            #   LOD score corresponding to LRT
    c(                                                               #   Return a named vector of results
      beta_hat = unname(co[2, 1]),                                   #     Estimated effect of X_v
      t = tval,                                                      #     t-statistic
      p = pval_t,                                                    #     p-value
      RSS0 = rss0,                                                   #     Null RSS
      RSS1 = rss1,                                                   #     Full-model RSS
      LRT = LRT,                                                     #     Likelihood-ratio statistic
      LOD = LOD                                                      #     LOD score
    )
  })
  as.data.frame(do.call(rbind, out),                                 # Bind per-marker results into a data frame
                row.names = x_cols,                                  # Use marker names as row names
                optional = TRUE)                                     # Prevent strings-as-factors behavior
}

compute_mmm_stats <- function(dat) {                                 # Compute per-marker stats from multi-marker model
  stopifnot("y" %in% names(dat))                                     # Guard: y must be present
  x_cols <- names(dat)[startsWith(names(dat), "X_")]                 # Collect all X_* predictors
  n <- nrow(dat)                                                     # Number of observations
  p <- length(x_cols)                                                # Number of predictors
  if (p == 0)                                                        # If no predictors exist
    return(data.frame())                                             #   Return empty data frame
  fit_full <- lm(reformulate(x_cols, "y"), data = dat)               # Fit full model: y ~ all X_*
  rss_full <- sum(residuals(fit_full)^2)                             # RSS of full model
  sm <- summary(fit_full)$coefficients                               # Coefficient table of full model
  
  out <- lapply(seq_along(x_cols), function(j) {                     # For each predictor index j
    v <- x_cols[j]                                                   #   Current predictor name
    red <- if (p == 1)                                               #   Build reduced-model RHS (drop v)
      "1"                                                            #   If only one predictor, reduced is intercept-only
    else                                                             #   Otherwise include all predictors except v
      paste(setdiff(x_cols, v), collapse = " + ")
    fit_red <- lm(as.formula(paste("y ~", red)), data = dat)         #   Fit reduced model without v
    rss_red <- sum(residuals(fit_red)^2)                             #   RSS of reduced model
    if (v %in% rownames(sm)) {                                       #   If coefficient for v exists in full model
      beta_hat <- sm[v, 1]                                           #     Estimated coefficient
      tval <- sm[v, 3]                                               #     t-statistic
      pval <- sm[v, 4]                                               #     p-value
    } else {                                                         #   Else (e.g., due to collinearity)
      beta_hat <- NA                                                 #     Mark beta as NA
      tval <- NA                                                     #     Mark t as NA
      pval <- NA                                                     #     Mark p as NA
    }
    LRT  <- n * log(rss_red / rss_full)                              #   LRT comparing reduced vs full
    LOD  <- (n / 2) * log10(rss_red / rss_full)                      #   Corresponding LOD score
    c(                                                                #   Return named vector for this marker
      beta_hat = beta_hat,                                           #     Beta estimate in full model
      t = tval,                                                      #     t-statistic in full model
      p = pval,                                                      #     p-value in full model
      RSS_reduced = rss_red,                                         #     RSS of reduced model
      RSS_full = rss_full,                                           #     RSS of full model
      LRT = LRT,                                                     #     Likelihood-ratio statistic
      LOD = LOD                                                      #     LOD score
    )
  })
  as.data.frame(do.call(rbind, out),                                 # Bind per-marker outputs into a data frame
                row.names = x_cols,                                  # Use marker names as row names
                optional = TRUE)                                     # Keep consistent types (no factors)
}

# -----------------------------                                   # Section divider
# Small helpers for output report                                  # Title: utility helpers for HTML report
# -----------------------------                                   # Section divider

`%||%` <- function(a, b) {                                          # Null-or-empty coalescing operator
  if (is.null(a) ||                                                 #   If a is NULL
      (is.character(a) && length(a) == 1 && nchar(a) == 0))         #   or a is a length-1 empty string
    b                                                               #   return b (fallback)
  else                                                              #   otherwise
    a                                                               #   return a as-is
}
num_fmt <- function(x, digits = 3) {                                # Format numerics with significant digits
  if (is.numeric(x))                                                #   Only format numeric inputs
    return(signif(x, digits))                                       #   Apply signif() and return
  x                                                                  #   Non-numeric: return unchanged
}

df_to_table <- function(df,                                         # Convert a small data.frame to a simple HTML table
                        max_rows = 50,                              #   Maximum number of rows to display
                        digits = 3) {                               #   Significant digits for numeric cells
  if (is.null(df) || !nrow(df))                                     # If df is NULL or has zero rows
    return(tags$em("No rows."))                                     #   Return an emphasized placeholder
  df2 <- head(df, max_rows)                                         # Keep only the first max_rows rows
  thead <- tags$thead(tags$tr(lapply(colnames(df2), tags$th)))      # Build header row with column names
  tbody <- tags$tbody(lapply(seq_len(nrow(df2)), function(i) {      # Build table body by rows
    tags$tr(lapply(df2[i, , drop = TRUE], function(val)             #   For each cell in row i
      tags$td(as.character(                                         #     Create <td> with formatted text
        num_fmt(val, digits)                                        #       Apply numeric formatting if needed
      ))))
  }))
  tags$table(class = "table", thead, tbody)                         # Return complete <table> element
}

plot_to_base64 <- function(plt_expr,                                # Render a ggplot to PNG and return data URI
                           width = 900,                             #   Output image width in pixels
                           height = 420,                            #   Output image height in pixels
                           res = 110) {                             #   PNG resolution (DPI)
  tf <- tempfile(fileext = ".png")                                  # Create a temporary file path
  png(tf,                                                           # Open PNG device to the temp file
      width = width,                                                #   Set width
      height = height,                                              #   Set height
      res = res)                                                    #   Set resolution
  print(plt_expr)                                                   # Draw the ggplot onto the device
  dev.off()                                                         # Close the graphics device (flush image)
  on.exit(try(unlink(tf), silent = TRUE))                           # Ensure temp file is cleaned up on exit
  base64enc::dataURI(file = tf, mime = "image/png")                 # Return image as base64 data URI
}

# -----------------------------                                   # Section divider for readability
# UI                                                              # Top-level UI definition
# -----------------------------                                   # Section divider

ui <- navbarPage(                                                  # Create a navbar-based multi-tab layout
  title = "QTL Design Sandbox",                                    # App title shown in the navbar
  theme = shinythemes::shinytheme("flatly"),                       # Use the "flatly" Bootstrap theme
  
  # ---- readme tab ----
  tabPanel(                                                        # Define the first tab
    "readme",                                                      # Tab label in the navbar
    fluidPage(                                                     # Standard fluid layout for content
      tags$head(tags$style(                                        # Inject custom CSS into the page head
        HTML(                                                      # CSS rules written as a single HTML string
          "
        .readme-body h1,h2,h3 { margin-top: 1rem; }                # Add spacing above headers inside readme body
        .readme-body ul, .readme-body ol { margin-left: 1.1rem; }  # Indent lists a bit for readability
        .readme-body pre { background:#f6f8fa; padding:10px;       # Style code blocks with light background
                             border-radius:6px; }                  # Rounded corners for code blocks
      "
        )
      )),
      h2("Project Readme"),                                        # Section header inside the tab
      div(class = "readme-body", readme_html),                     # Render the HTML README content
      tags$hr(),                                                   # Horizontal rule separator
      h4("Session info"),                                          # Small header for session info
      verbatimTextOutput("session_info")                           # Plain-text output of sessionInfo()
    )
  ),
  
  
  # ---- data ----
  tabPanel("data", sidebarLayout(                                  # Second tab: data; use sidebar layout
    sidebarPanel(                                                  # Left panel with controls
      h4("Active data source"),                                    # Section header for source selection
      radioButtons(                                                # Radio inputs to choose data source
        "active_source",                                           # Input id
        "Use data from",                                           # Input label
        choices = c("Generated" = "generated", "Uploaded" = "uploaded"), # Two choices with values
        selected = "generated"                                     # Default selection
      ),
      tags$hr(),                                                   # Divider line
      h4("Generate clean X"),                                      # Header for X-generation controls
      numericInput(                                                # Input: number of individuals (rows)
        "n_ind",
        "Number of individuals",
        200,
        min = 20,
        step = 20
      ),
      numericInput(                                                # Input: number of markers (columns)
        "n_mark",
        "Number of markers",
        50,
        min = 5,
        step = 5
      ),
      sliderInput(                                                 # Slider: minor allele frequency
        "maf",
        "Minor allele frequency",
        min = 0.05,
        max = 0.5,
        value = 0.3,
        step = 0.01
      ),
      numericInput("seed_x", "Seed for X", 1, step = 1),           # Input: RNG seed for X
      actionButton("btn_gen_X", "Generate X"),                     # Button: trigger X generation
      tags$hr(),                                                   # Divider
      h4("Generate y given X"),                                    # Header for y-generation controls
      numericInput(                                                # Input: number of causal markers
        "n_causal",
        "Number of causal markers",
        3,
        min = 0,
        step = 1
      ),
      numericInput("beta", "Effect size per causal marker (beta)", 0.2, step = 0.05), # Input: beta
      numericInput(                                                # Input: residual sigma
        "error_sd",
        "Error SD (sigma)",
        1,
        min = 1e-6,
        step = 0.05
      ),
      numericInput("seed_y", "Seed for y", 1, step = 1),           # Input: RNG seed for y
      selectInput(                                                 # Selector for PVE control mode
        "pve_mode",
        "PVE control mode",
        choices = c(
          "No control" = "none",
          "Fix sigma and scale beta (scale_beta)" = "scale_beta",
          "Fix beta and solve sigma (set_sigma)"   = "set_sigma"
        ),
        selected = "none"
      ),
      sliderInput(                                                 # Slider: target PVE value
        "pve_target",
        "Target PVE (variance explained)",
        min = 0.01,
        max = 0.99,
        value = 0.5,
        step = 0.01
      ),
      actionButton("btn_gen_Y", "Generate y"),                     # Button: trigger y generation
      tags$hr(),                                                   # Divider
      h4("Or upload CSV (clean)"),                                 # Header for upload option
      fileInput("upload", "Upload CSV", accept = ".csv"),          # File input: upload CSV
      helpText("Columns: y, X_1 ... X_m (no missing values).")     # Guidance text for CSV format
    ),
    mainPanel(                                                     # Right panel with outputs
      h4("Preview X (active source)"),                             # Header for X preview
      DTOutput("x_head"),                                          # DataTable showing head of X
      tags$hr(),                                                   # Divider
      h4("Preview data (y + X) — active source"),                  # Header for full data preview
      DTOutput("data_head"),                                       # DataTable with y and X
      tags$hr(),                                                   # Divider
      h4("y histogram (active source)"),                           # Header for histogram
      plotOutput("y_hist", height = 240),                          # Plot: histogram of y
      tags$hr(),                                                   # Divider
      h4("Realized PVE from simulation"),                          # Header for realized PVE text
      uiOutput("pve_text"),                                        # Dynamic UI text about PVE
      tags$hr(),                                                   # Divider
      verbatimTextOutput("data_info")                              # Plain-text summary of dataset
    )
  )),
  
  # ---- model ----
  tabPanel("model", sidebarLayout(                                 # Third tab: modeling; sidebar layout
    sidebarPanel(                                                  # Left panel with modeling controls
      radioButtons(                                                # Choose model kind (SIM vs MMM)
        "model_kind",
        "Model type",
        c(
          "Single-marker models (SIM: y ~ X_j)" = "single",
          "Multiple-marker model (MMM: y ~ all X_*)" = "multi"
        ),
        selected = "single"
      ),
      tags$hr(),                                                   # Divider
      h4("Selection threshold"),                                   # Header for threshold controls
      selectInput(                                                 # Choose statistic for thresholding
        "thr_stat",
        "Statistic",
        choices = c(
          "p-value" = "p",
          "LOD score" = "LOD",
          "LRT statistic" = "LRT"
        ),
        selected = "p"
      ),
      numericInput(                                                # Numeric threshold value
        "thr_value",
        "Threshold value",
        0.05,
        min = 0,
        step = 0.01
      ),
      checkboxInput("thr_bonf", "Bonferroni across markers (for p-value)", FALSE), # Toggle Bonferroni
      helpText(                                                    # Explain how thresholding works
        "For p: selected if p <= threshold (after Bonferroni if enabled). For LOD/LRT: selected if stat >= threshold."
      ),
      tags$hr(),                                                   # Divider
      actionButton("btn_fit", "Fit model")                         # Button: run the selected model(s)
    ),
    mainPanel(                                                     # Right panel with results/plots
      conditionalPanel(                                            # Show this block when SIM is selected
        condition = "input.model_kind == 'single'",                # JS expression toggling visibility
        h4("SIM results (per-marker)"),                            # Header for SIM table
        DTOutput("single_res"),                                    # DataTable of SIM stats
        uiOutput("sel_summary_single"),                            # Text summary of selection counts
        tags$hr(),                                                 # Divider
        h4("SIM visualization"),                                   # Header for SIM plot
        plotOutput("sim_vis", height = 260)                        # Plot: SIM statistics by marker
      ),
      conditionalPanel(                                            # Show this block when MMM is selected
        condition = "input.model_kind == 'multi'",                 # JS expression toggling visibility
        h4("MMM results (per-marker; reduced-vs-full LRT)"),       # Header for MMM table
        DTOutput("multi_res"),                                     # DataTable of MMM stats
        uiOutput("sel_summary_multi"),                             # Text summary for MMM selection
        tags$hr(),                                                 # Divider
        h4("MMM visualization"),                                   # Header for MMM plot
        plotOutput("mmm_vis", height = 260)                        # Plot: MMM statistics by marker
      ),
      tags$hr(),                                                   # Divider below model section
      uiOutput("model_hint"),                                      # Context/help text about scale and stats
      uiOutput("model_status")                                     # Staleness/warning status messages
    )
  )),
  
  # ---- power ----
  tabPanel("power", sidebarLayout(                                 # Fourth tab: power analysis
    sidebarPanel(                                                  # Left panel with power controls
      h4("Choose method"),                                         # Header for method choice
      radioButtons(                                                # Choose between closure and bootstrap methods
        "power_method",
        NULL,
        choices = c(
          "Closure (lambda + bisection)" = "closure",
          "Bootstrap (resampling)" = "bootstrap"
        ),
        selected = "closure"
      ),
      tags$hr(),                                                   # Divider
      conditionalPanel(                                            # Show closure method controls
        condition = "input.power_method == 'closure'",             # Visible only when closure selected
        helpText(                                                  # Short description of the lambda approach
          "NCP/lambda method (single coefficient): lambda_j ≈ n*(beta_j^2/sigma^2)/VIF_j; power via noncentral-t."
        ),
        numericInput(                                              # Given n for power check
          "n_fixed",
          "n to check (given n)",
          1000,
          min = 10,
          step = 10
        ),
        numericInput(                                              # Alpha level for the test
          "alpha_lambda",
          "Significance level (alpha)",
          0.05,
          min = 1e-8,
          step = 0.01
        ),
        numericInput("beta_j", "Assumed beta_j", 0.2, step = 0.01),# Effect size used in lambda
        numericInput(                                              # Sigma used in lambda
          "sigma_j",
          "Assumed sigma (error SD)",
          1,
          min = 1e-8,
          step = 0.01
        ),
        checkboxInput("manual_vif", "Enter VIF manually", FALSE),  # Toggle: manual VIF entry
        conditionalPanel(                                          # Manual VIF input box
          condition = "input.manual_vif == true",
          numericInput(
            "vif_manual",
            "VIF_j (manual)",
            1,
            min = 1,
            step = 0.1
          )
        ),
        conditionalPanel(                                          # Automatic VIF with marker chooser
          condition = "input.manual_vif == false",
          uiOutput("marker_choice_ui"),                            # Dropdown of markers (server-rendered)
          helpText(                                                # Explain how VIF is computed
            "If data is available, VIF_j is computed from regressing X_j on X_-j; otherwise switch to manual VIF."
          )
        ),
        numericInput(                                              # Target power to solve for n*
          "target_power_lambda",
          "Target power (default 0.90)",
          0.90,
          min = 0.5,
          max = 0.999,
          step = 0.01
        ),
        actionButton("btn_lambda", "Compute power / solve n (lambda method)") # Run lambda solver
      ),
      conditionalPanel(                                            # Show bootstrap method controls
        condition = "input.power_method == 'bootstrap'",           # Visible only when bootstrap selected
        helpText(                                                  # Short description of bootstrap scheme
          "Bootstrap resampling: estimate rejection proportion across an n-grid."
        ),
        selectInput(                                               # Choose SIM or MMM for bootstrap
          "boot_model",
          "Model",
          choices = c("single", "multi"),
          selected = "single"
        ),
        numericInput(                                              # Alpha for bootstrap tests
          "alpha_boot",
          "Significance level (alpha)",
          0.05,
          min = 1e-8,
          step = 0.01
        ),
        checkboxInput(                                            # Toggle Bonferroni in bootstrap
          "bonf_boot",
          "Bonferroni across markers (multiple testing)",
          FALSE
        ),
        uiOutput("marker_boot_ui"),                                # Marker selection (server-rendered)
        numericInput(                                              # n-grid lower bound
          "n_min_boot",
          "n grid: min",
          100,
          min = 10,
          step = 10
        ),
        numericInput(                                              # n-grid upper bound
          "n_max_boot",
          "n grid: max",
          1000,
          min = 10,
          step = 10
        ),
        numericInput(                                              # n-grid step size
          "n_step_boot",
          "n grid: step",
          100,
          min = 1,
          step = 1
        ),
        numericInput(                                              # Number of bootstrap repeats
          "B_boot",
          "Bootstrap repeats B",
          1000,
          min = 10,
          step = 10
        ),
        sliderInput(                                               # Target power line for plot/summary
          "target_power_boot",
          "Target power",
          min = 0.5,
          max = 0.99,
          value = 0.9,
          step = 0.01
        ),
        actionButton("btn_boot", "Run bootstrap power")            # Run bootstrap workflow
      ),
      tags$hr(),                                                   # Divider before FDR section
      h4("FDR curve (bootstrap)"),                                 # Header for FDR interface
      # Section title
      helpText(                                                     # Guidance on the FDR estimation approach
        "Estimate FDR(n) with BH q-values via resampling. If 'Use permutation' is checked, we estimate false discoveries from a null permutation (recommended to avoid a flat ~0 curve when signals are strong)."
      ),
      selectInput(                                                 # Choose model for FDR (SIM or MMM)
        "fdr_model",
        "Model",
        choices = c("single", "multi"),
        selected = "multi"
      ),
      # SIM or MMM
      numericInput(                                                # q-value threshold for BH selection
        "fdr_q",
        "q-value threshold (BH)",
        0.05,
        min = 1e-6,
        step = 0.01
      ),
      # BH q threshold
      numericInput("fdr_n_min", "n grid: min", 100, min = 10, step = 10), # n-grid minimum
      # n-grid min
      numericInput(                                                # n-grid maximum
        "fdr_n_max",
        "n grid: max",
        1000,
        min = 10,
        step = 10
      ),
      # n-grid max
      numericInput(                                                # n-grid step size
        "fdr_n_step",
        "n grid: step",
        100,
        min = 1,
        step = 1
      ),
      # n-grid step
      numericInput(                                                # Number of bootstrap repeats for FDR
        "fdr_B",
        "Bootstrap repeats B",
        1000,
        min = 20,
        step = 10
      ),
      # bootstrap B
      numericInput("fdr_seed", "Random seed", 1, min = 1, step = 1), # Seed for reproducibility
      # RNG seed
      checkboxInput(                                               # Toggle permutation-based FDP estimation
        "fdr_use_perm",
        "Use permutation to estimate false discoveries (ignore ground truth)",
        TRUE
      ),
      actionButton("btn_fdr", "Run FDR curve")                     # Run the FDR computation
    ),
    mainPanel(                                                     # Right panel with method outputs
      conditionalPanel(                                            # Show closure-method outputs
        condition = "input.power_method == 'closure'",
        h4("Required sample size (lambda method)"),                # Header for sample size result
        infoBoxOutput("ibox_n", width = 12),                       # Info box with n* and power
        tags$hr(),                                                 # Divider
        h4("NCP/lambda results"),                                  # Header for details text
        uiOutput("lambda_text")                                    # Explanatory text for lambda method
      ),
      conditionalPanel(                                            # Show bootstrap-method outputs
        condition = "input.power_method == 'bootstrap'",
        h4("Bootstrap power results"),                             # Header for power results
        infoBoxOutput("ibox_boot", width = 12),                    # Info box with n* or status
        tags$hr(),                                                 # Divider
        plotOutput("boot_plot", height = 300),                     # Plot of power vs n
        DTOutput("boot_table")                                     # Table of power grid values
      ),
      tags$hr(),                                                   # Divider before FDR outputs
      h4("FDR curve results (BH)"),                                # Header for FDR results
      uiOutput("fdr_status"),                                      # Text status summary for FDR run
      plotOutput("fdr_plot", height = 300),                        # Plot of estimated FDR vs n
      DTOutput("fdr_table")                                        # Table of FDR summary statistics
    )
  )),
  
  # ---- output (report) ----
  tabPanel("output", sidebarLayout(                                # Fifth tab: report export/preview
    sidebarPanel(                                                  # Left panel with report controls
      h4("Export HTML report"),                                    # Header for export section
      textInput("report_title", "Title", "QTL Design Report"),     # Input: report title
      checkboxGroupInput(                                          # Choose which sections to include
        "report_sections",
        "Include sections",
        choices = c("data", "SIM", "MMM", "power", "session"),
        selected = c("data", "SIM", "MMM", "power", "session")
      ),
      numericInput(                                                # Max rows per table in report
        "report_max_rows",
        "Max rows per table",
        50,
        min = 5,
        step = 5
      ),
      actionButton("btn_build_report", "Build report"),            # Button: generate the HTML report
      tags$hr(),                                                   # Divider
      downloadButton("btn_download_html", "Download HTML")         # Button: download the generated file
    ),
    mainPanel(uiOutput("report_preview_hint"), uiOutput("report_preview")) # Right panel: report preview iframe
  ))
)

# -----------------------------
# SERVER
# -----------------------------

server <- function(input, output, session) {                    # Main Shiny server function with I/O bindings
  report_dir <- file.path(tempdir(), "qtl_reports")             # Path to a temp folder where HTML reports will be saved
  dir.create(report_dir, showWarnings = FALSE, recursive = TRUE)# Create the folder if it does not exist
  addResourcePath("qtl_reports", report_dir)                    # Make the folder available to the browser as /qtl_reports
  
  output$session_info <- renderPrint({                          # Output block for the "Session info" area
    sessionInfo()                                               # Print current R session information
  })
  
  # --- reactive state holders ---
  X_gen   <- reactiveVal(NULL)                                  # Stores the generated genotype matrix X
  y_gen   <- reactiveVal(NULL)                                  # Stores the generated phenotype vector y
  df_gen  <- reactiveVal(NULL)                                  # Stores generated data.frame (y + X)
  df_upl  <- reactiveVal(NULL)                                  # Stores uploaded data.frame (y + X)
  stale_model <- reactiveVal(TRUE)                              # Flag: model results are stale and need recompute
  stale_power <- reactiveVal(TRUE)                              # Flag: power/FDR results are stale and need recompute
  rule_of_thumb_n <- function(p)                                # Helper: crude required n for linear models
    ceiling(10 * max(1, p))                                     # Return 10*p (at least 10)
  
  lambda_cache <- reactiveVal(NULL)                             # Cache for lambda (closure) power results
  boot_cache   <- reactiveVal(NULL)   # power bootstrap cache    # Cache for bootstrap power results
  sim_cache    <- reactiveVal(NULL)                              # Cache for SIM statistics table
  mmm_cache    <- reactiveVal(NULL)                              # Cache for MMM statistics table
  fdr_cache    <- reactiveVal(NULL)   # FDR curve cache          # Cache for FDR curve results
  
  # ---- unified bad-data message & checker ----
  bad_data_msg <- "Data do not meet the requirements; please replace the dataset."
  
  # Return TRUE when data are invalid: NA present, no 'y', or no 'X_*' columns
  is_invalid_data <- function(dat) {
    is.null(dat) ||
      !is.data.frame(dat) ||
      anyNA(dat) ||
      !("y" %in% names(dat)) ||
      !any(startsWith(names(dat), "X_"))
  }
  
  output$ibox_n <- renderInfoBox({                              # Default info box for lambda method
    shinydashboard::infoBox(
      "Sample size",                                            # Title shown in the box
      "—",                                                      # Placeholder value
      "Click 'Compute power / solve n' to compute.",            # Helper text
      icon = icon("calculator"),                                # Icon shown on the left
      color = "light-blue",                                     # Color theme
      fill = TRUE                                               # Use filled background style
    )
  })
  output$ibox_boot <- renderInfoBox({                           # Default info box for bootstrap power
    shinydashboard::infoBox(
      "Bootstrap power",                                        # Title shown in the box
      "—",                                                      # Placeholder value
      "Click 'Run bootstrap power' to compute.",                # Helper text
      icon = icon("chart-line"),                                # Icon shown on the left
      color = "light-blue",                                     # Color theme
      fill = TRUE                                               # Use filled background style
    )
  })
  
  use_generated <- function() {                                 # Switch active source to "generated"
    updateRadioButtons(session, "active_source", selected = "generated") # Update UI radio buttons
    stale_model(TRUE)                                           # Mark model outputs as stale
    stale_power(TRUE)                                           # Mark power/FDR outputs as stale
  }
  use_uploaded  <- function() {                                 # Switch active source to "uploaded"
    updateRadioButtons(session, "active_source", selected = "uploaded")  # Update UI radio buttons
    stale_model(TRUE)                                           # Mark model outputs as stale
    stale_power(TRUE)                                           # Mark power/FDR outputs as stale
  }
  observeEvent(input$active_source, {                           # When the user toggles active source
    stale_model(TRUE)                                           # Mark model outputs as stale
    stale_power(TRUE)                                           # Mark power/FDR outputs as stale
  }, ignoreInit = TRUE)                                         # Do not run on initial load
  
  observeEvent(input$upload, {
    file <- input$upload
    req(file)
    
    # Try reading CSV; if failed or invalid, notify and reset uploaded data
    dat <- try(read.csv(file$datapath, check.names = FALSE), silent = TRUE)
    if (inherits(dat, "try-error") || is_invalid_data(dat)) {
      showNotification(bad_data_msg, type = "error", duration = 6)
      df_upl(NULL)
      use_uploaded()
      return()
    }
    
    # If passed all checks, accept the data
    df_upl(dat)
    use_uploaded()
  }, ignoreInit = TRUE)
  
  observeEvent(input$btn_gen_X, {                               # When "Generate X" button is pressed
    G <- sim_X_clean(                                           # Simulate a clean genotype matrix
      n = input$n_ind,                                          # Number of individuals
      m = input$n_mark,                                         # Number of markers
      maf = input$maf,                                          # Minor allele frequency
      seed = input$seed_x                                       # Random seed for reproducibility
    )
    X_gen(G)                                                    # Store new X
    y_gen(NULL)                                                 # Clear any previously simulated y
    df_gen(NULL)                                                # Clear any combined data
    use_generated()                                             # Switch active source to generated; mark stale
  }, ignoreInit = TRUE)                                         # Do not run on initial load
  
  observeEvent(input$btn_gen_Y, {                               # When "Generate y" button is pressed
    req(X_gen())                                                # Require that X has been generated
    pve_target <- if (identical(input$pve_mode, "none"))        # If PVE control is off
      NULL                                                      #   -> no target PVE
    else
      input$pve_target                                          #   -> use user-provided target PVE
    pve_mode   <- if (identical(input$pve_mode, "none"))        # Choose PVE control mode
      "scale_beta"                                              # Default internal mode if disabled
    else
      input$pve_mode                                            # Use selected mode
    y <- sim_y_given_X(                                         # Simulate y from X with given settings
      X_gen(),                                                  # Genotype matrix
      n_causal = input$n_causal,                                # Number of causal markers
      beta = input$beta,                                        # Effect size per causal marker
      error_sd = input$error_sd,                                # Residual standard deviation
      seed = input$seed_y,                                      # Random seed for y
      pve_target = pve_target,                                  # Optional target PVE
      pve_mode = pve_mode                                       # How to enforce the target PVE
    )
    y_gen(y)                                                    # Store generated y
    df_gen(data.frame(y = as.numeric(y), X_gen(),               # Build and store combined (y + X) data.frame
                      check.names = FALSE))
    use_generated()                                             # Switch active source and mark outputs stale
  }, ignoreInit = TRUE)                                         # Do not run on initial load
  
  current_data <- reactive({                                    # Reactive getter for the current full dataset
    if (input$active_source == "uploaded"  &&                   # If uploaded is active and available
        !is.null(df_upl()))
      return(df_upl())                                          # Return uploaded data
    if (input$active_source == "generated" &&                   # If generated is active and available
        !is.null(df_gen()))
      return(df_gen())                                          # Return generated data
    NULL                                                        # Otherwise, no dataset available
  })
  current_X_only <- reactive({                                  # Reactive getter for X-only matrix
    if (input$active_source == "generated" &&                   # If generated X exists
        !is.null(X_gen()))
      return(X_gen())                                           # Return generated X
    if (input$active_source == "uploaded" && !is.null(df_upl())) { # If uploaded data exists
      x_cols <- names(df_upl())[startsWith(names(df_upl()), "X_")] # Identify marker columns
      if (length(x_cols) > 0)
        return(as.matrix(df_upl()[, x_cols, drop = FALSE]))     # Return X as a matrix
    }
    NULL                                                        # Otherwise, no X available
  })
  
  output$x_head <- renderDT({                                   # Render preview table for X
    if (input$active_source == "generated") {                   # For generated source
      if (is.null(X_gen()))                                     # If X not present yet
        return(datatable(                                       # Show guidance message
          data.frame(Info = "No generated X. Click 'Generate X'."),
          options = list(dom = 't')
        ))
      datatable(head(as.data.frame(X_gen()), 10),               # Otherwise show first 10 rows of X
                options = list(scrollX = TRUE, pageLength = 10))
    } else {                                                    # For uploaded source
      if (is.null(df_upl()))                                    # If nothing uploaded
        return(datatable(                                       # Show guidance message
          data.frame(Info = "No uploaded data. Upload a CSV."),
          options = list(dom = 't')
        ))
      x_cols <- names(df_upl())[startsWith(names(df_upl()), "X_")] # Get X_* columns
      if (!length(x_cols))                                      # If no X_* columns exist
        return(datatable(                                       # Show guidance message
          data.frame(Info = "Uploaded data has no X_* columns."),
          options = list(dom = 't')
        ))
      datatable(head(df_upl()[, x_cols, drop = FALSE], 10),     # Show first 10 rows of X columns
                options = list(scrollX = TRUE, pageLength = 10))
    }
  })
  output$data_head <- renderDT({                                # Render preview table for full data (y + X)
    dat <- current_data()                                       # Get active dataset
    if (is.null(dat))                                           # If not available
      return(datatable(                                         # Show guidance message
        data.frame(Info = "No full data (y + X). Generate y or upload CSV."),
        options = list(dom = 't')
      ))
    datatable(head(dat, 10), options = list(scrollX = TRUE, pageLength = 10)) # Show first 10 rows
  })
  output$y_hist <- renderPlot({                                 # Render histogram of y
    dat <- current_data()                                       # Get active dataset
    if (is.null(dat) ||                                         # If dataset missing
        !"y" %in% names(dat)) {                                 # or y column missing
      plot.new()                                                # Draw empty plot
      title("No y yet (generate y or upload CSV including y).") # Add informative title
      return(invisible())                                       # Quit the render function
    }
    hist(dat$y,                                                 # Plot histogram of y
         breaks = 30,
         main = "Histogram of y (active source)",
         xlab = "y")
  })
  output$pve_text <- renderUI({                                 # Render realized/estimated PVE message
    dat <- current_data()                                       # Get active dataset
    if (is.null(dat) ||                                         # If dataset missing
        !"y" %in% names(dat))                                   # or y column missing
      return(tags$em("No y."))                                  # Show placeholder text
    pve_attr <- attr(y_gen(), "pve")                            # Try to read realized PVE from simulated y
    sigma_attr <- attr(y_gen(), "sigma")                        # Try to read sigma used in simulation
    if (!is.null(pve_attr) && is.finite(pve_attr)) {            # If attributes exist and are valid
      HTML(                                                      # Print realized PVE/Sigma
        sprintf(
          "Realized sample PVE approx <b>%.3f</b> (sigma approx %.3f).",
          pve_attr,
          ifelse(is.null(sigma_attr), NA_real_, sigma_attr)
        )
      )
    } else {                                                    # Otherwise estimate a rough PVE using small model
      x_cols <- names(dat)[startsWith(names(dat), "X_")]        # Identify marker columns
      if (!length(x_cols))                                      # If none present
        return(tags$em("No X_* to estimate PVE."))              # Abort with message
      q <- min(3L, length(x_cols))                              # Use first q predictors for quick R^2
      qnames <- paste0("X_", seq_len(q))                        # Build their names
      fit <- try(lm(reformulate(qnames, "y"), data = dat),      # Fit small linear model
                 silent = TRUE)
      if (inherits(fit, "try-error"))                           # If the fit failed
        return(tags$em("Could not estimate PVE."))              # Abort with message
      r2 <- summary(fit)$r.squared                              # Extract R^2 as PVE proxy
      HTML(sprintf(                                             # Show estimated PVE
        "Estimated PVE (via R^2 of y ~ X_1..X_%d) approx <b>%.3f</b>.",
        q,
        r2
      ))
    }
  })
  output$data_info <- renderPrint({                              # Console-style summary of the active dataset
    src <- input$active_source                                  # Read which source is active
    if (src == "generated") {                                   # If using generated data
      if (is.null(X_gen()) &&                                   # If nothing generated yet
          is.null(df_gen())) {
        cat("Active source: Generated - No X or y yet.\n")      # Print message and stop
        return(invisible())
      }
      if (!is.null(df_gen())) {                                 # If full data exists
        dat <- df_gen()                                         # Get it
        p <- sum(startsWith(names(dat), "X_"))                  # Count markers
        cat(                                                    # Print summary lines
          "Active source: Generated\nRows:",
          nrow(dat),
          "  Markers:",
          p,
          "\nColumns:",
          paste(names(dat), collapse = ", "),
          "\n"
        )
      } else {                                                  # If only X exists
        cat(
          "Active source: Generated (X only). Rows:",
          nrow(X_gen()),
          "  Markers:",
          ncol(X_gen()),
          "\n"
        )
      }
    } else {                                                    # If using uploaded data
      if (is.null(df_upl())) {                                  # If nothing uploaded
        cat("Active source: Uploaded - none.\n")                # Print message and stop
        return(invisible())
      }
      dat <- df_upl()                                           # Get uploaded dataset
      p <- sum(startsWith(names(dat), "X_"))                    # Count markers
      cat(                                                      # Print summary lines
        "Active source: Uploaded\nRows:",
        nrow(dat),
        "  Markers:",
        p,
        "\nColumns:",
        paste(names(dat), collapse = ", "),
        "\n"
      )
    }
  })
  
  # threshold auto-suggestion when switching stat                       # Section header: auto-suggest thresholds when the chosen statistic changes
  observeEvent(input$thr_stat, {                                       # React to changes in the 'thr_stat' input selector
    if (identical(input$thr_stat, "p"))                                # If user selects p-value thresholding
      updateNumericInput(session, "thr_value", value = 0.05)           #   Suggest 0.05 as the default threshold
    if (identical(input$thr_stat, "LRT"))                              # If user selects LRT statistic
      updateNumericInput(session, "thr_value", value = 3.84)           #   Suggest 3.84 (≈ chi-square_1,0.05 cutoff)
    if (identical(input$thr_stat, "LOD"))                              # If user selects LOD score
      updateNumericInput(session, "thr_value", value = 0.83)           #   Suggest ~0.83 (conservative; ~3 can be used manually)
  })                                                                   
  
  sim_stats  <- eventReactive(input$btn_fit, {                          # Compute SIM (single-marker) stats when 'Fit model' is clicked
    dat <- current_data()                                               #   Grab the current active dataset
    validate(need(                                                      #   Ensure a dataset exists
      !is.null(dat),
      "No full data. Generate y or switch to an uploaded dataset with y."
    ))
    validate(need("y" %in% names(dat), "Missing column 'y'."))          #   Ensure response 'y' column exists
    validate(need(any(startsWith(                                       #   Ensure there are predictor columns X_*
      names(dat), "X_"
    )), "No X_* columns found."))
    validate(need(!anyNA(dat), bad_data_msg))
    stale_model(FALSE)                                                  #   Mark model outputs as up-to-date
    res <- compute_sim_stats(dat)                                       #   Run per-marker regressions y ~ X_j
    sim_cache(res)                                                      #   Cache results for later use (tables/plots/report)
    res                                                                 #   Return the results to the reactive
  }, ignoreInit = TRUE)                                                 # Do not run until the button is actually clicked
  
  mmm_stats <- eventReactive(input$btn_fit, {                           # Compute MMM (multiple-marker) stats on click
    dat <- current_data()                                               #   Grab the current active dataset
    validate(need(                                                      #   Ensure a dataset exists
      !is.null(dat),
      "No full data. Generate y or switch to uploaded data."
    ))
    validate(need(any(startsWith(                                       #   Ensure there are predictor columns X_*
      names(dat), "X_"
    )), "No X_* columns to fit."))
    validate(need(                                                      #   For full model, require n > p + 1
      nrow(dat) > sum(startsWith(names(dat), "X_")) + 1,
      "MMM requires n > p + 1. Increase n or reduce p."
    ))
    validate(need("y" %in% names(dat), bad_data_msg))
    validate(need(!anyNA(dat), bad_data_msg)) 
    stale_model(FALSE)                                                  #   Mark model outputs as up-to-date
    res <- compute_mmm_stats(dat)                                       #   Fit full model and reduced models to get LRT/LOD per marker
    mmm_cache(res)                                                      #   Cache results for later use
    res                                                                 #   Return the results
  }, ignoreInit = TRUE)                                                 # Do not run until the button is clicked
  
  apply_threshold <- function(df,                                       # Helper: flag selected markers given a thresholding rule
                              stat = c("p", "LOD", "LRT"),              #   Which statistic to threshold on
                              thr = 0.05,                               #   Threshold value
                              bonf = FALSE) {                           #   Whether to use Bonferroni for p-values
    stat <- match.arg(stat)                                             #   Resolve the chosen statistic
    if (nrow(df) == 0)                                                  #   If no rows, return as-is
      return(df)
    p <- nrow(df)                                                       #   Number of tests/markers (used for Bonferroni)
    if (stat == "p") {                                                  #   p-value case
      a <- if (bonf)                                                    #     Effective alpha: Bonferroni if requested
        thr / max(1L, p)
      else
        thr
      df$selected <- is.finite(df$p) & (df$p <= a)                      #     Select markers with p <= alpha_eff
      attr(df, "alpha_eff") <- a                                        #     Store effective alpha for display
    } else if (stat == "LOD") {                                         #   LOD case
      df$selected <- is.finite(df$LOD) & (df$LOD >= thr)                #     Select markers with LOD >= threshold
    } else {                                                            #   LRT case
      df$selected <- is.finite(df$LRT) & (df$LRT >= thr)                #     Select markers with LRT >= threshold
    }
    df                                                                   #   Return augmented data frame
  }
  
  output$single_res <- renderDT({                                       # Render SIM results table
    df <- sim_stats()                                                   #   Obtain SIM results reactive
    req(df)                                                             #   Require data to proceed
    df2 <- df                                                           
    if ("p" %in% names(df2))                                            #   Format p-values to 4 significant digits
      df2$p        <- signif(df2$p, 4)
    if ("beta_hat" %in% names(df2))                                     #   Format coefficient estimates
      df2$beta_hat <- signif(df2$beta_hat, 4)
    if ("t" %in% names(df2))                                            #   Round t statistics
      df2$t        <- round(df2$t, 3)
    if ("LRT" %in% names(df2))                                          #   Round LRT values
      df2$LRT      <- round(df2$LRT, 3)
    if ("LOD" %in% names(df2))                                          #   Round LOD values
      df2$LOD      <- round(df2$LOD, 3)
    df2 <- apply_threshold(                                             #   Apply user-selected thresholding rule
      df2,
      stat = input$thr_stat,
      thr = input$thr_value,
      bonf = isTRUE(input$thr_bonf)
    )
    df2 <- cbind(marker = rownames(df2), df2)                           #   Add marker name as a column
    datatable(df2, options = list(pageLength = 10))                     #   Show as a DataTable with 10 rows per page
  })
  output$sel_summary_single <- renderUI({                               # Render a summary text for SIM selections
    df <- sim_stats()                                                   #   Obtain SIM results
    req(df)                                                             #   Require data
    df2 <- apply_threshold(                                             #   Apply thresholding to compute selections
      df,
      stat = input$thr_stat,
      thr = input$thr_value,
      bonf = isTRUE(input$thr_bonf)
    )
    k  <- sum(df2$selected, na.rm = TRUE)  # FIXED                       #   Count selected markers safely
    a_eff <- attr(df2, "alpha_eff")                                      #   Retrieve effective alpha if present
    if (input$thr_stat == "p" && is.finite(a_eff)) {                     #   If p-based rule, include alpha_eff in the message
      HTML(
        sprintf(
          "Selected markers (SIM) = <b>%d</b> / %d (effective alpha = %.3g).",
          k,
          nrow(df2),
          a_eff
        )
      )
    } else {                                                             #   Otherwise, report count without alpha_eff
      HTML(sprintf("Selected markers (SIM) = <b>%d</b> / %d.", k, nrow(df2)))
    }
  })
  
  output$multi_res <- renderDT({                                        # Render MMM results table
    df <- mmm_stats()                                                   #   Obtain MMM results reactive
    req(df)                                                             #   Require data
    df2 <- df
    if ("p" %in% names(df2))                                            #   Format p-values to 4 significant digits
      df2$p        <- signif(df2$p, 4)
    if ("beta_hat" %in% names(df2))                                     #   Format coefficient estimates
      df2$beta_hat <- signif(df2$beta_hat, 4)
    if ("t" %in% names(df2))                                            #   Round t statistics
      df2$t        <- round(df2$t, 3)
    if ("LRT" %in% names(df2))                                          #   Round LRT values
      df2$LRT      <- round(df2$LRT, 3)
    if ("LOD" %in% names(df2))                                          #   Round LOD values
      df2$LOD      <- round(df2$LOD, 3)
    df2 <- apply_threshold(                                             #   Apply user-selected thresholding rule
      df2,
      stat = input$thr_stat,
      thr = input$thr_value,
      bonf = isTRUE(input$thr_bonf)
    )
    df2 <- cbind(marker = rownames(df2), df2)                           #   Add marker name as a column
    datatable(df2, options = list(pageLength = 10))                     #   Show as a DataTable with 10 rows per page
  })
  output$sel_summary_multi <- renderUI({                                # Render a summary text for MMM selections
    df <- mmm_stats()                                                   #   Obtain MMM results
    req(df)                                                             #   Require data
    df2 <- apply_threshold(                                             #   Apply thresholding to compute selections
      df,
      stat = input$thr_stat,
      thr = input$thr_value,
      bonf = isTRUE(input$thr_bonf)
    )
    k  <- sum(df2$selected, na.rm = TRUE)  # FIXED                       #   Count selected markers safely
    a_eff <- attr(df2, "alpha_eff")                                      #   Retrieve effective alpha if present
    if (input$thr_stat == "p" && is.finite(a_eff)) {                     #   If p-based rule, include alpha_eff in the message
      HTML(
        sprintf(
          "Selected markers (MMM) = <b>%d</b> / %d (effective alpha = %.3g).",
          k,
          nrow(df2),
          a_eff
        )
      )
    } else {                                                             #   Otherwise, report count without alpha_eff
      HTML(sprintf("Selected markers (MMM) = <b>%d</b> / %d.", k, nrow(df2)))
    }
  })
  
  plot_stat <- function(df,                                             # Helper: generic bar plot of a chosen statistic with threshold line
                        stat,
                        thr,
                        bonf = FALSE,
                        title_prefix = "") {
    df$marker <- factor(rownames(df), levels = rownames(df))            #   Keep marker order consistent on the x-axis
    if (stat == "p") {                                                  #   p-value case
      a_eff <- if (bonf)                                                #     Effective alpha with/without Bonferroni
        thr / max(1L, nrow(df))
      else
        thr
      df$yval <- -log10(df$p)                                           #     Plot -log10(p) so higher bars mean stronger evidence
      ylab <- "-log10(p)"                                               #     y-axis label
      thr_line <- -log10(a_eff)                                         #     Threshold line in -log10 scale
      sub <- sprintf("Dashed line at -log10(alpha%s)=%.3f",             #     Subtitle describing the threshold
                     if (bonf)
                       " (Bonferroni)"
                     else
                       "",
                     thr_line)
    } else if (stat == "LOD") {                                         #   LOD case
      df$yval <- df$LOD                                                 #     Use LOD as the bar height
      ylab <- "LOD"                                                     #     y-axis label
      thr_line <- thr                                                   #     Threshold line at user value
      sub <- sprintf("Dashed line at LOD = %.3f", thr)                  #     Subtitle describing the threshold
    } else {                                                            #   LRT case
      df$yval <- df$LRT                                                 #     Use LRT as the bar height
      ylab <- "LRT"                                                     #     y-axis label
      thr_line <- thr                                                   #     Threshold line at user value
      sub <- sprintf("Dashed line at LRT = %.3f", thr)                  #     Subtitle describing the threshold
    }
    df <- df[is.finite(df$yval), , drop = FALSE]                        #   Keep only finite values to plot
    if (!nrow(df)) {                                                    #   If nothing to plot, draw an empty panel with a note
      plot.new()
      title("No finite values to plot")
      return(invisible())
    }
    ggplot(df, aes(x = marker, y = yval)) +                             #   Create bar chart by marker
      geom_col() +                                                      #   Draw columns
      geom_hline(yintercept = thr_line, linetype = "dashed") +          #   Add dashed threshold line
      labs(                                                             #   Axes labels and title
        x = "Marker",
        y = ylab,
        title = sprintf("%s %s", title_prefix, switch(                  #   Title includes prefix and stat name
          stat,
          p = "p-value",
          LOD = "LOD",
          LRT = "LRT"
        )),
        subtitle = sub                                                  #   Subtitle shows threshold info
      ) +
      theme_minimal() +                                                 #   Minimal theme for clean look
      theme(axis.text.x = element_text(                                 #   Rotate x labels for readability
        angle = 90,
        vjust = 0.5,
        hjust = 1
      ))
  }
  output$sim_vis <- renderPlot({                                        # Render SIM plot panel
    df <- sim_stats()                                                   #   Obtain SIM results
    req(df)                                                             #   Require data
    plot_stat(                                                          #   Produce the plot with user-chosen statistic and threshold
      df,
      stat = input$thr_stat,
      thr = input$thr_value,
      bonf = isTRUE(input$thr_bonf),
      title_prefix = "SIM ·"
    )
  })
  output$mmm_vis <- renderPlot({                                        # Render MMM plot panel
    df <- mmm_stats()                                                   #   Obtain MMM results
    req(df)                                                             #   Require data
    plot_stat(                                                          #   Produce the plot with user-chosen statistic and threshold
      df,
      stat = input$thr_stat,
      thr = input$thr_value,
      bonf = isTRUE(input$thr_bonf),
      title_prefix = "MMM ·"
    )
  })
  
  output$model_hint <- renderUI({                                       # Render helper text with scale checks and tips
    if (is.null(current_data()))                                        #   If no data yet, show a simple note
      return(tagList(h4("Scale check / intuition"), p("No data yet.")))
    dat <- current_data()                                               #   Get active dataset
    p <- sum(startsWith(names(dat), "X_"))                              #   Count predictors
    n <- nrow(dat)                                                      #   Count rows (sample size)
    tagList(                                                             #   Build a small info block
      h4("Scale check / intuition"),
      p(sprintf("Current: n = %d, p = %d.", n, p)),                     #   Display current n and p
      p(                                                                
        sprintf(
          "Rule-of-thumb: n >= 10 * p approx %d (for linear regression).",
          rule_of_thumb_n(p)                                            #   Show a rough rule-of-thumb for n vs p
        )
      ),
      p(                                                                
        "SIM: per-marker null (y ~ 1) vs full (y ~ X_j): report t, p, LRT = n * log(RSS0/RSS1), LOD = (n/2) * log10(RSS0/RSS1)." #   Explain SIM outputs
      ),
      p(                                                                
        "MMM: full (y ~ all X) vs reduced (drop X_j): report coef t/p plus LRT_j, LOD_j." #   Explain MMM outputs
      )
    )
  })
  output$model_status <- renderUI({                                     # Render a small warning if model results are stale
    if (isTRUE(stale_model()))                                          #   If model outputs are out-of-date
      tags$em("Model outputs are outdated due to data changes. Click 'Fit model' to refresh.") #   Ask user to refit
    else                                                                #   Otherwise, show nothing
      NULL
  })
  
  # ---- closure (lambda) helpers & observers ----                      # Section header: analytical power / sample size helpers
  compute_vif_j <- function(X, j) {                                     # Compute VIF for predictor j in matrix X
    p <- ncol(X)                                                        #   Number of predictors
    if (p <= 1)                                                         #   If only one predictor, VIF is trivially 1
      return(1)
    others <- setdiff(seq_len(p), j)                                    #   Indices of all other predictors
    if (nrow(X) <= length(others) + 1)                                  #   Need enough rows to fit regression
      return(NA_real_)
    df_tmp <- data.frame(xj = X[, j], X[, others, drop = FALSE])        #   Data frame with xj and the others
    fit <- try(lm(xj ~ ., data = df_tmp), silent = TRUE)                #   Regress xj on remaining predictors
    if (inherits(fit, "try-error"))                                     #   If fit fails, return NA
      return(NA_real_)
    r2 <- max(0, min(1, summary(fit)$r.squared))                        #   Clamp R^2 into [0,1]
    1 / max(1e-12, 1 - r2)                                              #   VIF = 1 / (1 - R^2), protect against division by ~0
  }
  lambda_power_for_n <- function(n, beta, sigma, vif, p, alpha) {       # Approximate power for given n using NCP (noncentral t)
    n <- as.integer(n)                                                  #   Ensure integer n
    p <- as.integer(p)                                                  #   Ensure integer p
    if (n <= p + 1)                                                     #   Need n > p + 1 for valid t-test in multiple regression
      return(0)
    nu <- n - p - 1                                                     #   Degrees of freedom for t-statistic
    lambda <- n * (beta^2 / (sigma^2)) * (1 / vif)                      #   Noncentrality parameter λ ≈ n*(β^2/σ^2)/VIF
    delta  <- sqrt(lambda)                                              #   Noncentral t uses ncp = √λ
    tcrit  <- qt(1 - alpha / 2, df = nu)                                #   Two-sided t critical value at level α
    pow <- pt(-tcrit, df = nu, ncp = delta) +                           #   Power = P(T < -tcrit) + P(T > tcrit)
      (1 - pt(tcrit, df = nu, ncp = delta))
    as.numeric(max(0, min(1, pow)))                                     #   Clamp result to [0,1]
  }
  solve_n_lambda_halving_then_bisect <- function(target_power,          # Find minimal n meeting target power by bracketing + bisection
                                                 beta,
                                                 sigma,
                                                 vif,
                                                 p,
                                                 alpha,
                                                 n_init = 1000,
                                                 tol = 1,
                                                 max_iter = 60) {
    hi <- max(as.integer(n_init), p + 2)                                #   Start with upper bound ≥ p+2
    pow_hi <- lambda_power_for_n(hi, beta, sigma, vif, p, alpha)        #   Compute power at upper bound
    grow_guard <- 0                                                     #   Safety counter for growth loop
    while (pow_hi < target_power && hi < 1e6 && grow_guard < 30) {      #   Increase upper bound until power is reached
      hi <- as.integer(ceiling(hi * 1.5))                               #     Geometric growth of n
      pow_hi <- lambda_power_for_n(hi, beta, sigma, vif, p, alpha)      #     Recompute power at new hi
      grow_guard <- grow_guard + 1                                      #     Increment safety counter
    }
    if (pow_hi < target_power)                                          #   If even huge n fails, return NA bracket
      return(
        list(
          n_star = NA_integer_,
          lo = NA_integer_,
          hi = NA_integer_,
          pow_lo = NA_real_,
          pow_hi = NA_real_
        )
      )
    lo <- max(p + 2, floor(hi / 2))                                     #   Initialize lower bound
    pow_lo <- lambda_power_for_n(lo, beta, sigma, vif, p, alpha)        #   Power at lower bound
    shrink_guard <- 0                                                   #   Safety counter for shrink loop
    while (pow_lo >= target_power &&                                    #   If low bound already meets power,
           lo > (p + 2) && shrink_guard < 30) {                         #   try shrinking the bracket
      hi <- lo                                                          #     Move hi down
      pow_hi <- pow_lo                                                  #     Update power at hi
      lo <- max(p + 2, floor(lo / 2))                                   #     Shrink lo further
      pow_lo <- lambda_power_for_n(lo, beta, sigma, vif, p, alpha)      #     Recompute power at new lo
      shrink_guard <- shrink_guard + 1                                  #     Increment safety counter
    }
    if (!(pow_lo < target_power && pow_hi >= target_power))             #   Ensure bracket actually straddles the target
      return(list(
        n_star = hi,
        lo = lo,
        hi = hi,
        pow_lo = pow_lo,
        pow_hi = pow_hi
      ))
    it <- 0                                                             #   Bisection iteration counter
    while ((hi - lo) > tol && it < max_iter) {                          #   Standard integer bisection until tolerance
      it  <- it + 1                                                     #     Increment iteration
      mid <- as.integer(floor((lo + hi) / 2))                           #     Midpoint of bracket
      pow_mid <- lambda_power_for_n(mid, beta, sigma, vif, p, alpha)    #     Power at midpoint
      if (pow_mid >= target_power) {                                    #     If power at mid is enough, move hi down
        hi <- mid
        pow_hi <- pow_mid
      } else {                                                          #     Else move lo up
        lo <- mid
        pow_lo <- pow_mid
      }
    }
    list(                                                                #   Return minimal n and diagnostic values
      n_star = hi,
      lo = lo,
      hi = hi,
      pow_lo = pow_lo,
      pow_hi = pow_hi
    )
  }
  
  output$marker_choice_ui <- renderUI({                                                # Render the UI element used to choose a marker j for VIF/power calc
    Xonly <- current_X_only()                                                          # Get the current X-only matrix (predictors) from reactive state
    if (is.null(Xonly))                                                                # If we don't have X available
      return(tags$em("No X available; use manual VIF."))                               #   Inform the user to switch to manual VIF entry
    xnames <- colnames(Xonly)                                                          # Otherwise, obtain column (marker) names
    selectInput("marker_j",                                                            # Create a dropdown for picking a marker j
                "Choose marker j",                                                     #   Label for the dropdown
                choices = xnames,                                                      #   Available choices are the marker names
                selected = xnames[1])                                                  #   Default selection: the first marker
  })
  
  observeEvent(input$btn_lambda, {                                                     # When the "Compute power / solve n" (lambda method) button is clicked
    isolate({                                                                          # Use isolate so downstream reactives are not triggered unnecessarily
      Xonly <- current_X_only()                                                        # Fetch the current X-only design matrix
      if (is.null(Xonly) && !isTRUE(input$manual_vif)) {                               # If no X and manual VIF is not enabled
        output$ibox_n <- renderInfoBox({                                               #   Show an info box explaining the issue
          shinydashboard::infoBox(
            "Sample size (lambda)",                                                    #   Title of the info box
            "—",                                                                       #   Value placeholder (dash)
            "No X available to compute VIF. Use manual VIF.",                          #   Message instructing to use manual VIF
            icon = icon("triangle-exclamation"),                                       #   Warning icon
            color = "red",                                                             #   Red color to signal an issue
            fill = TRUE                                                                #   Fill the box background
          )
        })
        return()                                                                       #   Abort the handler early
      }
      if (!is.null(Xonly)) {                                                           # If X is available
        p <- ncol(Xonly)                                                               #   Count the number of predictors (p) for df and formulas
      } else {                                                                         # Else (no X and manual VIF must be used)
        output$ibox_n <- renderInfoBox({                                               #   Still cannot proceed because we also need p
          shinydashboard::infoBox(
            "Sample size (lambda)",                                                    #   Title
            "—",                                                                       #   Value placeholder
            "Please provide X to determine p (number of markers).",                    #   Ask user to provide X so we can know p
            icon = icon("triangle-exclamation"),                                       #   Warning icon
            color = "red",                                                             #   Red color
            fill = TRUE                                                                #   Fill background
          )
        })
        return()                                                                       #   Abort the handler
      }
      if (isTRUE(input$manual_vif)) {                                                  # If user selected manual VIF mode
        vif_j <- max(1, as.numeric(input$vif_manual))                                  #   Read VIF value (coerce numeric) and clamp to at least 1
        marker_lab <- "(manual VIF)"                                                   #   Label to indicate manual VIF was used
      } else {                                                                         # Otherwise, compute VIF from the data
        j <- match(input$marker_j, colnames(Xonly))                                    #   Find the column index of the chosen marker
        vif_j <- compute_vif_j(Xonly, j)                                               #   Compute VIF for that marker
        marker_lab <- if (is.finite(vif_j))                                            #   Build a label for display depending on success
          paste0(input$marker_j, sprintf(" (VIF = %.3f)", vif_j))                      #     If finite, show the numeric VIF
        else                                                                           #     Otherwise
          "(VIF unavailable)"                                                          #     Mark that VIF could not be computed
      }
      if (!is.finite(vif_j)) {                                                         # If we failed to compute a finite VIF
        output$ibox_n <- renderInfoBox({                                               #   Show an error info box suggesting manual VIF
          shinydashboard::infoBox(
            "Sample size (lambda)",                                                    #   Title
            "—",                                                                       #   Value placeholder
            "Could not compute VIF; try manual VIF.",                                  #   Error message
            icon = icon("triangle-exclamation"),                                       #   Warning icon
            color = "red",                                                             #   Red color
            fill = TRUE                                                                #   Fill background
          )
        })
        return()                                                                       #   Abort the handler
      }
      n_given <- as.integer(input$n_fixed)                                             # Read the user-provided n at which to compute power
      alpha <- input$alpha_lambda                                                      # Read significance level alpha
      beta  <- input$beta_j                                                            # Read assumed effect size beta_j
      sigma <- input$sigma_j                                                           # Read assumed error SD (sigma)
      target <- input$target_power_lambda                                              # Read target power level
      pow_n <- lambda_power_for_n(n_given, beta, sigma, vif_j, p, alpha)               # Compute power at the given n using lambda/NCP formula
      res   <- solve_n_lambda_halving_then_bisect(                                     # Solve for minimal n* that achieves target power
        target,                                                                        #   Target power
        beta,                                                                          #   Effect size
        sigma,                                                                         #   Error SD
        vif_j,                                                                         #   VIF for marker j
        p,                                                                             #   Number of predictors
        alpha,                                                                         #   Significance level
        n_init = n_given,                                                              #   Start search near the provided n
        tol = 1,                                                                       #   Tolerance on n (integer step)
        max_iter = 60                                                                  #   Max iterations for bisection
      )
      n_min <- res$n_star                                                               # Extract the solution n* (may be NA if unattainable)
      lo_fin <- res$lo                                                                  # Lower bracket end at termination
      hi_fin <- res$hi                                                                  # Upper bracket end at termination
      pow_lo <- res$pow_lo                                                              # Power at lower bracket
      pow_hi <- res$pow_hi                                                              # Power at upper bracket
      if (is.na(n_min)) {                                                               # If the target is unattainable under current settings
        output$ibox_n <- renderInfoBox({                                                #   Show a red warning box
          shinydashboard::infoBox(
            "Sample size (lambda)",                                                     #   Title
            "—",                                                                        #   Value placeholder
            "Target power may be unattainable under current settings.",                 #   Message
            icon = icon("triangle-exclamation"),                                        #   Warning icon
            color = "red",                                                              #   Red color
            fill = TRUE                                                                 #   Fill background
          )
        })
      } else {                                                                          # If a feasible n* was found
        meets <- if (n_given >= n_min)                                                  #   Determine if the provided n meets target power
          "meets target"                                                                #     Yes
        else                                                                            #     No
          "below target"                                                                #     Below target
        output$ibox_n <- renderInfoBox({                                                #   Display a blue info box with n* and power at n_given
          shinydashboard::infoBox(
            "Sample size (lambda)",                                                     #   Title
            paste0("n* approx ", n_min),                                                #   Show the computed n*
            sprintf(                                                                     #   Subtitle: power at given n and status
              "Given n = %d -> power approx %.3f (%s)",
              n_given,
              pow_n,
              meets
            ),
            icon = icon("calculator"),                                                  #   Calculator icon
            color = "blue",                                                             #   Blue color (informational)
            fill = TRUE                                                                 #   Fill background
          )
        })
      }
      stale_power(FALSE)                                                                # Mark power outputs as up-to-date
      lambda_cache(                                                                     # Cache all relevant values for reports/preview
        list(
          marker_label = marker_lab,                                                    #   Text label of marker/VIF used
          p = p,                                                                        #   Number of predictors
          alpha = alpha,                                                                #   Significance level
          beta = beta,                                                                  #   Effect size
          sigma = sigma,                                                                #   Error SD
          vif = vif_j,                                                                  #   VIF value
          n_given = n_given,                                                            #   Provided n
          power_given = pow_n,                                                          #   Power at provided n
          target = target,                                                              #   Target power
          n_star = n_min,                                                               #   Minimal n achieving target
          bracket = c(lo_fin, hi_fin),                                                  #   Final bracket [lo, hi]
          pow_lo = pow_lo,                                                              #   Power at lo
          pow_hi = pow_hi                                                               #   Power at hi
        )
      )
    })
  })
  output$lambda_text <- renderUI({                                                      # Render a textual explanation of the lambda/NCP calculation
    input$btn_lambda                                                                    # Depend on the button click (for re-rendering)
    isolate({                                                                           # Isolate to avoid unnecessary reactive churn
      Xonly <- current_X_only()                                                         # Get current X-only matrix
      if (is.null(Xonly) && !isTRUE(input$manual_vif))                                  # If no X and manual VIF not enabled
        return(tags$em("No X available to compute VIF. Check 'Enter VIF manually'."))   #   Ask user to switch to manual VIF
      if (!is.null(Xonly))                                                              # If X is available
        p <- ncol(Xonly)                                                                #   Set p to number of columns
      else                                                                              # Otherwise we still cannot proceed
        return(tags$em("Please provide X to determine p (number of markers)."))         #   Ask user to provide X
      if (isTRUE(input$manual_vif)) {                                                   # If manual VIF mode is on
        vif_j <- max(1, as.numeric(input$vif_manual))                                   #   Read and clamp manual VIF
        marker_lab <- "(manual VIF)"                                                    #   Label to indicate manual VIF
      }
      else {                                                                            # Otherwise, compute VIF from data
        j <- match(input$marker_j, colnames(Xonly))                                     #   Resolve marker index
        if (is.na(j))                                                                   #   If not found or invalid
          return(tags$em("Invalid marker selection."))                                  #     Inform the user
        vif_j <- compute_vif_j(Xonly, j)                                                #   Compute VIF
        if (!is.finite(vif_j))                                                          #   If VIF computation failed
          return(tags$em("Could not compute VIF; use manual VIF."))                     #     Suggest manual VIF
        marker_lab <- paste0(input$marker_j, sprintf(" (VIF = %.3f)", vif_j))           #   Label with the computed VIF value
      }
      n_given <- as.integer(input$n_fixed)                                              # Read provided n
      alpha <- input$alpha_lambda                                                       # Read alpha
      beta  <- input$beta_j                                                             # Read beta_j
      sigma <- input$sigma_j                                                            # Read sigma
      target <- input$target_power_lambda                                               # Read target power
      pow_n <- lambda_power_for_n(n_given, beta, sigma, vif_j, p, alpha)                # Compute power at given n
      res   <- solve_n_lambda_halving_then_bisect(                                      # Solve for minimal n*
        target,                                                                         #   Target power
        beta,                                                                           #   Effect size
        sigma,                                                                          #   Error SD
        vif_j,                                                                          #   VIF
        p,                                                                              #   Number of predictors
        alpha,                                                                          #   Significance level
        n_init = n_given,                                                               #   Start search near provided n
        tol = 1,                                                                        #   Tolerance on n
        max_iter = 60                                                                   #   Max iterations
      )
      n_min <- res$n_star                                                                # Extract n*
      if (is.null(n_min) || is.na(n_min)) {                                              # If unattainable
        tagList(                                                                         #   Build explanatory UI
          p(HTML(                                                                         #   Paragraph with parameter summary
            sprintf(
              "<b>Marker:</b> %s; <b>p</b> = %d; <b>alpha</b> = %.3g; <b>beta_j</b> = %.3f; <b>sigma</b> = %.3f.",
              marker_lab,
              p,
              alpha,
              beta,
              sigma
            )
          )),
          p(                                                                              #   Paragraph explaining the lambda formula
            HTML(
              "Using formula: lambda_j approx n * (beta_j^2/sigma^2) * 1/VIF_j."
            )
          ),
          tags$em(                                                                        #   Note that target may be unattainable
            "Target power may be unattainable under current settings."
          )
        )
      } else {                                                                            # If feasible
        tagList(                                                                          #   Build explanatory UI
          p(HTML(                                                                          #   Show the settings summary
            sprintf(
              "<b>Marker:</b> %s; <b>p</b> = %d; <b>alpha</b> = %.3g; <b>beta_j</b> = %.3f; <b>sigma</b> = %.3f.",
              marker_lab,
              p,
              alpha,
              beta,
              sigma
            )
          )),
          p(HTML(                                                                          #   Reiterate the lambda approximation formula
            "Using formula: lambda_j approx n * (beta_j^2/sigma^2) * 1/VIF_j."
          )),
          h4(                                                                              #   Display the computed minimal n*
            sprintf(
              "Minimal n to reach target power (%.2f): n* approx %d",
              target,
              n_min
            )
          ),
          p(                                                                               #   Show the final bracket and power at the bounds
            sprintf(
              "Final bracket: [%d, %d] with power(lo) = %.3f, power(hi) = %.3f",
              res$lo,
              res$hi,
              res$pow_lo,
              res$pow_hi
            )
          )
        )
      }
    })
  })
  
  # ---- bootstrap power ----
  output$marker_boot_ui <- renderUI({                                                     # Render UI for choosing a marker used in bootstrap power
    dat <- current_data()                                                                  # Get the currently active dataset (y + X_*)
    if (is.null(dat))                                                                      # If there is no data available yet
      return(tags$em("No data yet (need y + X_*)"))                                        #   Show a hint that both y and X_* are required
    x_cols <- names(dat)[startsWith(names(dat), "X_")]                                     # Extract all predictor column names starting with "X_"
    if (!length(x_cols))                                                                   # If there are no X_* columns
      return(tags$em("No X_* columns in current data."))                                   #   Inform the user
    selectInput("marker_boot",                                                             # Create a dropdown to select the marker tested in bootstrap
                "Marker to test",                                                          #   Label for the dropdown
                choices = x_cols,                                                          #   Choices are all X_* columns
                selected = x_cols[1])                                                      #   Default to the first marker
  })
  
  bootstrap_power <- function(dat,                                                         # Function to compute bootstrap-based power across an n grid
                              model = c("single", "multi"),                                #   Model type: single-marker (SIM) or multi-marker (MMM)
                              alpha = 0.05,                                                #   Nominal significance level
                              n_grid,                                                      #   Vector of sample sizes to evaluate
                              B = 1000,                                                     #   Number of bootstrap resamples per n
                              marker,                                                      #   Marker to test significance for
                              bonferroni = FALSE) {                                        #   Whether to use Bonferroni correction across p markers
    model <- match.arg(model)                                                              # Match and validate the model argument
    x_cols <- names(dat)[startsWith(names(dat), "X_")]                                     # Get X_* columns
    if (length(x_cols) == 0)                                                               # If none found
      stop("No X_* columns.", call. = FALSE)                                               #   Stop with a clear error
    if (!("y" %in% names(dat)))                                                            # Ensure response y is present
      stop("Missing column 'y'.", call. = FALSE)                                           #   Stop if missing
    if (is.null(marker) ||                                                                 # Validate the marker argument
        !(marker %in% x_cols))                                                             #   Must be one of X_* columns
      stop("Invalid marker selection.", call. = FALSE)                                     #   Stop if invalid
    p <- length(x_cols)                                                                    # Number of predictors
    alpha_eff <- if (bonferroni)                                                           # Effective alpha after multiple-testing correction (if any)
      alpha / max(1L, p)                                                                   #   Bonferroni: alpha / p
    else                                                                                   # Otherwise
      alpha                                                                                #   Use nominal alpha
    pow <- numeric(length(n_grid))                                                         # Initialize power vector across the n grid
    for (i in seq_along(n_grid)) {                                                         # Loop over each n in the grid
      n_i <- n_grid[i]                                                                     #   Current sample size to evaluate
      rej <- 0L                                                                            #   Counter of rejections across B resamples
      for (b in seq_len(B)) {                                                              #   Bootstrap loop
        idx <- sample.int(nrow(dat), size = n_i, replace = TRUE)                           #     Sample n_i rows with replacement
        d2  <- dat[idx, , drop = FALSE]                                                    #     Create the bootstrap sample
        if (model == "single") {                                                           #     If single-marker model
          fit <- try(lm(reformulate(marker, "y"), data = d2), silent = TRUE)               #       Fit y ~ marker
          if (!inherits(fit, "try-error")) {                                               #       If fit succeeds
            co <- summary(fit)$coefficients                                                #         Extract coefficient table
            pval <- if (nrow(co) >= 2)                                                     #         Get p-value for the marker term
              co[2, 4]                                                                      #           p-value in the second row (slope)
            else                                                                            #         If the term is missing
              NA_real_                                                                      #           Mark as NA
            if (is.finite(pval) && pval < alpha_eff)                                       #         If p-value is valid and below threshold
              rej <- rej + 1L                                                               #           Count a rejection
          }
        } else {                                                                           #     Else multi-marker model
          fit <- try(lm(reformulate(x_cols, "y"), data = d2), silent = TRUE)               #       Fit y ~ all X_* (full model)
          if (!inherits(fit, "try-error")) {                                               #       If fit succeeds
            sm <- summary(fit)$coefficients                                                #         Extract coefficient summary
            if (marker %in% rownames(sm)) {                                                #         Ensure marker row exists
              pval <- sm[marker, 4]                                                        #           Get p-value of the marker
              if (is.finite(pval) &&                                                       #           If p-value is valid and
                  pval < alpha_eff)                                                        #           below threshold
                rej <- rej + 1L                                                            #             Count a rejection
            }
          }
        }
      }
      pow[i] <- rej / B                                                                     #   Power at n_i = proportion of rejections
    }
    data.frame(n = n_grid, power = pow)                                                     # Return power curve as a data frame
  }
  
  observeEvent(input$btn_boot, {                                                            # Handler for "Run bootstrap power" button
    dat <- current_data()
    if (is.null(dat) || !("y" %in% names(dat)) ||
        !any(startsWith(names(dat), "X_")) || anyNA(dat)) {
      output$ibox_boot <- renderInfoBox({
        shinydashboard::infoBox(
          "Bootstrap power", "—",
          bad_data_msg,
          icon = icon("triangle-exclamation"),
          color = "red", fill = TRUE
        )
      })
      return()
    }
    x_cols <- names(dat)[startsWith(names(dat), "X_")]                                      # Extract X_* columns
    if (length(x_cols) == 0) {                                                              # If none present
      output$ibox_boot <- renderInfoBox({                                                   #   Display an error info box
        shinydashboard::infoBox(
          "Bootstrap power",                                                                #   Title
          "—",                                                                              #   Placeholder
          "No X_* columns in current data.",                                                #   Message
          icon = icon("triangle-exclamation"),                                              #   Warning icon
          color = "red",                                                                     #   Red color
          fill = TRUE                                                                        #   Filled style
        )
      })
      return()                                                                              #   Abort
    }
    if (!("y" %in% names(dat))) {                                                           # Ensure y exists
      output$ibox_boot <- renderInfoBox({                                                   #   If missing, show error info box
        shinydashboard::infoBox(
          "Bootstrap power",                                                                #   Title
          "—",                                                                              #   Placeholder
          "Current data has no 'y' column.",                                                #   Message
          icon = icon("triangle-exclamation"),                                              #   Warning icon
          color = "red",                                                                     #   Red color
          fill = TRUE                                                                        #   Filled style
        )
      })
      return()                                                                              #   Abort
    }
    req(input$marker_boot)                                                                  # Require that a marker has been selected
    if (!(input$marker_boot %in% x_cols)) {                                                 # Validate the selected marker
      output$ibox_boot <- renderInfoBox({                                                   #   If invalid, show error info box
        shinydashboard::infoBox(
          "Bootstrap power",                                                                #   Title
          "—",                                                                              #   Placeholder
          "Selected marker not found in X_* columns.",                                      #   Message
          icon = icon("triangle-exclamation"),                                              #   Warning icon
          color = "red",                                                                     #   Red color
          fill = TRUE                                                                        #   Filled style
        )
      })
      return()                                                                              #   Abort
    }
    n_min <- as.integer(input$n_min_boot)                                                   # Read minimum n for the grid
    n_max <- as.integer(input$n_max_boot)                                                   # Read maximum n for the grid
    n_step <- as.integer(input$n_step_boot)                                                 # Read step size for the grid
    if (is.na(n_min) ||                                                                     # Validate grid settings:
        is.na(n_max) || is.na(n_step) || n_max < n_min || n_step <= 0) {                    #   non-NA, monotonic, positive step
      output$ibox_boot <- renderInfoBox({                                                   #   If invalid, show error info box
        shinydashboard::infoBox(
          "Bootstrap power",                                                                #   Title
          "—",                                                                              #   Placeholder
          "Invalid n-grid settings.",                                                       #   Message
          icon = icon("triangle-exclamation"),                                              #   Warning icon
          color = "red",                                                                     #   Red color
          fill = TRUE                                                                        #   Filled style
        )
      })
      return()                                                                              #   Abort
    }
    n_grid <- seq(n_min, n_max, by = n_step)                                                # Construct the sequence of n values
    dfp <- try(bootstrap_power(                                                             # Try to compute the bootstrap power curve
      dat,                                                                                  #   Data
      model = input$boot_model,                                                             #   SIM or MMM
      alpha = input$alpha_boot,                                                             #   Alpha level
      n_grid = n_grid,                                                                      #   n grid
      B = as.integer(input$B_boot),                                                         #   Number of resamples
      marker = input$marker_boot,                                                           #   Marker under test
      bonferroni = isTRUE(input$bonf_boot)                                                  #   Apply Bonferroni if checked
    ),
    silent = TRUE)                                                                          # Suppress errors for custom handling
    if (inherits(dfp, "try-error")) {                                                       # If an error occurred in bootstrap_power
      msg <- conditionMessage(attr(dfp, "condition"))                                       #   Extract the error message
      output$ibox_boot <- renderInfoBox({                                                   #   Show it in a red info box
        shinydashboard::infoBox(
          "Bootstrap power",                                                                #   Title
          "—",                                                                              #   Placeholder
          paste("Error:", msg),                                                             #   Display the error text
          icon = icon("triangle-exclamation"),                                              #   Warning icon
          color = "red",                                                                     #   Red color
          fill = TRUE                                                                        #   Filled style
        )
      })
      return()                                                                              #   Abort
    }
    target <- input$target_power_boot                                                       # Read the target power threshold
    n_star <- NA_integer_                                                                   # Initialize minimal n achieving target
    for (i in seq_len(nrow(dfp)))                                                           # Iterate over the power rows
      if (is.finite(dfp$power[i]) &&                                                        #   If power is valid and
          dfp$power[i] >= target) {                                                         #   meets/exceeds the target
        n_star <- dfp$n[i]                                                                  #     Record the corresponding n
        break                                                                               #     Stop at the first success (smallest n in grid)
      }
    if (is.na(n_star)) {                                                                    # If target power was not reached
      output$ibox_boot <- renderInfoBox({                                                   #   Show a red info box explaining that
        shinydashboard::infoBox(
          "Bootstrap power",                                                                #   Title
          "—",                                                                              #   Placeholder
          "Target power not reached within the n-grid.",                                    #   Message
          icon = icon("triangle-exclamation"),                                              #   Warning icon
          color = "red",                                                                     #   Red color
          fill = TRUE                                                                        #   Filled style
        )
      })
    } else {                                                                                # Otherwise, target power was achieved
      output$ibox_boot <- renderInfoBox({                                                   #   Show a blue info box with the n*
        shinydashboard::infoBox(
          "Bootstrap power",                                                                #   Title
          paste0("n* approx ", n_star),                                                     #   Display minimal grid n meeting target
          sprintf("Target power %.2f reached within grid.", target),                        #   Confirm the target was reached
          icon = icon("chart-line"),                                                        #   Chart icon
          color = "blue",                                                                    #   Informational blue
          fill = TRUE                                                                        #   Filled style
        )
      })
    }
    output$boot_plot <- renderPlot({                                                        # Render the power curve plot
      ggplot(dfp, aes(n, power)) + geom_line() + geom_point() +                             #   Draw line + points for power vs n
        geom_hline(yintercept = target, linetype = "dashed") +                              #   Add dashed horizontal line at target power
        labs(                                                                               #   Axis labels and title
          y = "Power (reject proportion)",                                                  #     Y-axis label
          x = "n",                                                                          #     X-axis label
          title = sprintf(                                                                  #     Title including settings
            "Bootstrap power (%s; marker = %s; B = %d)",
            input$boot_model,
            input$marker_boot,
            input$B_boot
          )
        ) +
        theme_minimal()                                                                     #   Clean minimal theme
    })
    output$boot_table <- renderDT({                                                         # Render a table of the power results
      datatable(transform(dfp, power = round(power, 3)),                                    #   Round power values for display
                options = list(pageLength = 10))                                            #   Show 10 rows per page
    })
    boot_cache(                                                                             # Cache results for report/export
      list(
        data = dfp,                                                                         #   The power data frame
        target = target,                                                                    #   Target power used
        model = input$boot_model,                                                           #   Model type
        marker = input$marker_boot,                                                         #   Tested marker
        B = as.integer(input$B_boot)                                                        #   Number of bootstrap resamples
      )
    )
  })
  
  # ---- FDR helpers ----                                                                      # Section header for functions related to FDR computation
  
  # Compute p-values for each marker under given model                                          # Function comment: returns per-marker p-values for SIM/MMM
  pvals_for_data <- function(dat, model = c("single", "multi")) {                              # Define function with data frame and model type
    model <- match.arg(model)                                                                  # Validate and select the model argument
    x_cols <- names(dat)[startsWith(names(dat), "X_")]                                         # Collect predictor columns whose names start with "X_"
    if (length(x_cols) == 0 ||                                                                 # If there are no X_* columns
        !"y" %in% names(dat))                                                                  # or response y is missing
      stop("Data must contain y and X_*.", call. = FALSE)                                      #   Abort with a clear error message
    
    if (model == "single") {                                                                   # Branch: single-marker (SIM) model
      ps <- sapply(x_cols, function(v) {                                                       #   For each marker column compute its p-value
        fit <- try(lm(reformulate(v, "y"), data = dat), silent = TRUE)                         #     Try to fit y ~ X_v, suppressing errors
        if (inherits(fit, "try-error"))                                                        #     If the fit failed
          return(NA_real_)                                                                     #       Return NA for this marker
        co <- summary(fit)$coefficients                                                        #     Extract coefficient table
        if (nrow(co) < 2)                                                                      #     If slope row is not available
          return(NA_real_)                                                                     #       Return NA
        co[2, 4]                                                                               #     Return p-value of the marker (second row, 4th column)
      })
    } else {                                                                                   # Branch: multiple-marker (MMM) model
      if (nrow(dat) <= length(x_cols) + 1)                                                     #   Check that n > p + 1 for full model
        stop("For MMM, need n > p + 1.", call. = FALSE)                                        #   Abort if not enough observations
      fit <- try(lm(reformulate(x_cols, "y"), data = dat), silent = TRUE)                      #   Fit y ~ all X_* simultaneously
      if (inherits(fit, "try-error"))                                                          #   If the model fit failed
        stop("MMM fit failed.", call. = FALSE)                                                 #     Abort with error
      sm <- summary(fit)$coefficients                                                           #   Extract coefficient summary table
      ps <- rep(NA_real_, length(x_cols))                                                      #   Initialize p-value vector with NAs
      names(ps) <- x_cols                                                                      #   Name the vector by marker names
      ok <- intersect(rownames(sm), x_cols)                                                    #   Identify coefficient rows that match markers
      ps[ok] <- sm[ok, 4]                                                                      #   Fill p-values for available markers
    }
    as.numeric(ps)                                                                              # Return p-values as a numeric vector
  }
  
  # Bootstrap FDR curve estimator                                                               # Function comment: estimates FDR across n using bootstrap
  fdr_curve_bootstrap <- function(dat,                                                         # Define function
                                  model = c("single", "multi"),                                #   Model type (SIM/MMM)
                                  q = 0.05,                                                    #   BH q-value threshold
                                  n_grid,                                                      #   Vector of sample sizes to evaluate
                                  B = 1000,                                                     #   Number of bootstrap resamples per n
                                  seed = 1,                                                    #   RNG seed for reproducibility
                                  use_perm = TRUE) {                                           #   If TRUE, use permutation to estimate false discoveries
    # add option for permutation-based FDR                                                     # Note: permutation option added to avoid flat ~0 curves
    stopifnot(is.data.frame(dat))                                                              # Ensure input is a data frame
    model <- match.arg(model)                                                                  # Validate model argument
    set.seed(seed)                                                                             # Set RNG seed for reproducibility
    
    # --- Input validation                                                                     # Block: sanity checks for inputs
    x_cols <- names(dat)[startsWith(names(dat), "X_")]                                         # Get marker columns
    if (!length(x_cols))                                                                       # If none found
      stop("'X_*' columns required.", call. = FALSE)                                           #   Abort with message
    if (!("y" %in% names(dat)))                                                                # If response is missing
      stop("Column 'y' required.", call. = FALSE)                                              #   Abort with message
    if (model == "multi" && nrow(dat) <= length(x_cols) + 1)                                   # If MMM but n too small
      stop("For MMM, need n > p + 1.", call. = FALSE)                                          #   Abort with message
    
    # --- If user selects 'Use permutation', ignore ground truth                               # Decide whether to use true causal info or permutation
    truth <- if (isTRUE(use_perm))                                                             # If permutation estimation is requested
      NULL                                                                                     #   Do not use ground truth
    else                                                                                       # Otherwise
      get_truth_index()                                                                        #   Retrieve causal/noncausal indices from simulation
    
    # --- Output structure                                                                     # Prepare results container
    out <- data.frame(                                                                         # Create data frame to store summary per n
      n = n_grid,                                                                              #   Sample size values
      mean_selected = NA_real_,                                                                #   Mean number of discoveries S
      sd_selected = NA_real_,                                                                  #   SD of number of discoveries S
      mean_FDP = NA_real_,                                                                     #   Mean false discovery proportion
      sd_FDP = NA_real_                                                                        #   SD of false discovery proportion
    )
    
    # --- Helper: compute p-values for one dataset                                             # Local helper to compute p-values for a given dataset
    pvals_for_data <- function(d) {                                                            # Define helper
      if (model == "single") {                                                                 #   SIM case
        ps <- vapply(x_cols, function(v) {                                                     #     Loop over markers
          fit <- try(lm(reformulate(v, "y"), data = d), silent = TRUE)                         #       Fit y ~ X_v
          if (inherits(fit, "try-error"))                                                      #       If model fails
            return(NA_real_)                                                                   #         Return NA
          co <- summary(fit)$coefficients                                                      #       Extract coefficients
          if (nrow(co) < 2)                                                                    #       If slope row missing
            return(NA_real_)                                                                   #         Return NA
          co[2, 4]                                                                             #       Return marker p-value
        }, numeric(1))                                                                         #     Ensure numeric output of length 1 per marker
        names(ps) <- x_cols                                                                    #     Name vector by marker names
        ps                                                                                     #     Return SIM p-values
      } else {                                                                                 #   MMM case
        fit <- try(lm(reformulate(x_cols, "y"), data = d), silent = TRUE)                      #     Fit y ~ all X_* together
        if (inherits(fit, "try-error"))                                                        #     If fit fails
          return(rep(NA_real_, length(x_cols)))                                                #       Return all NA p-values
        sm <- summary(fit)$coefficients                                                        #     Extract coefficient summary
        ps <- rep(NA_real_, length(x_cols))                                                    #     Init p-vector with NAs
        names(ps) <- x_cols                                                                    #     Name by marker
        ok <- intersect(rownames(sm), x_cols)                                                  #     Identify overlapping coefficient rows
        ps[ok] <- sm[ok, 4]                                                                    #     Fill in available p-values
        ps                                                                                     #     Return MMM p-values
      }
    }
    
    # --- Bootstrap loop                                                                       # Main bootstrap over n values
    for (i in seq_along(n_grid)) {                                                             # Iterate across sample sizes
      n_i <- n_grid[i]                                                                         #   Current n
      sel_vec <- numeric(B)                                                                    #   Vector to store number of selections S per resample
      fdp_vec <- numeric(B)                                                                    #   Vector to store FDP per resample
      
      for (b in seq_len(B)) {                                                                  #   For each bootstrap resample
        idx <- sample.int(nrow(dat), size = n_i, replace = TRUE)                               #     Sample n_i rows with replacement
        d2 <- dat[idx, , drop = FALSE]                                                         #     Create bootstrap dataset
        
        ps <- pvals_for_data(d2)                                                               #     Compute marker p-values on bootstrap data
        qv <- p.adjust(ps, method = "BH")                                                      #     Convert to BH q-values
        keep <- which(is.finite(qv) & qv <= q)                                                 #     Indices of discoveries under threshold q
        S <- length(keep)                                                                       #     Number of discoveries
        sel_vec[b] <- S                                                                         #     Store S
        
        if (S == 0) {                                                                           #     If no discoveries
          fdp_vec[b] <- 0                                                                       #       FDP is 0 by convention
          next                                                                                  #       Continue to next resample
        }
        
        # Estimate false discoveries V                                                          #     Compute false discoveries
        if (!is.null(truth)) {                                                                  #     If ground truth is available
          sel_idx <- match(names(ps)[keep], x_cols)                                             #       Map selected names to column indices
          V <- sum(sel_idx %in% truth$noncausal)                                                #       Count how many selections are noncausal
          fdp_vec[b] <- V / S                                                                    #       FDP = V/S
        } else {                                                                                #     Otherwise use permutation-based estimate
          # Permutation-based FDR estimate                                                      #       Build a null dataset by permuting y
          d_null <- d2                                                                          #       Copy bootstrap dataset
          d_null$y <- sample(d_null$y)                                                          #       Permute y to break associations
          ps0 <- pvals_for_data(d_null)                                                         #       Compute p-values under the null
          qv0 <- p.adjust(ps0, method = "BH")                                                   #       Convert to BH q-values
          Vhat <- sum(is.finite(qv0) & qv0 <= q)                                                #       Estimated false positives at the same threshold
          fdp_vec[b] <- Vhat / S                                                                 #       FDP estimate = Vhat / S
        }
      }
      
      out$mean_selected[i] <- mean(sel_vec)                                                     # Store mean number of selections across resamples
      out$sd_selected[i]   <- stats::sd(sel_vec)                                                # Store SD of number of selections
      out$mean_FDP[i]      <- mean(fdp_vec)                                                     # Store mean FDP across resamples
      out$sd_FDP[i]        <- stats::sd(fdp_vec)                                                # Store SD of FDP across resamples
    }
    
    out                                                                                         # Return the summary data frame
  }
  
  # Observe FDR computation (integrated fully)                                                  # Observer: run FDR bootstrap when user clicks the button
  observeEvent(input$btn_fdr, {                                                                 # Trigger on "Run FDR curve"
    dat <- current_data()
    if (is.null(dat) || !("y" %in% names(dat)) ||
        !any(startsWith(names(dat), "X_")) || anyNA(dat)) {
      output$fdr_status <- renderUI(tags$div(class = "text-danger", bad_data_msg))
      output$fdr_plot  <- renderPlot({ plot.new(); title("No data") })
      output$fdr_table <- renderDT(DT::datatable(data.frame()))
      return()
    }
    
    
    # --- Validate input grid                                                                   # Validate the n-grid inputs
    n_min <- as.integer(input$fdr_n_min)                                                        # Minimum n from UI
    n_max <- as.integer(input$fdr_n_max)                                                        # Maximum n from UI
    n_step <- as.integer(input$fdr_n_step)                                                      # Step size from UI
    if (is.na(n_min) ||                                                                         # Check for NA values
        is.na(n_max) || is.na(n_step) || n_max < n_min || n_step <= 0) {                        # Check order and positivity
      output$fdr_status <- renderUI(tags$div(class = "text-danger", "FDR curve — invalid n-grid."))  #   Show validation error
      return()                                                                                  #   Abort
    }
    n_grid <- seq(n_min, n_max, by = n_step)                                                    # Build the sequence of n values
    
    withProgress(message = "Running FDR curve ...", value = 0, {                                # Show progress bar during computation
      incProgress(0.2)                                                                           #   Advance progress
      
      # --- Run FDR bootstrap                                                                    # Try to compute the FDR curve
      res <- try(fdr_curve_bootstrap(                                                            #   Call the bootstrap routine
        dat,                                                                                     #   Data
        model    = input$fdr_model,                                                              #   SIM or MMM
        q        = input$fdr_q,                                                                  #   BH threshold q
        n_grid   = n_grid,                                                                       #   n grid
        B        = as.integer(input$fdr_B),                                                      #   Number of resamples
        seed     = as.integer(input$fdr_seed),                                                   #   RNG seed
        use_perm = isTRUE(input$fdr_use_perm)                                                    #   Use permutation if checkbox is ticked
      ),
      silent = TRUE)                                                                             #   Suppress errors (we handle below)
      
      incProgress(0.7)                                                                           #   Advance progress
      
      if (inherits(res, "try-error")) {                                                          #   If the bootstrap routine failed
        msg <- conditionMessage(attr(res, "condition"))                                          #     Extract error message
        output$fdr_status <- renderUI(tags$div(class = "text-danger", paste("FDR curve — Error:", msg))) #   Show error in status
        output$fdr_plot <- renderPlot({                                                          #     Render placeholder plot
          plot.new()                                                                             #       Empty plot
          title("Error")                                                                         #       Title indicating error
        })
        output$fdr_table <- renderDT(DT::datatable(data.frame()))                                #     Empty results table
        return()                                                                                 #     Abort
      }
      
      # --- Display status info                                                                  # If success, show summary line above outputs
      output$fdr_status <- renderUI(tags$div(HTML(                                               #   Render small HTML status block
        sprintf(                                                                                 #   Compose text with model, B, q, and mode used
          "FDR curve — done.<br><small>Model = <b>%s</b>, B = <b>%d</b>, q = <b>%.3f</b>%s</small>",
          input$fdr_model,
          as.integer(input$fdr_B),
          input$fdr_q,
          if (isTRUE(input$fdr_use_perm))                                                        #   Indicate whether permutation mode was used
            " (permutation-based FDP)"
          else
            " (truth-based FDP if available)"
        )
      )))
      
      # --- Plot FDR curve                                                                       # Draw the FDR curve
      output$fdr_plot <- renderPlot({                                                            #   Render plot into UI
        ggplot(res, aes(n, mean_FDP)) +                                                          #   Map n to x-axis and mean_FDP to y-axis
          geom_line(color = "blue") +                                                            #   Line for the curve
          geom_point(color = "blue") +                                                           #   Points on the curve
          geom_hline(                                                                            #   Horizontal reference line at q
            yintercept = input$fdr_q,
            linetype = "dashed",
            color = "red"
          ) +
          labs(                                                                                  #   Titles and axis labels
            title = sprintf(
              "FDR curve (BH, %s model, B = %d, q = %.2f)",
              input$fdr_model,
              input$fdr_B,
              input$fdr_q
            ),
            x = "Sample size (n)",                                                               #   X-axis label
            y = "Estimated FDR (E[FDP])"                                                         #   Y-axis label
          ) +
          theme_minimal()                                                                        #   Minimal theme for clarity
      })
      
      # --- Show result table                                                                     # Render the numeric results table
      output$fdr_table <- renderDT({                                                             #   Create DT table
        DT::datatable(                                                                           #   Build datatable object
          transform(                                                                             #   Round numeric columns for readability
            res,
            mean_selected = round(mean_selected, 3),
            sd_selected   = round(sd_selected, 3),
            mean_FDP      = round(mean_FDP, 3),
            sd_FDP        = round(sd_FDP, 3)
          ),
          options = list(pageLength = 10)                                                        #   Show 10 rows per page
        )
      })
      
      fdr_cache(list(                                                                            # Cache results for the report tab
        data = res,                                                                              #   FDR curve data
        q = input$fdr_q,                                                                         #   q threshold used
        model = input$fdr_model,                                                                 #   model used
        B = as.integer(input$fdr_B)                                                              #   number of resamples
      ))
      
      incProgress(0.1)                                                                           # Finish progress bar
    })
  })
  
  # ---- report build & download ----                                                           # Section header for report generation and download
  report_path <- reactiveVal(NULL)                                                              # Reactive holder for the path to the generated HTML report
  
  observeEvent(input$btn_build_report, {                                                        # When the "Build report" button is clicked, run this block
    dat <- current_data()                                                                       #   Get the current active dataset (y + X_*)
    max_rows <- input$report_max_rows %||% 50                                                   #   Max rows to show in each table (fallback to 50)
    sections <- input$report_sections %||% character(0)                                         #   Sections selected to include in the report
    title <- input$report_title %||% "QTL Design Report"                                        #   Report title (fallback default)
    now <- Sys.time()                                                                           #   Timestamp for report generation time
    
    css <- HTML(                                                                                #   Inline CSS to style the exported HTML report
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
    
    blocks <- list(tags$h1(title), tags$p(class = "muted", sprintf(                            #   Start assembling a list of HTML blocks for the report body
      "Generated: %s", format(now, "%Y-%m-%d %H:%M:%S")                                        #     Add a subtitle with the generation time
    )))
    
    if ("data" %in% sections) {                                                                #   If the "data" section is requested
      blk <- list(tags$div(                                                                     #     Build the "Data" section block
        class = "section",                                                                      #     Apply section styling class
        tags$h2("Data"),                                                                        #     Section title
        if (is.null(dat))                                                                       #     If no data is available
          tags$em("No data available. Generate or upload on the 'data' tab.")                   #       Show hint to generate or upload
        else                                                                                    #     Otherwise show summary, preview, and y histogram
          tagList(
            tags$p(sprintf(                                                                     #       One-line summary of rows and number of markers
              "Rows: %d, Markers (p): %d", nrow(dat), sum(startsWith(names(dat), "X_"))
            )),
            tags$h3("Preview (head)"),                                                          #       Subheader for data preview
            df_to_table(dat, max_rows = min(max_rows, 10)),                                     #       Render first rows of the dataset as an HTML table
            {                                                                                   #       Inline plotting block (histogram of y)
              if ("y" %in% names(dat)) {                                                        #         Only if y exists
                plt <- ggplot(dat, aes(x = y)) + geom_histogram(bins = 30) + theme_minimal() +  #           Build histogram ggplot
                  labs(title = "Histogram of y (active)", x = "y")                              #           Add labels
                img <- plot_to_base64(plt)                                                      #           Render plot to base64 PNG
                tags$div(                                                                       
                  tags$img(src = img, style = "max-width:100%;height:auto;border:1px solid #eee") #         Embed the image into the report
                )
              } else                                                                            #         If no y column
                tags$em("No 'y' column to plot.")                                              #           Show an explanatory message
            },
            {                                                                                   #       Optional realized PVE text if present
              pve_attr <- attr(y_gen(), "pve")                                                  #         Get realized PVE from y_gen attributes
              sigma_attr <- attr(y_gen(), "sigma")                                              #         Get sigma attribute if present
              if (!is.null(pve_attr) && is.finite(pve_attr))                                    #         If PVE is available and finite
                tags$p(HTML(                                                                     #           Render a line with PVE and sigma
                  sprintf(
                    "Realized sample PVE approx <b>%.3f</b> (sigma approx %.3f).",
                    pve_attr,
                    ifelse(is.null(sigma_attr), NA_real_, sigma_attr)
                  )
                ))
            }
          )
      ))
      blocks <- c(blocks, list(blk))                                                            #     Append the "Data" section to blocks
    }
    
    if ("SIM" %in% sections) {                                                                  #   If the "SIM" section is requested
      df_sim <- sim_cache()                                                                     #     Try to use cached SIM results
      if (is.null(df_sim) &&                                                                    #     If not cached but data is available
          !is.null(dat) &&                                                                      #     and y and X_* exist
          "y" %in% names(dat) && any(startsWith(names(dat), "X_"))) {
        df_sim <- try(compute_sim_stats(dat), silent = TRUE)                                    #       Compute SIM stats on the fly
        ; if (inherits(df_sim, "try-error"))                                                    #       If computation fails
          df_sim <- NULL                                                                        #         Fallback to NULL
      }
      blk <- list(tags$div(                                                                     #     Build "SIM" section block
        class = "section",                                                                      #     Apply section style
        tags$h2("Single-marker (SIM) results"),                                                 #     Section title
        if (is.null(df_sim) ||                                                                  #     If no results or empty
            !nrow(df_sim))
          tags$em(                                                                              #       Show instruction to compute
            "Not computed yet. Please click 'Fit model' on the 'model' tab."
          )
        else                                                                                    #     Otherwise render table and plot
          tagList(
            tags$p(HTML(                                                                         #       Show the threshold settings used
              sprintf(
                "Threshold on <b>%s</b> = %.3g%s",
                input$thr_stat,
                input$thr_value,
                if (identical(input$thr_stat, "p") &&                                           #       Mention Bonferroni if applicable
                    isTRUE(input$thr_bonf))
                  " (Bonferroni)"
                else
                  ""
              )
            )),
            df_to_table(cbind(marker = rownames(df_sim), df_sim), max_rows = max_rows),         #       Table of SIM results
            {                                                                                   #       Inline SIM plot block
              plt <- plot_stat(                                                                  #         Build a plot using the selected stat and threshold
                df_sim,
                stat = input$thr_stat,
                thr = input$thr_value,
                bonf = isTRUE(input$thr_bonf),
                title_prefix = "SIM ·"
              )
              img <- plot_to_base64(plt)                                                        #         Render to base64 image
              tags$div(
                tags$img(src = img, style = "max-width:100%;height:auto;border:1px solid #eee") #         Embed in report
              )
            }
          )
      ))
      blocks <- c(blocks, list(blk))                                                            #     Append the "SIM" section
    }
    
    if ("MMM" %in% sections) {                                                                  #   If the "MMM" section is requested
      df_mmm <- mmm_cache()                                                                     #     Try to use cached MMM results
      if (is.null(df_mmm) &&                                                                    #     If not cached, compute if feasible
          !is.null(dat) && any(startsWith(names(dat), "X_")) &&
          nrow(dat) > sum(startsWith(names(dat), "X_")) + 1) {
        df_mmm <- try(compute_mmm_stats(dat), silent = TRUE)                                    #       Compute MMM stats
        ; if (inherits(df_mmm, "try-error"))                                                    #       On failure
          df_mmm <- NULL                                                                        #         Set to NULL
      }
      blk <- list(tags$div(                                                                     #     Build "MMM" section block
        class = "section",
        tags$h2("Multiple-marker (MMM) results"),                                               #     Section title
        if (is.null(df_mmm) ||                                                                  #     If unavailable or empty
            !nrow(df_mmm))
          tags$em("Not computed or insufficient n > p + 1.")                                   #       Explain why it's missing
        else                                                                                    #     Otherwise render content
          tagList(
            tags$p(HTML(                                                                         #       Show threshold settings used
              sprintf(
                "Threshold on <b>%s</b> = %.3g%s",
                input$thr_stat,
                input$thr_value,
                if (identical(input$thr_stat, "p") &&
                    isTRUE(input$thr_bonf))
                  " (Bonferroni)"
                else
                  ""
              )
            )),
            df_to_table(cbind(marker = rownames(df_mmm), df_mmm), max_rows = max_rows),         #       Table of MMM results
            {                                                                                   #       Inline MMM plot block
              plt <- plot_stat(
                df_mmm,
                stat = input$thr_stat,
                thr = input$thr_value,
                bonf = isTRUE(input$thr_bonf),
                title_prefix = "MMM ·"
              )
              img <- plot_to_base64(plt)                                                        #         Render to base64 image
              tags$div(
                tags$img(src = img, style = "max-width:100%;height:auto;border:1px solid #eee") #         Embed plot
              )
            }
          )
      ))
      blocks <- c(blocks, list(blk))                                                            #     Append the "MMM" section
    }
    
    if ("power" %in% sections) {                                                                #   If the "power" section is requested
      lam <- lambda_cache()                                                                     #     Retrieve cached lambda-method results
      boot <- boot_cache()                                                                       #     Retrieve cached bootstrap power results
      fdrc <- fdr_cache()                                                                        #     Retrieve cached FDR curve results
      blk <- list(                                                                               #     Build the "Power" section
        tags$div(
          class = "section",
          tags$h2("Power"),                                                                      #     Main power heading
          tags$h3("Lambda / NCP method"),                                                        #     Subsection for lambda method
          if (is.null(lam))                                                                      #     If lambda results absent
            tags$em(
              "Not computed yet. Use the 'closure' method on the 'power' tab."
            )
          else                                                                                   #     Otherwise, show parameter summary and results
            tagList(
              tags$p(HTML(
                sprintf(
                  "Marker: %s <span class='pill'>p = %d</span>",                                 #       Show marker label and p (number of predictors)
                  lam$marker_label,
                  lam$p
                )
              )),
              tags$p(HTML(
                sprintf(                                                                          #       Show alpha, beta, sigma, VIF used
                  "alpha = %.3g; beta_j = %.3f; sigma = %.3f; VIF = %.3f",
                  lam$alpha,
                  lam$beta,
                  lam$sigma,
                  lam$vif
                )
              )),
              tags$p(                                                                             #       Show power at the given n
                sprintf(
                  "Given n = %d -> power approx %.3f",
                  lam$n_given,
                  lam$power_given
                )
              ),
              if (is.na(lam$n_star))                                                              #       If minimal n* not found
                tags$em("Target power unattainable under settings.")                              #         Show a warning
              else
                tags$p(HTML(                                                                       #       Otherwise show n* and bracket
                  sprintf(
                    "Target power = %.2f -> minimal n*: <b>%d</b> (bracket [%d, %d])",
                    lam$target,
                    lam$n_star,
                    lam$bracket[1],
                    lam$bracket[2]
                  )
                ))
            ),
          tags$h3("Bootstrap method"),                                                            #     Subsection for bootstrap power
          if (is.null(boot))                                                                      #     If no bootstrap results
            tags$em("Not computed yet. Run bootstrap on the 'power' tab.")                        #       Show hint
          else
            tagList(
              tags$p(                                                                             #       Summary line for bootstrap settings
                sprintf(
                  "Model = %s; Marker = %s; B = %d; target = %.2f",
                  boot$model,
                  boot$marker,
                  boot$B,
                  boot$target
                )
              ),
              df_to_table(transform(boot$data, power = round(power, 3)), max_rows = max_rows),    #       Table of (n, power)
              {                                                                                   #       Inline bootstrap power plot
                plt <- ggplot(boot$data, aes(n, power)) + geom_line() + geom_point() +            #         Line + points
                  geom_hline(yintercept = boot$target, linetype = "dashed") +                     #         Target power reference line
                  labs(
                    y = "Power",
                    x = "n",
                    title = sprintf("Bootstrap power (%s; %s)", boot$model, boot$marker)
                  ) +
                  theme_minimal()
                img <- plot_to_base64(plt)                                                        #         Render to base64 image
                tags$div(
                  tags$img(src = img, style = "max-width:100%;height:auto;border:1px solid #eee") #         Embed plot
                )
              }
            ),
          tags$h3("FDR curve (BH)"),                                                              #     Subsection for FDR curve
          if (is.null(fdrc))                                                                      #     If FDR results absent
            tags$em("Not computed yet.")                                                          #       Show hint
          else
            tagList(df_to_table(                                                                  #       Table of FDR curve summary
              transform(
                fdrc$data,
                mean_selected = round(mean_selected, 3),
                sd_selected   = round(sd_selected, 3),
                mean_FDP      = round(mean_FDP, 3),
                sd_FDP        = round(sd_FDP, 3)
              ),
              max_rows = max_rows
            ), {                                                                                  #       Inline FDR curve plot
              plt <- ggplot(fdrc$data, aes(n, mean_FDP)) + geom_line() + geom_point() +           #         Curve of E[FDP] vs n
                geom_hline(yintercept = (fdrc$q %||% 0.05),
                           linetype = "dashed") +
                labs(
                  y = "Estimated E[FDP]",
                  x = "n",
                  title = sprintf("FDR curve (BH, %s, B=%d)", fdrc$model, fdrc$B)
                ) +
                theme_minimal()
              img <- plot_to_base64(plt)                                                          #         Render to base64
              tags$div(
                tags$img(src = img, style = "max-width:100%;height:auto;border:1px solid #eee")   #         Embed image
              )
            })
        )
      )
      blocks <- c(blocks, list(blk))                                                             #     Append the "Power" section
    }
    
    if ("session" %in% sections) {                                                               #   If the "session" section is requested
      ses <- paste(capture.output(utils::sessionInfo()), collapse = "\n")                        #     Capture sessionInfo() output as text
      blk <- list(tags$div(                                                                      #     Build the "Session info" section
        class = "section",
        tags$h2("Session info"),
        tags$div(class = "code", ses)                                                            #     Render in monospaced code-styled block
      ))
      blocks <- c(blocks, list(blk))                                                             #     Append the "Session info" section
    }
    
    html <- tags$html(tags$head(                                                                 #   Assemble the full HTML document
      tags$meta(charset = "utf-8"),                                                              #     UTF-8 meta tag
      tags$title(title),                                                                          #     Document title
      tags$style(css)                                                                             #     Inject CSS styles
    ),
    tags$body(do.call(tagList, blocks)))                                                         #   Body contains all collected blocks
    
    tf <- tempfile(fileext = ".html")                                                             #   Create a temporary file path for the report
    htmltools::save_html(html, file = tf, background = "white")                                   #   Write the HTML to disk
    report_path(tf)                                                                               #   Store the path in the reactive value
    file.copy(tf, file.path(report_dir, basename(tf)), overwrite = TRUE)                          #   Copy into the served static folder for preview
    output$report_preview_hint <- renderUI(tags$p(class = "muted", "Preview below shows the generated HTML.")) #   Small hint above preview
    output$report_preview <- renderUI({                                                           #   Render an iframe that loads the generated HTML
      tags$iframe(src = file.path("qtl_reports", basename(tf)),
                  style = "width:100%;height:800px;border:1px solid #ddd;")
    })
  })
  
  output$btn_download_html <- downloadHandler(                                                    # Download handler for the HTML report
    filename = function() {                                                                       #   Function to generate the download file name
      paste0((input$report_title %||% "QTL-Design-Report"),                                       #     Use title or default prefix
             "_",                                                                                 #     Separator
             format(Sys.time(), "%Y%m%d-%H%M%S"),                                                 #     Timestamp (YYYYMMDD-HHMMSS)
             ".html"                                                                              #     File extension
      )
    },
    content  = function(file) {                                                                    #   Function that writes the file to the download path
      pth <- report_path()                                                                         #     Get the last-built report path
      validate(need(                                                                               #     Ensure a report exists
        !is.null(pth) &&
          file.exists(pth),
        "Please click 'Build report' first."
      ))
      file.copy(pth, file, overwrite = TRUE)                                                       #     Copy the report to the requested download path
    }
  )
  
}

shinyApp(ui, server)                                                                            # Launch the Shiny app with the defined UI and server



