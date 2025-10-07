# app.R — Shiny
# model 页：单标记（single marker）与多标记（multiple marker）两种线性模型。
# power 页：根据 individuals 与 markers（自变量个数）计算所需样本量（Bonferroni）。
# 参考《Mastering Shiny》基本结构（Ch.1）、输入/输出（Ch.2）、tab 布局（Ch.6）、文件上传（Ch.9）。

suppressPackageStartupMessages({
  #Starts a code block that suppresses startup messages from packages, preventing the usual “Loading package …” text from appearing in the console.
  library(shiny)
  #Load the Shiny package, which provides all the function needed to build interactive web application in R
  library(DT)
  #Load the DT package, which allow you to create and display interactive data tables in Shiny app
})

# -----------------------------
# Data generation
# -----------------------------
# Define a function called sim_clean_data that generates a clean, synthetic dataset.
sim_clean_data <- function(n = 200, # number of individuals, representing for row number
                           m = 50, # number of markers, representing for column number
                           maf = 0.3, # minor allele frequency, which means the probability of "success" when generating genotype values 0,1 or 2
                           n_causal = 3, # number of causal markers that truly affect the outcome
                           beta = 0.2, # effect size, which means the coefficient for causal markers
                           seed = 1) { # set random seed for reproductibility
  set.seed(seed) # set the random number seed, ensuring that the same “random” data can be regenerated every time the function is run with the same seed.
  G <- matrix(rbinom(n * m, 2, maf), nrow = n, ncol = m)
  # create a numeric matrix G of size n × m that represents genotype data: 
  # rbinom(n * m, 2, maf) generates n*m random numbers from a binomial distribution with two trials (0, 1, or 2) and probability maf
  # matrix(..., nrow = n, ncol = m) reshapes the vector into a matrix.
  colnames(G) <- paste0("X_", seq_len(m)) 
  # name the columns of G as "X_1", "X_2", ..., "X_m".These names identify the predictor variables.
   b <- rep(0, m) # create a coefficient vector b of length m, filled with zeros — initially, all markers have no effect.
  if (n_causal > 0)
    b[seq_len(min(n_causal, m))] <- beta
  # If the number of causal markers (n_causal) is greater than 0, assign the value beta to the first few elements of b.
  # which means that only the first n_causal markers are “truly causal” and have a real effect on the outcome y.
  
   y <- as.numeric(G %*% b + rnorm(n, 0, 1))
  # generate the outcome variable y:
  # G %*% b computes the linear combination of genotype matrix and effect sizes.
  # rnorm(n, 0, 1) adds normally distributed random noise (mean = 0, sd = 1).
  # as.numeric(...) converts the result into a numeric vector. y is the simulated phenotype or response variable.
  
  data.frame(y = y, G, check.names = FALSE)
  # combine y and all columns of G into a single data frame and returns it.
  # The first column is the outcome variable y.
  # The remaining columns are predictors X_1, X_2, …, X_m.
  # check.names = FALSE ensures the column names remain exactly as defined (X_1, etc.) without R changing them.
}

# Single-marker scanning（y ~ X_j）
single_marker_scan <- function(dat) {
  # define a function single_marker_scan that takes a data frame dat and performs a separate regression of y on each single predictor X_j.
  x_cols <- names(dat)[startsWith(names(dat), "X_")]
  # extract the names of all columns in dat that start with "X_"; these are treated as the marker columns.
  out <- lapply(x_cols, function(v) {
    # create a list by applying a function to each predictor name v in x_cols.
    fit <- lm(dat$y ~ dat[[v]])
    # fit a simple linear regression for the outcome y on the single predictor v.
    co <- summary(fit)$coefficients
    # compute the model summary and extracts the coefficient table (estimate, standard error, t-value, p-value).
    c(beta_hat = unname(co[2, 1]),
      t = unname(co[2, 3]),
      p = unname(co[2, 4]))
  })
  # build a named numeric vector containing:
  # beta_hat: the estimated coefficient for the predictor (row 2, column 1),
  # t: the t-statistic (row 2, column 3),
  # p: the p-value (row 2, column 4).
  # unname() drops row/column names from these extracted values.
  as.data.frame(do.call(rbind, out), row.names = x_cols)
}
# row-bind all the small vectors from out into a matrix, convert it to a data frame, 
# and set the row names to the corresponding predictor names x_cols. 
# return the table of results for all single-marker regressions.

# Sample size per group with Bonferroni + two-sample t approximation (effect size = Cohen’s d)
required_n_per_group <- function(markers, # marker, which is the number of independent statistical tests.
                                          # This determines how strong the Bonferroni correction is
# define required_n_per_group, which computes the required per-group sample size for a two-sample t-test
                                 alpha = 0.05, # The overall significance level for all tests combined.
                                              # After Bonferroni correction, each individual test uses alpha / markers.
                                 power = 0.8, # The desired statistical power (1 − β), the probability of correctly 
                                              # rejecting the null hypothesis when there truly is an effect.
                                 d = 0.3) { # The assumed effect size (Cohen’s d), representing the standardized mean difference between two groups
  markers <- max(1L, as.integer(markers)) # Let the markers to be an integer and ensures it’s at least 1.
  alpha_eff <- alpha / markers # apply Bonferroni: divide the overall alpha by the number of markers to get the per-test significance level.
  alpha_eff <- max(alpha_eff, 1e-12) # guard against extremely tiny or zero alpha due to large markers; enforces a lower bound.
  res <- power.t.test( # call power.t.test to compute the required per-group sample size for a two-sample, two-sided t-test:
    delta = d, # delta = d uses Cohen’s d as the mean difference in SD units
    sd = 1, # sd = 1, which is consistent with that scaling
    sig.level = alpha_eff, # sig.level = alpha_eff uses the Bonferroni-adjusted alpha
    power = power, # power = power targets the desired statistical power
    type = "two.sample", # specify a two-sample t-test, comparing means from two independent groups
    alternative = "two.sided" # specify a two-tailed test — detects differences in either direction (μ₁ > μ₂ or μ₁ < μ₂).
  )
  ceiling(max(2, res$n)) # take the result res$n, which is the computed per-group sample size
} # ensure it’s at least 2 participants per group to avoid nonsensical values
# round up to the next whole number using ceiling() because sample size must be an integer.

# -----------------------------
# UI
# -----------------------------
ui <- navbarPage(
  # define the UI object using navbarPage, which creates a top navigation bar with multiple tabs.
  title = "QTL Design Sandbox",
  # title: The title shown on the navbar
  
  # ---- data ----
  tabPanel("data", sidebarLayout(
    # creates a tab named “data”. 
    # sidebarLayout splits the page into a left sidebar and a right main area.
    sidebarPanel( # left sidebar starts.
      h4("Generate clean demo"), # A small header indicating the demo-data generator section.
      numericInput( #add a numeric input
        "n_ind", # input id n_ind
        "Number of individuals", # label “Number of individuals”
        200, # default value 200
        min = 20, # minimum 20
        step = 20 # step 20 (how much it changes per click)
      ),
      numericInput( #add a numeric input
        "n_mark", # input id n_mark
        "Number of markers", # label “Number of markers”
        50, # default value 50
        min = 5, # minimum 5
        step = 5 # step 5 (how much it changes per click)
      ),
      sliderInput( #add a slider input
        "maf", # input id n_maf
        "MAF (demo)", # label “MAF (demo)”
        min = 0.05, # minimum 0.05
        max = 0.5, # maximum 0.5
        value = 0.3, # default value 0.3
        step = 0.01 # step 0.01
      ),
      numericInput( #add a numeric input
        "n_causal", # input id n_causal
        "# causal markers", # label “causal markers”
        3, # default value 3
        min = 0, # minimum 0
        step = 1 #step 1 (how much it changes per click)
      ),
      sliderInput( #add a slider input
        "beta", # input id beta
        "Effect size per causal marker (β)", # label “Effect size per causal marker (β)”
        min = 0, # minimum 0
        max = 1, # maximum 1
        value = 0.2, # default value 0.2
        step = 0.05 #step 0.05 (how much it changes per click)
      ),
      actionButton("btn_gen", "Generate demo data"), # action button to trigger demo data generation
      tags$hr(), # horizontal rule separator
      h4("Or upload CSV (clean)"), # a header for the upload option
      fileInput("upload", "Upload CSV", accept = ".csv"), # file upload control (id upload) that accepts .csv. Uploaded data overrides demo data.
      helpText(
        "Columns: y, X_1 ... X_m (no missing values). Upload overrides demo."
      ) # help text explaining the required structure of uploaded data: one y column and predictors named X_1 … X_m, with no missing values.
    ), # end the sidebarPanel
    mainPanel(
      h4("Preview"), # header “Preview”.
      DTOutput("data_head"), # interactive preview table (top rows) of current data.
      tags$hr(), # create a separator line.
      verbatimTextOutput("data_info") # print summary like number of individuals, number of markers, and column names.
    )
  )), # close the data tab
  
  # ---- model ----
  tabPanel("model", sidebarLayout( # create a tab named “model” with a sidebar layout.
    sidebarPanel( # create a sidebarPanel
      radioButtons( # radioButtons with id model_kind to choose between two model types:
        "model_kind", # the internal input ID you read in the server (e.g., input$model_kind) to get or update the selection.
        "Model type", # user-visible label
        c(
          "Single-marker models (y ~ X_j)" = "single", # “single” = run one regression per marker (y ~ X_j)
          "Multiple-marker model (y ~ all X_*)" = "multi" # “multi” = run one regression with all markers (y ~ X_1 + … + X_m)
        ),
        selected = "single" # default selection is single.
      ),
      actionButton("btn_fit", "Fit model") # trigger the modeling action on the server.
    ),
    mainPanel( # create the main content area of the sidebarLayout on the model tab, which is the right-hand pane next to the sidebar.
      conditionalPanel(
        condition = "input.model_kind == 'single'", #start a conditional container that is shown only when the client-side expression is true.
        h4("Single-marker scan results"), # display a level-4 heading as a section title: “Single-marker scan results.”
        DTOutput("single_res") # creaet a placeholder for an interactive DataTable
      ),
      conditionalPanel( # start a conditional container that is shown only when the client-side expression is true.
        condition = "input.model_kind == 'multi'", # start a second conditional container that is shown when the radio selection is "multi".
        h4("Multiple-marker model summary: y ~ all X_*"), # display a heading for the multiple-marker model section, clarifying the model form (y ~ all X_*).
        verbatimTextOutput("multi_summary") # allocate a text area that renders preformatted, monospaced output
      ), # The server should provide output$multi_summary via renderPrint(...) or renderText(...).
      tags$hr(), # insert a horizontal rule, a dividing line to visually separate the model results from the hint section below.
      uiOutput("model_hint") # create a placeholder for dynamically generated UI.
    ) # close the mainPanel.
  )), # close the surrounding sidebarLayout,and the enclosing tabPanel("model", ...)
  
  # ---- power ----
  tabPanel("power", sidebarLayout( # creates a new tab in the navbar labeled “power”, and inside that tab uses a sidebarLayout 
    sidebarPanel( # start the left sidebar area that will contain user inputs and brief guidance.
      helpText(
        "Sample size by Bonferroni-adjusted two-sample t-test per marker."
      ),
      # display gray, non-interactive helper text explaining the approach: the app will compute sample size per marker using a two-sample t-test with Bonferroni correction.
      numericInput("alpha", "Overall α", 0.05, min = 1e-6, step = 0.01), # add a numeric input for the overall significance level:ID: "alpha"; Default: 0.05; Minimum: 1e-6;Step: 0.01
      sliderInput( # Adds a slider for desired statistical power (1−β)
        "target_power", # input id target_power
        "Desired power", # label "Desired power"
        min = 0.5, # minimum 0.5
        max = 0.99, # maximum 0.99
        value = 0.8, # default value 0.8
        step = 0.01 # step 0.01 (how much it changes per click)
      ),
      sliderInput( # add a slider for the assumed effect size (Cohen’s d):
        "cohen_d", # input id cohen_d
        "Assumed effect size (Cohen's d)", # label "Assumed effect size (Cohen's d)"
        min = 0.1, # minimum 0.1
        max = 1.0, # maximum 1.0
        value = 0.3, # default value 0.3
        step = 0.05 # step 0.05 (how much it changes per click)
      ),
      actionButton("btn_power", "Compute required n") # add a button to trigger the calculation: ID: "btn_power";Label: “Compute required n”.
    ),
    mainPanel( # begin the main panel (right side), where results will be displayed.
      h4("Required sample size (given current markers M and overall α)"), # Shows a header explaining the result
      uiOutput("power_text") # reserve a dynamic output region with ID "power_text"
    ) # close the mainPanel.
  )) # close the sidebarLayout
) # close the outer navbarPage if this is the final tab in the UI definition.



# -----------------------------
# SERVER # mark the beginning of the server logic
# -----------------------------
server <- function(input, output, session) { # define the server function that Shiny will run.
  # input: a read-only list of values coming from UI controls
  # output: a list you assign to, so UI outputs can display results.
  # session: an object representing the current user session
  rule_of_thumb_n <- function(p) ceiling(10 * max(1, p)) # define a small helper inside server:
  # return the heuristic n ≥ 10 × p (rounded up), with p forced to be at least 1.
  demo_data <- eventReactive(input$btn_gen, { # Creates a reactive expression that re-computes 
    # only when the button btn_gen is clicked. The result it returns will be cached until the next click.
    sim_clean_data( # when triggered, calls sim_clean_data() to generate clean demo data, pulling parameters directly from the UI controls:
      n = input$n_ind, # n: number of individuals
      m = input$n_mark, # m: number of markers
      maf = input$maf, # MAF slider
      n_causal = input$n_causal, # n_causal: number of causal markers,
      beta = input$beta, # effect size per causal marker
      seed = 1 # fixed to 1 for reproducibility
    )
  }, ignoreInit = TRUE) # prevent running this code on app startup; it will only run after the first button click.
  
  uploaded <- reactive({ #define a reactive that represents the uploaded CSV. It re-evaluates whenever input$upload changes.
    file <- input$upload # grab the upload object
    if (is.null(file)) # if nothing is uploaded
      return(NULL) # return NULL
    read.csv(file$datapath, check.names = FALSE)
  }) # read the uploaded CSV from its temporary path, preserving original column names.
  
  current_data <- reactive({ # create a unified data source reactive that the rest of the app will use.
    up <- uploaded() # if a CSV has been uploaded, prefer using it.
    if (!is.null(up))
      return(up)
    req(demo_data()) # otherwise, require that demo data exists and return it
  }) #  req() will show a “waiting” state if demo data hasn’t been generated yet
  
  # data page output
  output$data_head <- renderDT({ # define the DT table output for the preview area
    dat <- current_data() # name current_data() briefly as dat
    datatable(head(dat, 10), options = list(scrollX = TRUE, pageLength = 10))
  }) # fetch the current dataset and renders the first 10 rows as an interactive table with horizontal scrolling and a 10-row page length.
  
  output$data_info <- renderPrint({ # define a verbatim print output for the info area (verbatimTextOutput("data_info") in UI).
    dat <- current_data() # name current_data() briefly as dat
    p <- sum(startsWith(names(dat), "X_")) # get the current dataset and counts how many columns start with "X_", treated as markers.
    cat("Rows (individuals):", nrow(dat), "\n") # print the number of rows (n)
    cat("Markers (predictors):", p, "\n") # print the number of markers (p)
    cat("Columns:", paste(names(dat), collapse = ", ")) # print a comma-separated list of column names.
  })
  
  # -------- model --------
  single_res <- eventReactive(input$btn_fit, { # define a reactive that runs when the user clicks the “Fit model” button. 
    # It will compute single-marker scan results.
    req(current_data()) # ensure data is available
    dat <- current_data() # pull the data
    validate(need("y" %in% names(dat), "Missing column 'y'.")) # validation: require a y column; if missing, show the message and stops this render.
    validate(need(any(startsWith( # require at least one predictor named X_*
      names(dat), "X_"
    )), "No X_* columns found.")) # otherwise show the message and stops.
    single_marker_scan(dat) # if validations pass, run the single-marker scan function and returns its results
  }, ignoreInit = TRUE) # ignoreInit = TRUE prevents running at startup.\
  
  output$single_res <- renderDT({ # define the DT output for the single-marker results table (DTOutput("single_res") in UI).
    res <- single_res() # pull the scan results
    req(res) # require they exist.
    datatable( # create an interactive DataTable widget to display a data frame in the UI.
      transform( # begin transforming a data frame by modifying or adding columns
        cbind(marker = rownames(res), res), # column-bind a new marker column (from res row names) to the front of the res results table.
        p = signif(p, 4), # format the p column to four significant digits for cleaner display.
        beta_hat = signif(beta_hat, 4), # format the beta_hat (coefficient estimate) column to four significant digits.
        t = round(t, 3) # rounds the t-statistic column to three decimal places.
      ), # close the transform() call returning the formatted data frame.
      options = list(pageLength = 10) # set DataTable options so the table shows ten rows per page.
    ) # close the datatable() call and returns the widget to the output.
  }) # close the renderDT({ ... }) expression for the single-marker results output.
  
  output$multi_summary <- renderPrint({ # define a text/console-style output called multi_summary that will print model summary text.
    req(input$btn_fit > 0) # require that the “Fit model” button has been clicked at least once before proceeding.
    dat <- current_data() # retrieve the current dataset
    x_cols <- names(dat)[startsWith(names(dat), "X_")] # collect the names of all predictor columns that start with "X_"
    fit <- lm(reformulate(x_cols, "y"), data = dat) # fit a multiple linear regression y ~ X_1 + ... + X_p using all X_* columns as predictors.
    summary(fit) # print the regression summary (coefficients, standard errors, t-values, p-values, R², etc.) to the output.
  }) #close the renderPrint({ ... }) call for the multi-marker summary output.
  
  output$model_hint <- renderUI({ # define a dynamic UI output called model_hint that will render a small block of explanatory UI.
    dat <- current_data() # fetch the current dataset to compute scale information.
    req(dat) # ensure the dataset exists before proceeding
    p <- sum(startsWith(names(dat), "X_")) # count the number of predictors/markers by tallying X_* columns.
    n <- nrow(dat) # get the current sample size as the number of rows in the dataset.
    tagList( # start building a list of UI elements to return together.
      h4("Scale check / intuition"), # add a small section header describing what this block shows.
      p(sprintf("Current: n = %d, p = %d.", n, p)), # display the current scale of the model by printing sample size n and number of predictors p.
      p( # begin a paragraph containing a formatted rule-of-thumb message.
        sprintf( # prepare a formatted string with values interpolated into the text.
          "Rule-of-thumb: n \u2265 10 \u00D7 p \u2248 %d (for linear regression).",
          # specifie the message string using Unicode symbols for ≥, ×, and ≈.
          rule_of_thumb_n(p) # insert the heuristic total sample size suggestion computed as ceiling(10 * p).
        ) # Closes the sprintf() call that created the formatted string.
      ), # Closes the paragraph element that contained the rule-of-thumb text.
      p( # Starts another paragraph for a short disclaimer.
        "This is only a heuristic; see power page for Bonferroni-based sample size."
        # Explains that the rule-of-thumb is informal and points users to the power page for a formal calculation.
      ) # Closes the disclaimer paragraph
    ) # Closes the tagList() of UI elements to be rendered.
  }) # Closes the renderUI({ ... }) call for the model hint output.
  
  # -------- power --------
  output$power_text <- renderUI({ # Defines a dynamic UI output called power_text that will show the computed sample-size recommendation.
    input$btn_power # Registers a reactive dependency on the “Compute required n” button so this code runs after clicks.
    isolate({ # Runs the enclosed code without tracking additional reactive dependencies (besides the button).
      dat <- current_data() # Gets the current dataset for counting markers and current sample size.
      req(dat) # Ensures the dataset exists; otherwise, do not render.
      p <- sum(startsWith(names(dat), "X_")) # Counts the number of markers/predictors by detecting X_* columns.
      n_curr <- nrow(dat) # Stores the current total sample size.
      d <- input$cohen_d # Reads the assumed effect size (Cohen’s d) from the power tab controls.
      pow <- input$target_power # Reads the desired statistical power from the power tab controls.
      alp <- input$alpha # Reads the overall alpha from the power tab controls.
      n_per_group <- required_n_per_group( # Starts computing the per-group required sample size using your helper function.
        markers = p, # Supplies the number of multiple tests (markers) for Bonferroni correction.
        alpha = alp, # Supplies the overall alpha (before Bonferroni)
        power = pow, # Supplies the desired power.
        d = d # Supplies the assumed effect size (Cohen’s d).
      ) # Closes the required_n_per_group() call.
      n_total <- 2 * n_per_group # Calculates the total sample size assuming two groups of equal size.
      alp_eff <- alp / max(1L, p) # Computes the per-test Bonferroni-adjusted alpha (guarding against zero markers).
      tagList(p( # Begins a UI block and starts the first paragraph.
        sprintf( # Prepares a formatted string for the first line of the summary.
          "Markers (M) = %d; overall α = %.3g; Bonferroni-adjusted α' = %.3g.",
          # Specifies the template text including placeholders for M, α, and the adjusted α′.
          p, # Inserts the number of markers.
          alp, # Inserts the overall alpha.
          alp_eff # Inserts the Bonferroni-adjusted alpha.
        ) # Closes the sprintf() for the first summary line.
      ), p( # Closes the first paragraph and begins the second paragraph.
        sprintf("Assumed Cohen's d = %.2f; target power = %.2f.", d, pow)
        # Prints the assumed effect size and the desired power with two decimal places.
      ), h4( # Closes the second paragraph and begins a header (to emphasize the required n).
        sprintf( # Prepares a formatted string for the required sample sizes.
          "Required n (per group) ≈ %d; total n ≈ %d.",
          # Specifies the message showing per-group and total sample sizes.
          n_per_group, # Inserts the computed required per-group sample size.
          n_total # Inserts the computed total sample size.
        ) # Closes the sprintf() for the header line
      ), p( # Closes the header and opens a final paragraph.
        sprintf( # Starts formatting the line comparing required total n to current n.
          "Current data has n = %d. %s",
          # Template text: current sample size followed by a recommendation sentence.
          n_curr, # Inserts the current total sample size.
          if (n_curr >= n_total)
            "It meets the total-n suggestion."
          else
            "Consider increasing sample size."
          # Chooses a message based on whether the current total n meets the recommended total n.
        ) # Closes the sprintf() for the recommendation line.
      )) # Closes the final paragraph and the tagList() of UI elements.
    }) # Closes the isolate({ ... }) block.
  }) # Closes the renderUI({ ... }) for the power text output.
} # Closes the server function definition.

shinyApp(ui, server) # Launches the Shiny app by pairing the ui object with the server function.
