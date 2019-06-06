if (!require("drake")) {install.packages("drake"); library(drake)}

pkgconfig::set_config("drake::strings_in_dots" = "literals")

source("functions.R")

set.seed(1692)

N <- 25e5

plan <- drake_plan(

  full_dat = rbind(make_data(N = N, loglin = FALSE),
                   make_data(N = round(N / .75), loglin = TRUE)),

  res = sum_data(full_dat),

  vals_4_rmd = vals(full_dat),

  table1 = overall_res(res),

  table2a = full_dat %>%
    filter(loglin == "Unrestricted probabilities") %>%
    quant_rat(),
  table2b = full_dat %>%
    filter(loglin == "Unrestricted probabilities",
           qual_int == "No qualitative interaction") %>%
    quant_rat(),
  table2c = full_dat %>%
    filter(loglin != "Unrestricted probabilities") %>%
    quant_rat(),
  table2d = full_dat %>%
    filter(loglin != "Unrestricted probabilities",
           qual_int == "No qualitative interaction") %>%
    quant_rat(),

  table3 = resid_bias_add(full_dat),

  compare_correctedA = ratio_true_plot(
    full_dat,
    "Unrestricted probabilities",
    c("No interaction restriction", "No qualitative interaction"),
    "eFigure 1A. Probabilities drawn from uniform distribution"
  ),
  compare_correctedB = ratio_true_plot(
    full_dat,
    "Unrestricted probabilities",
    "No qualitative interaction",
    "eFigure 1A. Uniform probabilities, restricted to no qualitative interaction"
  ),
  compare_correctedC = ratio_true_plot(
    full_dat,
    "Log-linear model",
    c("No interaction restriction", "No qualitative interaction"),
    "eFigure 1C. Log-linear model for M"
  ),
  compare_correctedD = ratio_true_plot(
    full_dat,
    "Log-linear model",
    "No qualitative interaction",
    "eFigure 1D. Log-linear model, restricted to no qualitative interaction"
  ),

  resid_plotA = resid_bias(
    full_dat,
    "Unrestricted probabilities",
    c("No interaction restriction", "No qualitative interaction"),
    "eFigure 2A. Probabilities drawn from uniform distribution"
  ),
  resid_plotB = resid_bias(
    full_dat,
    "Unrestricted probabilities",
    "No qualitative interaction",
    "eFigure 2B. Uniform probabilities, restricted to no qualitative interaction"
  ),
  resid_plotC = resid_bias(
    full_dat,
    "Log-linear model",
    c("No interaction restriction", "No qualitative interaction"),
    "eFigure 2C. Log-linear model for M"
  ),
  resid_plotD = resid_bias(
    full_dat,
    "Log-linear model",
    "No qualitative interaction",
    "eFigure 2D. Log-linear model, restricted to no qualitative interaction"
  ),
  eAppendix = rmarkdown::render(
    knitr_in("eAppendix.Rmd"),
    output_file = file_out("eAppendix.pdf"),
    quiet = TRUE
  )
)

make(plan)
