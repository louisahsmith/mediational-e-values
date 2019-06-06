if (!require("tidyverse")) {install.packages("tidyverse"); library(tidyverse)}
if (!require("xtable")) {install.packages("xtable"); library(xtable)}


# function for bounding factor from any two RRs
BF_func <- function(RRa, RRb) {
  RRa * RRb / (RRa + RRb - 1)
}

# function to create dataset with NIE, NDE (true and observed), bias, etc.
# N: number of samples to draw
# loglin: whether M should be drawn according to loglinear model
make_data <- function(N, loglin = FALSE) {

  # draw N sets of 14 probabilities from uniform(0, 1)
  dat <- data.frame(matrix(runif(N * 14), ncol = 14))

  # letter names correspond to probabilities
  names(dat) <- letters[1:14]

  dat$loglin <- "unrestricted"

  # force to follow a loglinear model
  if (loglin) {
    dat <- dat %>%
      mutate(
        # randomly choose one probability to throw out
        keep = runif(N),
        b = case_when(keep < .25 ~ d * c / e, TRUE ~ b),
        c = case_when(between(keep, .25, .5) ~ b * e / d, TRUE ~ c),
        d = case_when(between(keep, .5, .75) ~ b * e / c, TRUE ~ d),
        e = case_when(keep > .75 ~ d * c / b, TRUE ~ e),
        loglin = "loglinear"
      ) %>%
      select(-keep) %>%
      # remove any impossible probabilities
      filter(b < 1 & c < 1 & d < 1 & e < 1)
  }

  # make new variables
  dat %>%

    # 1 - prob for P(var = 0)
    mutate_if(is.numeric, funs(comp = 1 - .)) %>%

    mutate(

      # ratio of obs/true NDE
      obs_true = ((1 / (b * a + c * a_comp)) * (f * b * a + h * c * a_comp) *
        (d * a + e * a_comp) +
        (1 / (b_comp * a + c_comp * a_comp)) *
          (g * b_comp * a + i * c_comp * a_comp) * (d_comp * a + e_comp * a_comp)) /
        (f * d * a + g * d_comp * a + h * e * a_comp + i * e_comp * a_comp),

      # values of true and observed effects
      NIE_true = (f * b * a + g * b_comp * a + h * c * a_comp + i * c_comp * a_comp) /
        (f * d * a + g * d_comp * a + h * e * a_comp + i * e_comp * a_comp),
      NIE_obs = NIE_true / obs_true,
      NDE_true = (f * d * a + g * d_comp * a + h * e * a_comp + i * e_comp * a_comp) /
        (j * d * a + k * d_comp * a + l * e * a_comp + m * e_comp * a_comp),
      NDE_obs = NDE_true * obs_true,

      # make sure the bias > 1 so bound is sensical
      bias = (obs_true > 1) * obs_true + (obs_true < 1) / obs_true,

      # make loglin into a factor for later grouping
      loglin = fct_recode(loglin,
        "Unrestricted probabilities" = "unrestricted",
        "Log-linear model" = "loglinear"
      ),

      # create possible RRau values
      RRau_m1u1 = b * (d * a + e * a_comp) / (d * (b * a + c * a_comp)),
      RRau_m1u0 = c * (d * a + e * a_comp) / (e * (b * a + c * a_comp)),
      RRau_m0u1 = b_comp * (d_comp * a + e_comp * a_comp) /
        (d_comp * (b_comp * a + c_comp * a_comp)),
      RRau_m0u0 = c_comp * (d_comp * a + e_comp * a_comp) /
        (e_comp * (b_comp * a + c_comp * a_comp)),

      # create possible RRuy values
      RRuy_m1u1 = f / h,
      RRuy_m1u0 = 1 / RRuy_m1u1,
      RRuy_m0u1 = g / i,
      RRuy_m0u0 = 1 / RRuy_m0u1,

      # create possible RRam values
      RRamm_a1m1 = (b * a + c * a_comp) / (d * a + e * a_comp),
      RRamm_a0m1 = 1 / RRamm_a1m1,
      RRamm_a1m0 = (b_comp * a + c_comp * a_comp) / (d_comp * a + e_comp * a_comp),
      RRamm_a0m0 = 1 / RRamm_a1m0,

      # create possible RRam|u values
      RRamc_a1m1u1 = b / d,
      RRamc_a0m1u1 = 1 / RRamc_a1m1u1,
      RRamc_a1m0u1 = b_comp / d_comp,
      RRamc_a0m0u1 = 1 / RRamc_a1m0u1,
      RRamc_a1m1u0 = c / e,
      RRamc_a0m1u0 = 1 / RRamc_a1m1u0,
      RRamc_a1m0u0 = c_comp / e_comp,
      RRamc_a0m0u0 = 1 / RRamc_a1m0u0,

      # create possible RRum values
      RRumm_u1m1 = (b * n + d * n_comp) / (c * n + e * n_comp),
      RRumm_u0m1 = 1 / RRumm_u1m1,
      RRumm_u1m0 = (b_comp * n + d_comp * n_comp) / (c_comp * n + e_comp * n_comp),
      RRumm_u0m0 = 1 / RRumm_u1m0,

      # create possible RRum|a values
      RRumc_a1m1u1 = b / c,
      RRumc_a1m1u0 = 1 / RRumc_a1m1u1,
      RRumc_a1m0u1 = b_comp / c_comp,
      RRumc_a1m0u0 = 1 / RRumc_a1m0u1,
      RRumc_a0m1u1 = d / e,
      RRumc_a0m1u0 = 1 / RRumc_a0m1u1,
      RRumc_a0m0u1 = d_comp / e_comp,
      RRumc_a0m0u0 = 1 / RRumc_a0m0u1,

      # look at interaction direction
      AMdirection = case_when(
        c / e > 1 & b / d > 1 ~ "A-M positive",
        c / e < 1 & b / d < 1 ~ "A-M negative",
        c / e > 1 & b / d < 1 ~ "A-M positive U0, negative U1",
        c / e < 1 & b / d > 1 ~ "A-M negative U0, positive U1"
      ),
      UMdirection = case_when(
        d / e > 1 & b / c > 1 ~ "U-M positive",
        d / e < 1 & b / c < 1 ~ "U-M negative",
        d / e > 1 & b / c < 1 ~ "U-M positive A0, negative A1",
        d / e < 1 & b / c > 1 ~ "U-M negative A0, positive A1"
      ),
      MYdirection = case_when(
        h / i > 1 & f / g > 1 ~ "M-Y positive",
        h / i < 1 & f / g < 1 ~ "M-Y negative",
        h / i > 1 & f / g < 1 ~ "M-Y positive U0, negative U1",
        h / i < 1 & f / g > 1 ~ "M-Y negative U0, positive U1"
      ),
      UYdirection = case_when(
        g / i > 1 & f / h > 1 ~ "U-Y positive",
        g / i < 1 & f / h < 1 ~ "U-Y negative",
        g / i > 1 & f / h < 1 ~ "U-Y positive M0, negative M1",
        g / i < 1 & f / h > 1 ~ "U-Y negative M0, positive M1"
      ),
      M1int = case_when(
        d / e > b / c ~ "A0 > A1",
        d / e < b / c ~ "A0 < A1"
      ),
      M0int = case_when(
        d_comp / e_comp > b_comp / c_comp ~ "A0 > A1",
        d_comp / e_comp < b_comp / c_comp ~ "A0 < A1"
      ),
      Mint_directions = case_when(
        M1int == "A0 > A1" & M0int == "A0 < A1" ~
          "Not consistent",
        M1int == "A0 < A1" & M0int == "A0 > A1" ~
          "Not consistent",
        TRUE ~ "Consistent interaction"
      ),
      # summarize interaction directions
      qual_int = case_when(
        AMdirection %in% c("A-M positive", "A-M negative") &
          UMdirection %in% c("U-M positive", "U-M negative") &
          MYdirection %in% c("M-Y positive", "M-Y negative") &
          UYdirection %in% c("U-Y positive", "U-Y negative") ~ "no qual",
        TRUE ~ "qual"
      )
    ) %>%

    select(-(a:n)) %>%
    select(-(a_comp:n_comp)) %>%
    
    mutate(

      # choose RRs that match criteria
      RRa1u = pmap_dbl(select(., starts_with("RRau")), ~max(c(...))),
      RRa0u = pmap_dbl(select(., starts_with("RRau")), ~1 / min(c(...))),
      RRau = (obs_true > 1) * RRa1u + (obs_true < 1) * RRa0u,

      RRuy = pmap_dbl(select(., starts_with("RRuy")), ~max(c(...))),

      RRam = pmap_dbl(select(., starts_with("RRamm")), ~max(c(...))),
      RRam_no0 = pmap_dbl(select(., matches("RRamm_a\\dm1")), ~max(c(...))),

      RRam_strat = pmap_dbl(select(., starts_with("RRamc")), ~max(c(...))),
      RRam_strat_no0 = pmap_dbl(select(., matches("RRamc_a\\dm1u\\d")), ~max(c(...))),

      RRum = pmap_dbl(select(., starts_with("RRumm")), ~max(c(...))),
      RRum_no0 = pmap_dbl(select(., matches("RRumm_u\\dm1")), ~max(c(...))),

      RRum_strat = pmap_dbl(select(., starts_with("RRumc")), ~max(c(...))),
      RRum_strat_no0 = pmap_dbl(select(., matches("RRumc_a\\dm1u\\d")), ~max(c(...))),

      # Create bounding factors with new RRs
      BF_true = BF_func(RRau, RRuy),
      BF_am = BF_func(RRam, RRuy),
      BF_am_no0 = BF_func(RRam_no0, RRuy),
      BF_am_strat = BF_func(RRam_strat, RRuy),
      BF_am_strat_no0 = BF_func(RRam_strat_no0, RRuy),
      BF_um = BF_func(RRum, RRuy),
      BF_um_no0 = BF_func(RRum_no0, RRuy),
      BF_um_strat = BF_func(RRum_strat, RRuy),
      BF_um_strat_no0 = BF_func(RRum_strat_no0, RRuy)
    ) %>%

    select(-(RRau_m1u1:Mint_directions)) %>%

    mutate_at(
      # for all the bounding factors, create indicators if they succeed/fail
      # to bound bias, as well as correct NIE and NDE by each
      vars(starts_with("BF")),
      funs(
        fail = bias > .,
        GT = . > BF_true,
        NIE_correct = ifelse(obs_true > 1, NIE_obs * ., NIE_obs / .),
        NDE_correct = ifelse(obs_true < 1, NDE_obs * ., NDE_obs / .)
      )
    ) %>%

    select(-(RRa1u:BF_um_strat_no0)) %>%

    mutate(

      # create some variables for visualization/tables
      NIE_direction = ifelse(obs_true > 1, "NIE underestimated", "NIE overestimated"),
      NDE_direction = ifelse(obs_true < 1, "NDE underestimated", "NDE overestimated"),
      Bias = cut(bias, c(1, 1.05, 1.15, 1.3, 1.5, 2, Inf), right = F),
      qual_int = fct_recode(qual_int,
        "No interaction restriction" = "qual",
        "No qualitative interaction" = "no qual"
      )
    )
}

# function to get frequencies of failure
sum_data <- function(dat) {

  full_res <- dat %>%
    group_by(loglin) %>%
    summarise_at(vars(ends_with("GT"), ends_with("fail")), mean) %>%
    mutate(int = "any_int")

  no_int_res <- dat %>%
    filter(qual_int == "No qualitative interaction") %>%
    group_by(loglin) %>%
    summarise_at(vars(ends_with("GT"), ends_with("fail")), mean) %>%
    mutate(int = "no_int")

  rbind(full_res, no_int_res)
}

# function to make table 1
overall_res <- function(res) {
  GT <- res %>%
    arrange(loglin, int) %>%
    select(ends_with("GT")) %>%
    t()

  fail <- res %>%
    arrange(loglin, int) %>%
    select(ends_with("fail")) %>%
    t()

  tab <- cbind(GT, fail)
  tab[1, 1:4] <- NA

  addtorow <- list()
  addtorow$pos <- list(0, 0, 0)
  addtorow$command <- c(
    "& \\multicolumn{4}{c}{Alternative $\\text{RR }> \\RRau$} & \\multicolumn{4}{c}{Bias \\textgreater\\, bounding factor} \\\\ \\cmidrule(lr){2-5} \\cmidrule(lr){6-9} \n",
    "& \\multicolumn{2}{c}{$\\pr(M | A, U)$ uniform} & \\multicolumn{2}{c}{$\\pr(M | A, U)$ log-linear} & \\multicolumn{2}{c}{$\\pr(M | A, U)$ uniform} & \\multicolumn{2}{c}{$\\pr(M | A, U)$ log-linear}\\\\ \\cmidrule(lr){2-3} \\cmidrule(lr){4-5} \\cmidrule(lr){6-7} \\cmidrule(lr){8-9} \n",
    "Alternative & Full & Restricted & Full & Restricted & Full & Restricted & Full & Restricted \\\\\n"
  )
  rownames(tab) <- c(
    "$\\RRau$", "$\\text{RR}_{AM}$", "$\\text{RR}_{AM=1}$",
    "$\\text{RR}_{AM|U}$",
    "$\\text{RR}_{AM=1|U}$", "$\\text{RR}_{UM}$",
    "$\\text{RR}_{UM=1}$", "$\\text{RR}_{UM|A}$",
    "$\\text{RR}_{UM=1|A}$"
  )

  tab <- xtable(tab,
    digits = 6,
    align = c("r", rep("c", 8))
  )

  list(
    x = tab,
    comment = FALSE,
    booktabs = TRUE,
    include.rownames = TRUE,
    add.to.row = addtorow,
    include.colnames = FALSE,
    table.placement = "H",
    sanitize.text.function = function(x) x
  )
}

# function to make figure 2
ratio_true_plot <- function(dat, loglin_val, qualint_val, tag) {
  dat %>%
    gather(key = BF, value = corrected, BF_um_strat_no0_NIE_correct, BF_true_NIE_correct) %>%
    mutate(
      ratio = corrected / NIE_true,
      BF = fct_recode(BF,
        "True" = "BF_true_NIE_correct",
        "Alternative" = "BF_um_strat_no0_NIE_correct"
      )
    ) %>%
    filter(
      loglin == loglin_val,
      qual_int %in% qualint_val
    ) %>%
    ggplot() +
    geom_density(aes(ratio, col = BF)) +
    scale_x_log10(limits = c(.1, 10)) +
    scale_y_continuous(expand = c(0, 0)) +
    facet_grid(NIE_direction ~ Bias,
      labeller = labeller(.rows = label_value, .cols = label_both)
    ) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      strip.background.x = element_blank(),
      strip.background.y = element_rect(fill = "lightgrey", linetype = "blank")
    ) +
    labs(
      x = "Ratio of corrected NIE to true NIE (log scale)",
      subtitle = tag
    )
}

# function to make table 2
quant_rat <- function(dat) {

  p <- c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)

  newdat <- dat %>%
    gather(key = BF, value = corrected, BF_um_strat_no0_NIE_correct, BF_true_NIE_correct) %>%
    mutate(
      ratio = corrected / NIE_true,
      BF = fct_recode(BF,
        "True" = "BF_true_NIE_correct",
        "Alternative" = "BF_um_strat_no0_NIE_correct"
      )
    ) %>%
    mutate(ratiopos = case_when(obs_true > 1 ~ ratio, TRUE ~ 1 / ratio))

  ratio_vals <- newdat %>%
    nest(-BF) %>%
    mutate(x = map(.$data, ~tibble(
      q = paste0("q", p),
      val = quantile(.$ratiopos, p)
    ))) %>%
    select(-data) %>%
    unnest() %>%
    spread(q, val)

  ratio_vals_bias <- newdat %>%
    nest(-Bias, -BF) %>%
    mutate(x = map(.$data, ~tibble(
      q = paste0("q", p),
      val = quantile(.$ratiopos, p)
    ))) %>%
    select(-data) %>%
    unnest() %>%
    spread(q, val)

  tab <- cbind(t(ratio_vals[, -1]), t(ratio_vals_bias[, -(1:2)]))

  addtorow <- list()
  addtorow$pos <- list(0, 0)
  addtorow$command <- c(
    "& \\multicolumn{2}{c}{Overall} & \\multicolumn{2}{c}{[1, 1.05)} & \\multicolumn{2}{c}{[1.05, 1.15)} & \\multicolumn{2}{c}{[1.15, 1.3)} & \\multicolumn{2}{c}{[1.3, 1.5)}& \\multicolumn{2}{c}{[1.5, 2)}& \\multicolumn{2}{c}{[2, $\\infty$)}\\\\ \\cmidrule(lr){2-3} \\cmidrule(lr){4-5} \\cmidrule(lr){6-7} \\cmidrule(lr){8-9}\\cmidrule(lr){10-11} \\cmidrule(lr){12-13} \\cmidrule(lr){14-15} \n",
    "\\%ile & True & Alt. & True & Alt. & True & Alt. & True & Alt. & True & Alt. & True & Alt. & True & Alt.  \\\\\n"
  )
  rownames(tab) <- c("1st", "5th", "10th", "25th", "50th", "75th", "90th", "95th", "99th")

  tab <- xtable(tab,
    digits = 2,
    align = c("r", rep("c", ncol(tab)))
  )

  list(
    x = tab,
    comment = FALSE,
    booktabs = TRUE,
    include.rownames = TRUE,
    add.to.row = addtorow,
    include.colnames = FALSE,
    table.placement = "H",
    sanitize.text.function = function(x) x
  )
}

# function to make figure 2
resid_bias <- function(dat, loglin_val, qualint_val, tag) {

  dat %>%
    filter(
      BF_um_strat_no0_fail,
      loglin == loglin_val,
      qual_int %in% qualint_val
    ) %>%
    mutate(
      ratio = BF_um_strat_no0_NIE_correct / NIE_true,
      residual_bias = ratio * (obs_true < 1) + (1 / ratio) * (obs_true > 1),
      residual_quantile = cut(residual_bias, quantile(residual_bias, p = seq(0, 1, .1)), include.lowest = T)
    ) %>%
    ggplot() +
    geom_point(aes(bias, residual_bias, col = residual_quantile)) +
    scale_y_log10(breaks = c(1, 2, 5, 10), labels = c(1, 2, 5, 10), expand = c(0, 0)) +
    scale_x_log10(breaks = c(1, 2, 5, 10), labels = c(1, 2, 5, 10), expand = c(0, 0)) +
    geom_abline(slope = 1, intercept = 0) +
    scale_color_brewer(name = "Residual bias decile", palette = "Spectral") +
    theme_classic() +
    theme(legend.position = "bottom") +
    labs(
      x = "Bias (log scale)", y = "Residual bias (log scale)",
      subtitle = tag
    )
}

# function to make table 3
resid_bias_add <- function(dat){

  tab_any <- dat %>%
    filter(BF_um_strat_no0_fail) %>%
    mutate(
      diff = BF_um_strat_no0_NIE_correct - NIE_true,
      residual_bias = diff * (obs_true < 1) + (-diff) * (obs_true > 1)
    ) %>%
    nest(-loglin) %>%
    mutate(x = map(.$data, ~tibble(
      q = paste0("q", seq(.1, 1, .1)),
      val = quantile(.$residual_bias, seq(.1, 1, .1))
    ))) %>%
    select(-data) %>%
    unnest() %>%
    spread(q, val) %>%
    mutate_if(is.numeric, round, 3) %>%
    t()
  
  tab_no_int <- dat %>%
    filter(BF_um_strat_no0_fail, qual_int == "No qualitative interaction") %>%
    mutate(
      diff = BF_um_strat_no0_NIE_correct - NIE_true,
      residual_bias = diff * (obs_true < 1) + (-diff) * (obs_true > 1)
    ) %>%
    nest(-loglin) %>%
    mutate(x = map(.$data, ~tibble(
      q = paste0("q", seq(.1, 1, .1)),
      val = quantile(.$residual_bias, seq(.1, 1, .1))
    ))) %>%
    select(-data) %>%
    unnest() %>%
    spread(q, val) %>%
    mutate_if(is.numeric, round, 3) %>%
    t()

  tab <- cbind(tab_any[,1], tab_no_int[,1], tab_any[,2], tab_no_int[,2])[-1,]

  addtorow <- list()
  addtorow$pos <- list(0, 0)
  addtorow$command <- c(
    "& \\multicolumn{2}{c}{Uniform probabilities} & \\multicolumn{2}{c}{Log-linear model} \\\\ \\cmidrule(lr){2-3} \\cmidrule(lr){4-5} \n",
    "\\%ile & No restriction & No qualitative interaction & No restriction & No qualitative interaction  \\\\\n"
  )
  rownames(tab) <- c("10th", "20th", "30th", "40th", "50th", "60th", "70th", "80th", "90th", "Maximum")

  tab <- xtable(tab,
                digits = 2,
                align = c("r", rep("c", ncol(tab)))
  )

  list(
    x = tab,
    comment = FALSE,
    booktabs = TRUE,
    include.rownames = TRUE,
    add.to.row = addtorow,
    include.colnames = FALSE,
    table.placement = "H",
    sanitize.text.function = function(x) x
  )
}

# function to create other values for paper (how many had interaction)
vals <- function(dat) {
  dat %>%
    group_by(loglin) %>%
    summarise(
      qual_int = "No interaction restriction",
      n_noqual = n(),
      resUM = mean(BF_um_strat_no0_fail)
    ) %>% rbind(
  dat %>%
    filter(qual_int == "No qualitative interaction") %>%
    group_by(loglin) %>%
    summarise(
      qual_int = "No qualitative restriction",
      n_noqual = n(),
      resUM = mean(BF_um_strat_no0_fail)
    )
    ) %>% arrange(loglin, desc(qual_int))
}
