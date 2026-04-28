# analysis_strat_comparison.R
#
# Compares unstratified vs stratified Cox models for each endpoint (CAP, OXA,
# CAPOX). Presents the liver HU hazard ratio side-by-side to show what
# ignoring hospital does to the point estimate and its precision.
#
# Output: reports/strat_comparison_table.png (dot plot) and console table.
# Source of truth: presentation_2.Rmd pipeline.

library(survival)
library(dplyr)
library(ggplot2)
library(broom)

# ── Data & endpoint construction ─────────────────────────────────────────────

load("data/Dataset_liverfat.RData")
df <- dataset_leiden

make_endpoint <- function(data, endpoint) {
  if (endpoint == "CAP") {
    event <- data$TIMT_CAP2
    time  <- data$cycle_till_timt
    dose  <- data$starting_doseCAP
  }
  if (endpoint == "OXA") {
    event <- data$TIMT_OXA2
    time  <- data$cycle_till_timtOXA
    dose  <- data$starting_doseOXA
  }
  if (endpoint == "CAPOX") {
    event   <- data$TIMT_CAPOX2
    time    <- data$cycle_till_timtCAPOX
    doseCAP <- data$starting_doseCAP
    doseOXA <- data$starting_doseOXA
  }

  out <- data.frame(
    time, event,
    liverHU_manual = data$liverHU_manual,
    age      = data$age,
    sex      = data$sex,
    BBMI     = data$BBMI,
    BSMOKER  = data$BSMOKER,
    hospital = data$hospital
  )

  if (endpoint == "CAP" || endpoint == "OXA") out$starting_dose <- dose
  if (endpoint == "CAPOX") { out$doseCAP <- doseCAP; out$doseOXA <- doseOXA }

  out[complete.cases(out), ]
}

cap   <- make_endpoint(df, "CAP")
oxa   <- make_endpoint(df, "OXA")
capox <- make_endpoint(df, "CAPOX")

# ── Fit unstratified and stratified models ────────────────────────────────────
# Unstratified: hospital enters as a fixed-effect covariate.
# Stratified: hospital gets its own baseline hazard (presentation_2 approach).

fit_pair <- function(data, dose_terms, endpoint_label) {
  base_formula <- paste0(
    "Surv(time, event) ~ liverHU_manual + age + sex + BBMI + BSMOKER + ",
    dose_terms
  )

  unstrat <- coxph(as.formula(paste0(base_formula, " + hospital")),   data = data)
  strat   <- coxph(as.formula(paste0(base_formula, " + strata(hospital)")), data = data)

  bind_rows(
    tidy(unstrat, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term == "liverHU_manual") %>%
      mutate(Endpoint = endpoint_label, Model = "Unstratified"),
    tidy(strat, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term == "liverHU_manual") %>%
      mutate(Endpoint = endpoint_label, Model = "Stratified")
  )
}

results <- bind_rows(
  fit_pair(cap,   "starting_dose", "CAP"),
  fit_pair(oxa,   "starting_dose", "OXA"),
  fit_pair(capox, "doseCAP + doseOXA", "CAPOX")
)

# ── Clean comparison table ────────────────────────────────────────────────────

comparison_table <- results %>%
  mutate(
    HR    = round(estimate, 3),
    CI    = paste0(round(conf.low, 3), "–", round(conf.high, 3)),
    pval  = round(p.value, 3)
  ) %>%
  select(Endpoint, Model, HR, `95% CI` = CI, `p-value` = pval)

print(comparison_table)

# ── Dot plot ──────────────────────────────────────────────────────────────────

p <- results %>%
  mutate(
    Endpoint = factor(Endpoint, levels = c("CAP", "OXA", "CAPOX")),
    Model    = factor(Model, levels = c("Unstratified", "Stratified"))
  ) %>%
  ggplot(aes(x = estimate, y = Endpoint, color = Model, shape = Model)) +
  geom_point(size = 3, position = position_dodge(width = 0.4)) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 height = 0.2,
                 position = position_dodge(width = 0.4)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  scale_x_log10() +
  scale_color_manual(values = c("Unstratified" = "#D73027", "Stratified" = "#2166AC")) +
  labs(
    title    = "Liver HU hazard ratio: unstratified vs stratified by hospital",
    subtitle = "Stratified model allows hospital-specific baseline hazards",
    x        = "Hazard ratio (95% CI, log scale)",
    y        = NULL,
    color    = "Model",
    shape    = "Model"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave("reports/strat_comparison_plot.png", p, width = 7, height = 4, dpi = 150)
message("Saved: reports/strat_comparison_plot.png")
