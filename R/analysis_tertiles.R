# analysis_tertiles.R
#
# Tertile analysis of liver radiodensity (HU) and time to chemotherapy toxicity.
# Splits liverHU_manual into three equal-sized groups (low / medium / high),
# refits stratified Cox models for CAP, OXA, and CAPOX using tertiles instead
# of continuous HU, and saves one forest plot per endpoint to reports/.
#
# Source of truth: presentation_2.Rmd pipeline (make_endpoint function below
# mirrors the fixed version in that file).

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

# ── Tertile assignment ────────────────────────────────────────────────────────
# Tertile breaks computed per-endpoint so each group has equal N within that
# risk set (patients differ across endpoints due to complete-case filtering).

add_tertile <- function(data) {
  breaks <- quantile(data$liverHU_manual, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
  data$hu_tertile <- cut(data$liverHU_manual, breaks = breaks,
                         labels = c("Low HU", "Medium HU", "High HU"),
                         include.lowest = TRUE)
  data$hu_tertile <- relevel(data$hu_tertile, ref = "Low HU")
  data
}

cap   <- add_tertile(cap)
oxa   <- add_tertile(oxa)
capox <- add_tertile(capox)

# ── Tertile Cox models ────────────────────────────────────────────────────────

cox_cap_tert <- coxph(
  Surv(time, event) ~ hu_tertile + age + sex + BBMI +
    BSMOKER + starting_dose + strata(hospital),
  data = cap
)

cox_oxa_tert <- coxph(
  Surv(time, event) ~ hu_tertile + age + sex + BBMI +
    BSMOKER + starting_dose + strata(hospital),
  data = oxa
)

cox_capox_tert <- coxph(
  Surv(time, event) ~ hu_tertile + age + sex + BBMI +
    BSMOKER + doseCAP + doseOXA + strata(hospital),
  data = capox
)

# ── Forest plot helper ────────────────────────────────────────────────────────

plot_tertile_forest <- function(model, endpoint_label) {
  tidy(model, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(grepl("hu_tertile", term)) %>%
    mutate(
      label = gsub("hu_tertile", "", term),
      label = factor(label, levels = c("Medium HU", "High HU"))
    ) %>%
    ggplot(aes(x = estimate, y = label)) +
    geom_point(size = 3, color = "#2166AC") +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                   height = 0.2, color = "#2166AC") +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
    scale_x_log10() +
    labs(
      title    = paste0("Liver HU tertiles and toxicity risk: ", endpoint_label),
      subtitle = "Reference = Low HU tertile | Stratified Cox model",
      x        = "Hazard ratio (95% CI, log scale)",
      y        = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank())
}

# ── Save plots ────────────────────────────────────────────────────────────────

p_cap   <- plot_tertile_forest(cox_cap_tert,   "CAP")
p_oxa   <- plot_tertile_forest(cox_oxa_tert,   "OXA")
p_capox <- plot_tertile_forest(cox_capox_tert, "CAPOX")

ggsave("reports/tertile_forest_CAP.png",   p_cap,   width = 7, height = 4, dpi = 150)
ggsave("reports/tertile_forest_OXA.png",   p_oxa,   width = 7, height = 4, dpi = 150)
ggsave("reports/tertile_forest_CAPOX.png", p_capox, width = 7, height = 4, dpi = 150)

message("Saved: reports/tertile_forest_CAP.png, _OXA.png, _CAPOX.png")

# ── Print summary tables ──────────────────────────────────────────────────────

tert_results <- bind_rows(
  tidy(cox_cap_tert,   exponentiate = TRUE, conf.int = TRUE) %>%
    filter(grepl("hu_tertile", term)) %>% mutate(Endpoint = "CAP"),
  tidy(cox_oxa_tert,   exponentiate = TRUE, conf.int = TRUE) %>%
    filter(grepl("hu_tertile", term)) %>% mutate(Endpoint = "OXA"),
  tidy(cox_capox_tert, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(grepl("hu_tertile", term)) %>% mutate(Endpoint = "CAPOX")
) %>%
  mutate(
    Tertile = gsub("hu_tertile", "", term),
    HR      = round(estimate, 3),
    CI_low  = round(conf.low, 3),
    CI_high = round(conf.high, 3),
    p       = round(p.value, 3)
  ) %>%
  select(Endpoint, Tertile, HR, CI_low, CI_high, p)

print(tert_results)
