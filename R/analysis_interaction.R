# analysis_interaction.R
#
# Likelihood ratio test for a liver HU x hospital interaction in stratified
# Cox models. Tests whether the association between liver HU and toxicity
# differs across hospitals.
#
# Statistical validity note:
#   With 11 hospitals the interaction term adds 10 extra parameters (one slope
#   per hospital minus the reference). Rule of thumb for Cox models is ~10
#   events per parameter. Event counts are CAP=138, OXA=189, CAPOX=182, giving
#   roughly 14/19/18 events per interaction parameter — borderline but
#   acceptable for an exploratory test. Interpret p-values cautiously; the test
#   is underpowered to detect small interactions.
#
# Output: console table + reports/interaction_lrt.png
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

# ── LRT helper ───────────────────────────────────────────────────────────────
# Fits main-effects model and model with liverHU_manual:hospital interaction,
# both stratified by hospital. Returns a one-row tibble with LRT results.
#
# Note: in a stratified model strata(hospital) absorbs baseline hazard
# differences; the interaction term liverHU_manual:hospital tests whether the
# *slope* of HU differs across hospitals, which is the scientifically
# meaningful question.

run_lrt <- function(data, dose_terms, endpoint_label) {
  base <- paste0(
    "Surv(time, event) ~ liverHU_manual + age + sex + BBMI + BSMOKER + ",
    dose_terms, " + strata(hospital)"
  )
  interaction_f <- paste0(
    "Surv(time, event) ~ liverHU_manual * hospital + age + sex + BBMI + BSMOKER + ",
    dose_terms, " + strata(hospital)"
  )

  m_main <- coxph(as.formula(base),        data = data)
  m_int  <- coxph(as.formula(interaction_f), data = data)

  lrt     <- anova(m_main, m_int)
  chisq   <- lrt[["Chisq"]][2]
  df_diff <- lrt[["Df"]][2]
  pval    <- lrt[["Pr(>|Chi|)"]][2]

  n_hospitals <- length(unique(data$hospital))

  tibble(
    Endpoint    = endpoint_label,
    N           = nrow(data),
    Events      = sum(data$event),
    Hospitals   = n_hospitals,
    Extra_params = df_diff,
    Chisq       = round(chisq, 3),
    df          = df_diff,
    p_value     = round(pval, 4)
  )
}

lrt_results <- bind_rows(
  run_lrt(cap,   "starting_dose",    "CAP"),
  run_lrt(oxa,   "starting_dose",    "OXA"),
  run_lrt(capox, "doseCAP + doseOXA", "CAPOX")
)

print(lrt_results)

# ── Plot ──────────────────────────────────────────────────────────────────────

p <- lrt_results %>%
  mutate(
    Endpoint  = factor(Endpoint, levels = c("CAP", "OXA", "CAPOX")),
    label     = paste0("χ²=", Chisq, "\np=", p_value),
    sig_color = ifelse(p_value < 0.05, "#D73027", "#4D9221")
  ) %>%
  ggplot(aes(x = Endpoint, y = Chisq, fill = sig_color)) +
  geom_col(width = 0.5) +
  geom_text(aes(label = label), vjust = -0.4, size = 3.5) +
  scale_fill_identity() +
  labs(
    title    = "LRT: liver HU × hospital interaction",
    subtitle = "Compares stratified Cox model with and without HU × hospital term",
    x        = NULL,
    y        = "Chi-square statistic"
  ) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank())

ggsave("reports/interaction_lrt.png", p, width = 6, height = 4, dpi = 150)
message("Saved: reports/interaction_lrt.png")
