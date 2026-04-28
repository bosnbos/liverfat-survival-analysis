# analysis_hospital_events.R
#
# Visualises hospital-level toxicity event rates with 95% Wilson confidence
# intervals for each endpoint (CAP, OXA, CAPOX). One plot per endpoint saved
# to reports/. Demonstrates the hospital-level variation that motivates
# stratification in the Cox models.
#
# Source of truth: presentation_2.Rmd pipeline.

library(survival)
library(dplyr)
library(ggplot2)

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

# ── Hospital event rate helper ────────────────────────────────────────────────
# Wilson score interval for a proportion: more reliable than Wald CI when
# cell counts are small, as they are for some hospitals here.

wilson_ci <- function(x, n, z = 1.96) {
  p     <- x / n
  denom <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / denom
  half   <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denom
  list(low = pmax(0, centre - half), high = pmin(1, centre + half))
}

hospital_rates <- function(data, endpoint_label) {
  data %>%
    group_by(hospital) %>%
    summarise(n = n(), events = sum(event), .groups = "drop") %>%
    mutate(
      rate   = events / n,
      ci     = wilson_ci(events, n),
      ci_low = ci$low,
      ci_high = ci$high,
      Endpoint = endpoint_label
    ) %>%
    select(-ci)
}

rates_all <- bind_rows(
  hospital_rates(cap,   "CAP"),
  hospital_rates(oxa,   "OXA"),
  hospital_rates(capox, "CAPOX")
)

# ── Plot helper ───────────────────────────────────────────────────────────────

plot_hospital_events <- function(rates, endpoint_label) {
  rates %>%
    filter(Endpoint == endpoint_label) %>%
    arrange(rate) %>%
    mutate(hospital = factor(hospital, levels = hospital)) %>%
    ggplot(aes(x = rate, y = hospital)) +
    geom_point(size = 3, color = "#2166AC") +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high),
                   height = 0.3, color = "#2166AC") +
    geom_text(aes(label = paste0(events, "/", n)),
              hjust = -0.2, size = 3, color = "grey30") +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                       limits = c(0, 1.05)) +
    labs(
      title    = paste0("Toxicity event rate per hospital: ", endpoint_label),
      subtitle = "95% Wilson confidence intervals | sorted by event rate",
      x        = "Event rate",
      y        = "Hospital"
    ) +
    theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank())
}

p_cap   <- plot_hospital_events(rates_all, "CAP")
p_oxa   <- plot_hospital_events(rates_all, "OXA")
p_capox <- plot_hospital_events(rates_all, "CAPOX")

ggsave("reports/hospital_events_CAP.png",   p_cap,   width = 7, height = 5, dpi = 150)
ggsave("reports/hospital_events_OXA.png",   p_oxa,   width = 7, height = 5, dpi = 150)
ggsave("reports/hospital_events_CAPOX.png", p_capox, width = 7, height = 5, dpi = 150)

message("Saved: reports/hospital_events_CAP.png, _OXA.png, _CAPOX.png")

# ── Print summary ─────────────────────────────────────────────────────────────

print(rates_all %>% select(Endpoint, hospital, n, events, rate, ci_low, ci_high) %>%
        mutate(across(c(rate, ci_low, ci_high), round, 3)))
