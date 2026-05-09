# Canonical analysis for liver radiodensity & chemotherapy toxicity.
# Sourced by paper.Rmd, slides.Rmd, and analysis.Rmd — keep deterministic.
# Experimental code, sensitivity checks, and exploratory work belong in
# analysis.Rmd or in standalone scripts under R/, not here.

required_packages <- c("survival")
missing <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
if (length(missing) > 0) install.packages(missing)

library(survival)

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

cox_cap <- coxph(
  Surv(time, event) ~ liverHU_manual + age + sex + BBMI +
    BSMOKER + starting_dose + strata(hospital),
  data = cap
)
cox_oxa <- coxph(
  Surv(time, event) ~ liverHU_manual + age + sex + BBMI +
    BSMOKER + starting_dose + strata(hospital),
  data = oxa
)
cox_capox <- coxph(
  Surv(time, event) ~ liverHU_manual + age + sex + BBMI +
    BSMOKER + doseCAP + doseOXA + strata(hospital),
  data = capox
)

ph_cap   <- cox.zph(cox_cap)
ph_oxa   <- cox.zph(cox_oxa)
ph_capox <- cox.zph(cox_capox)
