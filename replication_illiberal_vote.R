################################################################################
## Replication script:
## Attitudinal predictors of illiberal vote choice in the 2023 Slovak election
## Author:       Jakub Szabó
## Last updated: 2025-05-29  (patched 2026-07-22)
## R version:    >= 4.2.0 (developed and tested under 4.x)
##
## ----------------------------------------------------------------------------
## OVERVIEW
## ----------------------------------------------------------------------------
## Reproduces the models, tables, and figures in the manuscript and appendices.
## Main text:  Figure 1 (odds ratios), Figure 2 (multinomial AMEs).
## Appendix A: Figure A1 (CHES party categorisation), Figure A2 (populism
##             distributions), Table A1 (descriptive statistics).
## Appendix B: Table B1 / Figure B1 (populism EFA), Table B2 (trust EFA).
## Appendix C: Tables C1-C4 and Figures C1-C4 (binary-model robustness).
## Appendix D: Table D1 and Figures D1-D4 (multinomial model and robustness).
## ----------------------------------------------------------------------------
## DIRECTORY STRUCTURE (expected)
## ----------------------------------------------------------------------------
##   <project_root>/
##     replication_illiberal_vote.R   <- this file
##     data/
##       Datafile_ISSP.sav            <- survey data (not redistributed)
##       CHES_2024.dta                <- CHES 2024 expert survey (Appendix A1)
##     output/                        <- created automatically; tables + figures
##
## Paths are resolved with `here`; no machine-specific paths are required.
##
## ----------------------------------------------------------------------------
## REPRODUCIBILITY NOTES
## ----------------------------------------------------------------------------
## Survey weights (`weight_FINAL`) are applied to all estimated models.
##
## Most estimators here are deterministic, but two are NOT:
##   - psych::fa.parallel (Appendix Figures B1) is simulation-based;
##   - nnet::multinom uses random starting weights.
## Both are given a local `set.seed()` immediately before they run, and the
## multinomial models are checked for convergence. For a fully frozen
## dependency set, wrap this project in `renv` (renv::init(); renv::snapshot());
## `sessionInfo()` is written to output/ at the end of this script.
################################################################################


## ============================================================================
## 0. SETUP
## ============================================================================

rm(list = ls())
set.seed(2023)   # global seed; local seeds are also set before stochastic steps

required_packages <- c(
  "here",            # project-relative paths
  "haven",           # read SPSS (.sav) and Stata (.dta) files
  "tidyverse",       # data wrangling and plotting
  "psych",           # alpha, fa, fa.parallel, principal, polychoric, KMO
  "nnet",            # multinomial logit
  "sandwich",        # robust (HC0) variance-covariance
  "marginaleffects", # average marginal effects
  "emmeans",         # estimated marginal means / predicted probabilities  ## PATCH: de-duplicated comment
  "ggeffects",       # predicted probabilities
  "DescTools",       # pseudo-R2 and information criteria
  "stargazer",       # regression tables
  "patchwork"        # figure composition
)
install.packages(setdiff(required_packages, rownames(installed.packages())))
invisible(lapply(required_packages, library, character.only = TRUE))

data_path  <- here::here("data", "Datafile_ISSP.sav")
ches_path  <- here::here("data", "CHES_2024.dta")
output_dir <- here::here("output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


## ----------------------------------------------------------------------------
## 0.1 Small helpers used throughout
## ----------------------------------------------------------------------------

## Anchor a PCA/factor score to a stable sign so it is comparable across
## sub-samples (principal()'s sign convention is otherwise arbitrary): flip the
## score if it correlates negatively with the mean of its input items.
anchor_sign <- function(scores, reference) {
  if (isTRUE(cor(scores, reference, use = "complete.obs") < 0)) -scores else scores
}

## Robust (HC0) point estimates, SEs, p-values, and 95% CIs for a fitted glm.
## `or`/`ci` are on the odds-ratio scale; `se` is kept on the log-odds scale
## (that is the scale on which p-values are computed).
robust_stats <- function(model) {
  se <- sqrt(diag(vcovHC(model, type = "HC0")))
  b  <- coef(model)
  list(
    or = exp(b),
    se = se,
    p  = 2 * pnorm(abs(b / se), lower.tail = FALSE),
    ci = cbind(exp(b - 1.96 * se), exp(b + 1.96 * se))
  )
}


## ============================================================================
## 1. DATA IMPORT
## ============================================================================
data <- read_sav(data_path, encoding = "latin1")

## ============================================================================
## 2. SAMPLE RESTRICTION
## ============================================================================
## Drop respondents below voting age at the 2023 election (q10lh_a == 3).

data <- data %>%
  filter(q10lh_a != 3) %>%
  mutate(vote_particip = q10lh_a)


## ============================================================================
## 3. DEPENDENT VARIABLES
## ============================================================================
## 3.1 Recode reported 2023 lower-house vote into named parties.

data <- data %>%
  mutate(vote_outcome = case_when(
    q10lh_b %in% c(6, 7, 9, 10, 14, 20, 21, 96) ~ "other",
    q10lh_a == 1 & is.na(q10lh_b)               ~ "other",
    q10lh_a %in% c(-8, -9)                       ~ "other",
    q10lh_b == 3  ~ "PS",
    q10lh_b == 5  ~ "OLANO",
    q10lh_b == 12 ~ "SAS",
    q10lh_b == 13 ~ "Sme Rodina",
    q10lh_b == 15 ~ "SNS",
    q10lh_b == 16 ~ "Smer",
    q10lh_b == 17 ~ "Hlas",
    q10lh_b == 18 ~ "Aliancia",
    q10lh_b == 22 ~ "Demokrati",
    q10lh_b == 23 ~ "KDH",
    q10lh_b == 25 ~ "Republika",
    q10lh_a == 2  ~ "non-voter",
    TRUE          ~ NA_character_
  ))

## vote so the decision below is explicit rather than silent.
message("Unclassifiable vote_outcome (NA) rows: ", sum(is.na(data$vote_outcome)))

## 3.2 Binary outcomes.
##     vote_illiberal      : Smer-SD, SNS, Republika  (main definition)
##     vote_illiberal_hlas : adds Hlas                (Appendix Table C3)
data <- data %>%
  mutate(
    vote_illiberal = case_when(
      is.na(vote_outcome)                                ~ NA_real_,
      vote_outcome %in% c("Smer", "Republika", "SNS")    ~ 1,
      TRUE                                               ~ 0
    ),
    vote_illiberal_hlas = case_when(
      is.na(vote_outcome)                                        ~ NA_real_,
      vote_outcome %in% c("Smer", "Republika", "SNS", "Hlas")    ~ 1,
      TRUE                                                       ~ 0
    )
  )


## ============================================================================
## 4. INDEPENDENT VARIABLES AND ANALYTIC SAMPLES
## ============================================================================

## ----------------------------------------------------------------------------
## 4.1 Socio-demographic controls and the base sample (`df_base`)
## ----------------------------------------------------------------------------
## Age is computed relative to 2024 (fielding year)

df_base <- data %>%
  filter(b1 %in% 1:2) %>%
  mutate(female = factor(if_else(b1 == 1, 0, 1), levels = c(0, 1))) %>%
  filter(between(b2, 1900, 2010)) %>%
  mutate(age = (2024 - b2) / 10) %>%
  filter(b3 > 0) %>%
  mutate(edu = as.numeric(b3)) %>%
  filter(between(qq19, 0, 10)) %>% mutate(left_right = as.numeric(qq19)) %>%
  filter(between(q01, 1, 4))   %>% mutate(interest = 5 - q01) %>%
  filter(between(qq22, 1, 4))  %>% mutate(satisfaction_democracy = 5 - qq22)

## ----------------------------------------------------------------------------
## 4.2 Political trust (latent index; diagnostics for Appendix Table B2)
## ----------------------------------------------------------------------------

data %>%                                   # reliability of the seven trust items
  select(q07a, q07b, q07c, q07d, q07e, q07f, q07g) %>%
  filter(if_all(everything(), ~ .x %in% 1:4)) %>%
  psych::alpha()

efa_trust <- data %>%
  select(q07a, q07b, q07c, q07d, q07e, q07f, q07g) %>%
  filter(if_all(everything(), ~ .x %in% 1:4))

poly_trust <- polychoric(efa_trust)
set.seed(2023)                         
fa.parallel(poly_trust$rho, fa = "fa", n.obs = nrow(efa_trust), fm = "ml")
fa_trust <- fa(poly_trust$rho, nfactors = 3, rotate = "oblimin", fm = "ml")
print(fa_trust)                            # Appendix Table B2

## Trust in Parliament/Government/Justice/Parties (q07a/b/c/e) -> single index
df_base <- df_base %>%
  filter(if_all(c(q07a, q07b, q07c, q07e), ~ .x %in% 1:4)) %>%
  mutate(across(c(q07a, q07b, q07c, q07e), ~ 5 - .x))      # higher = more trust

pca_trust <- principal(df_base[, c("q07a", "q07b", "q07c", "q07e")],
                       nfactors = 1, scores = TRUE)

## principal()'s arbitrary sign convention.
trust_ref <- rowMeans(df_base[, c("q07a", "q07b", "q07c", "q07e")], na.rm = TRUE)
df_base <- df_base %>% mutate(trust = anchor_sign(pca_trust$scores[, 1], trust_ref))

## ----------------------------------------------------------------------------
## 4.3 Illiberal attitudes (`df_illib`)
## ----------------------------------------------------------------------------
df_illib <- df_base %>%
  filter(between(q04b, 1, 5)) %>% mutate(government_constraints = q04b) %>%
  filter(between(q19,  1, 7)) %>% mutate(political_rights  = 8 - q19) %>%
  filter(between(q20,  1, 7)) %>% mutate(political_rights2 = 8 - q20) %>%   # validation item
  filter(between(q7_5, 1, 5)) %>% mutate(state_neutral  = 6 - q7_5) %>%
  filter(between(q6_2, 1, 5)) %>% mutate(closed_society = q6_2) %>%
  filter(between(q04a, 1, 5)) %>% mutate(democracy = q04a)

## ----------------------------------------------------------------------------
## 4.4 Authoritarianism (`df_auth`)
## ----------------------------------------------------------------------------
add_authoritarianism <- function(df) {
  df %>%
    filter(between(q04c, 1, 5)) %>%
    mutate(authoritarianism = 6 - q04c)
}

df_auth <- df_base %>% add_authoritarianism()

## ----------------------------------------------------------------------------
## 4.5 Populist attitudes (`df_pop`; diagnostics for Appendix Table B1/Figure B1)
## ----------------------------------------------------------------------------
## Helper that, given a data frame with reverse-coded q28-q31/q05c, adds every
## populism operationalisation used in the paper (keeps df_pop and df_full in
## sync without duplicating the logic).

add_populism <- function(df) {
  pca1 <- principal(df[, c("q28", "q29", "q30")],  nfactors = 1, scores = TRUE)
  pca2 <- principal(df[, c("q05c", "q31", "q30")], nfactors = 1, scores = TRUE)
  ref1 <- rowMeans(df[, c("q28", "q29", "q30")],  na.rm = TRUE)
  ref2 <- rowMeans(df[, c("q05c", "q31", "q30")], na.rm = TRUE)
  df %>%
    mutate(
      people_center  = q28,
      anti_elit_pol  = q29,
      anti_elit_bus  = q31,
      manichean      = q30,
      pca_populism1  = anchor_sign(pca1$scores[, 1], ref1),
      pca_populism2  = anchor_sign(pca2$scores[, 1], ref2),
      populism_additive = rowMeans(select(., q28, q29, q30), na.rm = TRUE),
      populism_dummy = if_else(q28 %in% c(4, 5) &
                                 q29 %in% c(4, 5) &
                                 q30 %in% c(4, 5), 1, 0)
    )
}

efa_pop <- data %>%
  select(q28, q29, q30, q31) %>%
  filter(if_all(everything(), ~ .x %in% 1:5))

KMO(efa_pop)
cortest.bartlett(efa_pop)
poly_pop <- polychoric(efa_pop)$rho
set.seed(2023)                                                         ## PATCH: seed parallel analysis
fa.parallel(poly_pop, fa = "fa", n.obs = nrow(efa_pop))                 # Figure B1
print(fa(poly_pop, nfactors = 1, rotate = "oblimin", fm = "ml"))       # Table B1

df_pop <- df_base %>%
  filter(if_all(c(q28, q29, q30, q31, q05c), ~ .x %in% 1:5)) %>%
  mutate(across(c(q28, q29, q30, q31, q05c), ~ 6 - .x)) %>%             # higher = more populist
  add_populism()

## ----------------------------------------------------------------------------
## 4.6 Nativist attitudes (`df_nat`)
## ----------------------------------------------------------------------------
add_nativism <- function(df) {
  df %>%
    filter(if_all(c(q7_1, q7_2, q7_3, q7_4), ~ .x %in% 1:5)) %>%
    mutate(across(c(q7_1, q7_2, q7_3, q7_4), as.numeric),
           q7_1 = 6 - q7_1, q7_3 = 6 - q7_3) %>%
    mutate(anti_immigration = rowMeans(select(., q7_1, q7_2, q7_3, q7_4), na.rm = TRUE))
}

data %>%                                   # reliability of the four-item battery
  select(q7_1, q7_2, q7_3, q7_4) %>%
  filter(if_all(everything(), ~ .x %in% 1:5)) %>%
  mutate(q7_1 = 6 - q7_1, q7_3 = 6 - q7_3) %>%
  na.omit() %>%
  psych::alpha()

df_nat <- df_base %>% add_nativism()

## ----------------------------------------------------------------------------
## 4.7 Full sample (`df_full`) and national-purity variant (`df_full_np`)
## ----------------------------------------------------------------------------
## df_full = illiberal-attitudes sample augmented with authoritarianism,
## populism, and nativism (listwise-complete on all attitudinal predictors).

df_full <- df_illib %>%
  add_authoritarianism() %>%
  filter(if_all(c(q28, q29, q30, q31, q05c), ~ .x %in% 1:5)) %>%
  mutate(across(c(q28, q29, q30, q31, q05c), ~ 6 - .x)) %>%
  add_populism() %>%
  add_nativism()

## Adds the national-purity item (q38) for the nativism robustness checks.
df_full_np <- df_full %>%
  filter(between(q38, 1, 5)) %>%
  mutate(national_purity = 6 - q38)


## ============================================================================
## 5. DESCRIPTIVE OUTPUTS
## ============================================================================

## ----------------------------------------------------------------------------
## 5.1 CHES party categorisation (Appendix Figure A1)
## ----------------------------------------------------------------------------
## NB: the CHES file is read into its OWN object (`ches`) so that the survey
## `data` object is not overwritten. (Verify country == 28 is Slovakia in your
## CHES release.)

ches <- read_dta(ches_path) %>%
  filter(country == 28) %>%
  select(party, eu_position, civlib_laworder, galtan, nationalism,
         executive_power, judicial_independence) %>%
  mutate(across(c(eu_position, civlib_laworder, galtan, nationalism,
                  executive_power, judicial_independence), as.numeric)) %>%
  mutate(party = if_else(party == "Smer-SD", "Smer", party))

ches_bar <- function(var, ylab) {
  ggplot(ches, aes(x = reorder(party, .data[[var]], decreasing = TRUE),
                   y = .data[[var]], fill = party)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_grey() +
    theme_bw() +
    theme(legend.position = "none") +
    labs(x = "", y = ylab)
}

plot_ideology <-
  (ches_bar("eu_position",           "Position on the EU\n(7 - Strongly in Favor)") |
   ches_bar("civlib_laworder",       "Liberty vs. Order\n(10 - Authoritarian)")) /
  (ches_bar("galtan",                "Democratic Freedoms\n(10 - TAN)") |
   ches_bar("nationalism",           "Cosmopol. vs. National.\n(10 - National.)")) /
  (ches_bar("executive_power",       "Executive Power\n(10 - Unconstrained Leaders)") |
   ches_bar("judicial_independence", "Judicial Independence\n(10 - Controlled Judiciary)"))

ggsave(file.path(output_dir, "fig_A1_ches_categorisation.png"),
       plot_ideology, width = 11, height = 12, dpi = 300)

## ----------------------------------------------------------------------------
## 5.2 Populism distributions (Appendix Figure A2)
## ----------------------------------------------------------------------------

p_pca <- df_pop %>%
  select(pca_populism1, pca_populism2) %>%
  pivot_longer(everything(), names_to = "dimension", values_to = "value") %>%
  ggplot(aes(x = value)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  facet_wrap(~ dimension, scales = "free") +
  theme_minimal() +
  labs(title = "PCA-based populism measures", x = "PCA score", y = "Density")

p_items <- df_pop %>%
  select(people_center, anti_elit_pol, anti_elit_bus, manichean) %>%
  pivot_longer(everything(), names_to = "dimension", values_to = "value") %>%
  ggplot(aes(x = value)) +
  geom_histogram(binwidth = 1, boundary = 0.5, fill = "grey70", color = "black") +
  facet_wrap(~ dimension) +
  scale_x_continuous(breaks = 1:5) +
  theme_minimal() +
  labs(title = "Individual populism dimensions", x = "Response", y = "Count")

ggsave(file.path(output_dir, "fig_A2_populism_distributions.png"),
       p_pca / p_items, width = 12, height = 8, dpi = 300)

## ----------------------------------------------------------------------------
## 5.3 Descriptive statistics (Appendix Table A1)
## ----------------------------------------------------------------------------
## Computed on the full analytic sample (df_full), unweighted. Skewness and
## kurtosis are moment-based; kurtosis is in EXCESS form (normal = 0).

.skewness <- function(x) {
  x <- x[!is.na(x)]; n <- length(x); m <- mean(x)
  (sum((x - m)^3) / n) / (sum((x - m)^2) / n)^(3 / 2)
}
.kurtosis_excess <- function(x) {
  x <- x[!is.na(x)]; n <- length(x); m <- mean(x)
  (sum((x - m)^4) / n) / (sum((x - m)^2) / n)^2 - 3
}

descriptive_vars <- c(
  government_constraints = "Unconstrained Government",
  political_rights       = "Political Rights Curtailment",
  state_neutral          = "State Neutral Opposition",
  closed_society         = "Closed Society",
  democracy              = "Anti-democracy",
  authoritarianism       = "Authoritarianism",
  pca_populism1          = "Populism",
  anti_immigration       = "Nativism",
  female                 = "Female",
  age                    = "Age (in decades)",
  edu                    = "Education",
  left_right             = "Left-right Ideology",
  interest               = "Interest in Politics",
  satisfaction_democracy = "Satisfaction with Democracy",
  trust                  = "Trust"
)

descriptive_table <- df_full %>%
  mutate(female = as.numeric(as.character(female))) %>%
  select(all_of(names(descriptive_vars))) %>%
  summarise(across(everything(), list(
    Mean = ~ mean(.x, na.rm = TRUE), Median = ~ median(.x, na.rm = TRUE),
    SD = ~ sd(.x, na.rm = TRUE), Min = ~ min(.x, na.rm = TRUE),
    Max = ~ max(.x, na.rm = TRUE), Skewness = ~ .skewness(.x),
    Kurtosis = ~ .kurtosis_excess(.x)
  ), .names = "{.col}__{.fn}")) %>%
  pivot_longer(everything(), names_to = c("variable", "stat"), names_sep = "__") %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  mutate(Variable = factor(descriptive_vars[variable], levels = descriptive_vars)) %>%
  arrange(Variable) %>%
  select(Variable, Mean, Median, SD, Min, Max, Skewness, Kurtosis) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

print(descriptive_table)
write.csv(descriptive_table,
          file.path(output_dir, "table_A1_descriptives.csv"), row.names = FALSE)
stargazer(as.data.frame(descriptive_table), type = "html", summary = FALSE,
          rownames = FALSE, title = "Table A1. Descriptive Statistics of the Data",
          notes = "Source: Own elaboration based on CSES-ISSP data.",
          out = file.path(output_dir, "table_A1_descriptives.html"))


## ============================================================================
## 6. BINARY LOGISTIC MODELS (main specification)
## ============================================================================
ml1 <- glm(vote_illiberal ~ government_constraints + political_rights +
             state_neutral + closed_society + democracy +
             female + age + edu + left_right + interest + satisfaction_democracy + trust,
           weights = weight_FINAL, data = df_illib, family = binomial)

ml2 <- glm(vote_illiberal ~ authoritarianism +
             female + age + edu + left_right + interest + satisfaction_democracy + trust,
           weights = weight_FINAL, data = df_auth, family = binomial)

ml3 <- glm(vote_illiberal ~ pca_populism1 +
             female + age + edu + left_right + interest + satisfaction_democracy + trust,
           weights = weight_FINAL, data = df_pop, family = binomial)

ml4 <- glm(vote_illiberal ~ anti_immigration +
             female + age + edu + left_right + interest + satisfaction_democracy + trust,
           weights = weight_FINAL, data = df_nat, family = binomial)

ml5 <- glm(vote_illiberal ~ government_constraints + political_rights +
             state_neutral + closed_society + democracy + authoritarianism +
             pca_populism1 + anti_immigration +
             female + age + edu + left_right + interest + satisfaction_democracy + trust,
           weights = weight_FINAL, data = df_full, family = binomial)


## ============================================================================
## 7. REGRESSION TABLES AND BINARY-MODEL ROBUSTNESS (Appendix C)
## ============================================================================

covariate_labels_main <- c(
  "Unconstrained Government", "Political Rights Curtailment",
  "State Neutral Opposition", "Closed Society", "Anti-Democracy",
  "Authoritarianism", "Populism", "Nativism",
  "Female", "Age", "Education", "Left-right Ideology",
  "Interest in Politics", "Satisfaction with Democracy", "Trust"
)

or_table <- function(models, out, dep_label,
                     covariate.labels = covariate_labels_main,
                     notes = NULL) {
  rs <- lapply(models, robust_stats)
  stargazer(
    models,
    type             = "html",
    dep.var.labels   = dep_label,
    covariate.labels = covariate.labels,
    coef             = lapply(rs, `[[`, "or"),
    ci               = TRUE,
    ci.custom        = lapply(rs, `[[`, "ci"),
    p                = lapply(rs, `[[`, "p"),
    add.lines = list(
      c("n",                    vapply(models, nobs, numeric(1))),
      c("Pseudo R2 (McFadden)", vapply(models, function(m) round(PseudoR2(m, which = "McFadden"), 3), numeric(1))),
      c("AIC",                  vapply(models, function(m) round(AIC(m), 2), numeric(1))),
      c("BIC",                  vapply(models, function(m) round(BIC(m), 2), numeric(1)))
    ),
    omit.stat = c("aic", "n"),
    notes     = notes,
    out       = out
  )
}

## Note appended to the cross-sample tables (C1-C3): the five columns are fit on
## DIFFERENT analytic samples (df_illib/df_auth/df_pop/df_nat/df_full), so their
## AIC/BIC/n are NOT comparable across columns. (In C4 all six models share
## df_full, so there the information criteria ARE comparable.)
note_diff_samples <- paste(
  "Cells are odds ratios with 95% robust (HC0) CIs in brackets.",
  "Columns are estimated on different analytic samples;",
  "AIC/BIC/n are not comparable across columns."
)

## ----------------------------------------------------------------------------
## 7.1a Table C1 - all respondents (non-voters coded 0)
## ----------------------------------------------------------------------------
or_table(
  list(ml1, ml2, ml3, ml4, ml5),
  out       = file.path(output_dir, "table_C1_all_respondents.html"),
  dep_label = "Vote Choice (Smer, SNS, Republika)",
  notes     = note_diff_samples
)

## ----------------------------------------------------------------------------
## 7.1b Table C2 - voters only (non-voters excluded)
## ----------------------------------------------------------------------------
ml1_v <- update(ml1, data = filter(df_illib, vote_outcome != "non-voter"))
ml2_v <- update(ml2, data = filter(df_auth,  vote_outcome != "non-voter"))
ml3_v <- update(ml3, data = filter(df_pop,   vote_outcome != "non-voter"))
ml4_v <- update(ml4, data = filter(df_nat,   vote_outcome != "non-voter"))
ml5_v <- update(ml5, data = filter(df_full,  vote_outcome != "non-voter"))

or_table(
  list(ml1_v, ml2_v, ml3_v, ml4_v, ml5_v),
  out       = file.path(output_dir, "table_C2_voters_only.html"),
  dep_label = "Vote Choice (Smer, SNS, Republika) - voters only",
  notes     = note_diff_samples
)

## ----------------------------------------------------------------------------
## 7.1c Table C3 - Hlas counted as an illiberal party
## ----------------------------------------------------------------------------
## Refits each model on the same samples with the Hlas-inclusive outcome.
ml1_h <- update(ml1, vote_illiberal_hlas ~ .)
ml2_h <- update(ml2, vote_illiberal_hlas ~ .)
ml3_h <- update(ml3, vote_illiberal_hlas ~ .)
ml4_h <- update(ml4, vote_illiberal_hlas ~ .)
ml5_h <- update(ml5, vote_illiberal_hlas ~ .)

or_table(
  list(ml1_h, ml2_h, ml3_h, ml4_h, ml5_h),
  out       = file.path(output_dir, "table_C3_hlas_included.html"),
  dep_label = "Vote Choice (Smer, SNS, Republika, Hlas)",
  notes     = note_diff_samples
)

## ----------------------------------------------------------------------------
## Helper: odds-ratio coefficient plot from a fitted glm (Figures 1, C1-C3).
## ----------------------------------------------------------------------------
or_plot <- function(model, label_map, level_order) {
  se <- sqrt(diag(vcovHC(model, type = "HC0")))
  b  <- coef(model)
  data.frame(
    term = names(b), estimate = exp(b),
    conf.low = exp(b - 1.96 * se), conf.high = exp(b + 1.96 * se)
  ) %>%
    filter(term != "(Intercept)", term %in% names(label_map)) %>%
    mutate(predictor = factor(label_map[term], levels = level_order)) %>%
    ggplot(aes(x = predictor, y = estimate)) +
    geom_point(shape = 5) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                  width = 0.2, linewidth = 0.3) +
    geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.3) +
    coord_flip(ylim = c(0, 3)) +
    labs(x = NULL, y = expression(italic("Odds Ratio (95% Robust CI)"))) +
    theme_bw() +
    theme(
      text        = element_text(family = "Times New Roman", size = 12),
      axis.text.y = element_text(face = "italic")
    )
}

## ----------------------------------------------------------------------------
## 7.2a Figure 1 (main) - full model odds ratios
## ----------------------------------------------------------------------------
labels_main <- c(
  government_constraints = "Unconstrained Government",
  political_rights       = "Political Rights Curtailment",
  state_neutral          = "State Neutral Opposition",
  closed_society         = "Closed Society",
  democracy              = "Anti-democracy",
  authoritarianism       = "Authoritarianism",
  pca_populism1          = "Populism",
  anti_immigration       = "Nativism"
)
order_main <- rev(c(
  "Unconstrained Government", "Political Rights Curtailment",
  "State Neutral Opposition", "Closed Society", "Anti-democracy",
  "Authoritarianism", "Populism", "Nativism"
))

figure1 <- or_plot(ml5, labels_main, order_main)
ggsave(file.path(output_dir, "fig_1_odds_ratios.png"), figure1, width = 7, height = 5, dpi = 300)

## ----------------------------------------------------------------------------
## 7.2b Figure C1 - political-rights validation item (political_rights2)
## ----------------------------------------------------------------------------
ml_val <- glm(vote_illiberal ~ government_constraints + political_rights2 +
                state_neutral + closed_society + democracy + authoritarianism +
                pca_populism1 + anti_immigration +
                female + age + edu + left_right + interest + satisfaction_democracy + trust,
              weights = weight_FINAL, data = df_full, family = binomial)

labels_val <- labels_main
names(labels_val)[names(labels_val) == "political_rights"] <- "political_rights2"
labels_val["political_rights2"] <- "Political Rights Curtailment (Validation)"
order_val <- rev(c(
  "Unconstrained Government", "Political Rights Curtailment (Validation)",
  "State Neutral Opposition", "Closed Society", "Anti-democracy",
  "Authoritarianism", "Populism", "Nativism"
))

figure_C1 <- or_plot(ml_val, labels_val, order_val)
ggsave(file.path(output_dir, "fig_C1_political_rights_validation.png"),
       figure_C1, width = 7, height = 5, dpi = 300)

## ----------------------------------------------------------------------------
## 7.2c Table C4 - populism operationalisation robustness
## ----------------------------------------------------------------------------
## Each model swaps in a different populism operationalisation; all other terms
## are identical. All six share df_full, so their AIC/BIC ARE comparable.
mp1 <- glm(vote_illiberal ~ government_constraints + political_rights + state_neutral +
             closed_society + democracy + authoritarianism + pca_populism1 + anti_immigration +
             female + age + edu + left_right + interest + satisfaction_democracy + trust,
           weights = weight_FINAL, data = df_full, family = binomial)

mp2 <- glm(vote_illiberal ~ government_constraints + political_rights + state_neutral +
             closed_society + democracy + authoritarianism + pca_populism2 + anti_immigration +
             female + age + edu + left_right + interest + satisfaction_democracy + trust,
           weights = weight_FINAL, data = df_full, family = binomial)

mp3 <- glm(vote_illiberal ~ government_constraints + political_rights + state_neutral +
             closed_society + democracy + authoritarianism + populism_additive + anti_immigration +
             female + age + edu + left_right + interest + satisfaction_democracy + trust,
           weights = weight_FINAL, data = df_full, family = binomial)

mp4 <- glm(vote_illiberal ~ government_constraints + political_rights + state_neutral +
             closed_society + democracy + authoritarianism + populism_dummy + anti_immigration +
             female + age + edu + left_right + interest + satisfaction_democracy + trust,
           weights = weight_FINAL, data = df_full, family = binomial)

mp5 <- glm(vote_illiberal ~ government_constraints + political_rights + state_neutral +
             closed_society + democracy + authoritarianism +
             people_center + anti_elit_bus + manichean + anti_immigration +
             female + age + edu + left_right + interest + satisfaction_democracy + trust,
           weights = weight_FINAL, data = df_full, family = binomial)

mp6 <- glm(vote_illiberal ~ government_constraints + political_rights + state_neutral +
             closed_society + democracy + authoritarianism +
             people_center + anti_elit_pol + manichean + anti_immigration +
             female + age + edu + left_right + interest + satisfaction_democracy + trust,
           weights = weight_FINAL, data = df_full, family = binomial)

or_table(
  list(mp1, mp2, mp3, mp4, mp5, mp6),
  out       = file.path(output_dir, "table_C4_populism_robustness.html"),
  dep_label = "Vote Choice (Smer, SNS, Republika)",
  covariate.labels = c(
    "Unconstrained Government", "Political Rights Curtailment",
    "State Neutral Opposition", "Closed Society", "Anti-Democracy",
    "Authoritarianism", "PCA Populism", "PCA Populism (robustness)",
    "Populism (additive index)", "Populism (dummy)", "People-centeredness",
    "Anti-elitism (anti-business elite)", "Anti-elitism (anti-politician)",
    "Manichean", "Nativism", "Female", "Age", "Education", "Left-right Ideology",
    "Interest in Politics", "Satisfaction with Democracy", "Trust"
  ),
  notes = paste(
    "Cells are odds ratios with 95% robust (HC0) CIs in brackets.",
    "All columns share the df_full sample; AIC/BIC are comparable across columns."
  )
)

## ----------------------------------------------------------------------------
## 7.2d Figure C2 - national purity replacing nativism
## ----------------------------------------------------------------------------
ml_np1 <- glm(vote_illiberal ~ government_constraints + political_rights +
                state_neutral + closed_society + democracy + authoritarianism +
                pca_populism1 + national_purity +
                female + age + edu + left_right + interest + satisfaction_democracy + trust,
              weights = weight_FINAL, data = df_full_np, family = binomial)

labels_np <- c(labels_main[c("government_constraints", "political_rights",
                             "state_neutral", "closed_society", "democracy",
                             "authoritarianism", "pca_populism1")],
               national_purity = "National Purity")
order_np <- rev(c(
  "Unconstrained Government", "Political Rights Curtailment",
  "State Neutral Opposition", "Closed Society", "Anti-democracy",
  "Authoritarianism", "Populism", "National Purity"
))

figure_C2 <- or_plot(ml_np1, labels_np, order_np)
ggsave(file.path(output_dir, "fig_C2_national_purity_only.png"),
       figure_C2, width = 7, height = 5, dpi = 300)

## ----------------------------------------------------------------------------
## 7.2e Figure C3 - nativism AND national purity jointly
## ----------------------------------------------------------------------------
ml_np2 <- glm(vote_illiberal ~ government_constraints + political_rights +
                state_neutral + closed_society + democracy + authoritarianism +
                pca_populism1 + anti_immigration + national_purity +
                female + age + edu + left_right + interest + satisfaction_democracy + trust,
              weights = weight_FINAL, data = df_full_np, family = binomial)

labels_np2 <- c(labels_main, national_purity = "National Purity")
order_np2 <- rev(c(
  "Unconstrained Government", "Political Rights Curtailment",
  "State Neutral Opposition", "Closed Society", "Anti-democracy",
  "Authoritarianism", "Populism", "Nativism", "National Purity"
))

figure_C3 <- or_plot(ml_np2, labels_np2, order_np2)
ggsave(file.path(output_dir, "fig_C3_nativism_plus_purity.png"),
       figure_C3, width = 7, height = 5, dpi = 300)

## ----------------------------------------------------------------------------
## 7.2f Figure C4 - predicted probabilities from the binary full model (ml5)
## ----------------------------------------------------------------------------
## Predicts P(illiberal vote) across each focal attitude, other predictors held
## at their means (modal category for the binary control `female`).

control_means <- df_full %>%
  summarise(across(c(government_constraints, political_rights, state_neutral,
                     closed_society, democracy, authoritarianism, pca_populism1,
                     anti_immigration, age, edu, left_right, interest,
                     satisfaction_democracy, trust),
                   ~ mean(.x, na.rm = TRUE))) %>%
  mutate(female = factor(names(which.max(table(df_full$female))),
                         levels = levels(df_full$female)))

predict_over <- function(model, focal) {
  n_lev <- length(unique(df_full[[focal]]))
  grid <- tibble(value = sort(unique(df_full[[focal]]))) %>%
    bind_cols(control_means[rep(1, n_lev), ]) %>%
    mutate(!!focal := value)

  X <- model.matrix(delete.response(terms(model)), data = grid,
                    xlev = model$xlevels)
  X <- X[, names(coef(model)), drop = FALSE]      # align to coefficient order
  V <- vcovHC(model, type = "HC0")

  eta <- as.vector(X %*% coef(model))
  se  <- sqrt(rowSums((X %*% V) * X))             # robust SE on the link scale

  grid %>% mutate(predicted_prob = plogis(eta),
                  lower = plogis(eta - 1.96 * se),
                  upper = plogis(eta + 1.96 * se))
}

pp_panel <- function(df, focal, xlab) {
  ggplot(df, aes(x = .data[[focal]], y = predicted_prob)) +
    geom_point(color = "black", size = 2, shape = 5) +
    geom_line(color = "black", linetype = "dashed") +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, linewidth = 0.3) +
    labs(x = xlab, y = expression(italic("Predicted Probability"))) +
    coord_cartesian(ylim = c(0, 0.5)) +
    theme_bw() +
    theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10),
          axis.text.y = element_text(face = "italic"))
}

pred_plot1 <- pp_panel(predict_over(ml5, "authoritarianism"),
                       "authoritarianism", expression(italic("Authoritarianism Attitudes")))
pred_plot2 <- pp_panel(predict_over(ml5, "anti_immigration"),
                       "anti_immigration", expression(italic("Nativist Attitudes")))

ggsave(file.path(output_dir, "fig_C4_predicted_probabilities.png"),
       pred_plot1 | pred_plot2, width = 10, height = 5, dpi = 300)


## ============================================================================
## 8. MULTINOMIAL MODEL (Appendix Table D1)
## ============================================================================

df_mnl <- df_full %>%
  mutate(vote_outcome_multinom = case_when(
    vote_outcome == "Smer"                                        ~ "Smer",
    vote_outcome == "Republika"                                   ~ "Republika",
    vote_outcome == "SNS"                                         ~ "SNS",
    vote_outcome == "Hlas"                                        ~ "Hlas",
    vote_outcome %in% c("PS", "KDH", "SAS", "OLANO", "Demokrati") ~ "Opposition",
    vote_outcome == "non-voter"                                   ~ "non-voter",
    vote_outcome %in% c("other", "Aliancia", "Sme Rodina")        ~ "other",
    TRUE                                                          ~ NA_character_
  )) %>%
  filter(vote_outcome_multinom != "other") %>%
  mutate(vote_outcome_multinom = relevel(factor(vote_outcome_multinom), ref = "non-voter"))

set.seed(2023)
mult_model <- multinom(
  vote_outcome_multinom ~ government_constraints + political_rights +
    state_neutral + closed_society + democracy + authoritarianism +
    pca_populism1 + anti_immigration +
    female + age + edu + left_right + interest + satisfaction_democracy + trust,
  weights = weight_FINAL, data = df_mnl, model = TRUE, maxit = 1000
)
if (!identical(mult_model$convergence, 0L))
  warning("mult_model did not converge (convergence code ",
          mult_model$convergence, "); inspect before trusting Table D1.")

mult_fit <- list(
  pseudoR = PseudoR2(mult_model, which = "McFadden"),
  AIC     = PseudoR2(mult_model, which = "AIC"),
  BIC     = PseudoR2(mult_model, which = "BIC"),
  logLik  = PseudoR2(mult_model, which = "logLik")
)

mult_se <- summary(mult_model)$standard.errors
mult_p  <- 2 * pnorm(abs(coef(mult_model) / mult_se), lower.tail = FALSE)

stargazer(
  mult_model, type = "html",
  covariate.labels = c(
    "Unconstrained Government", "Political Rights Curtailment",
    "State Neutral Opposition", "Closed Society", "Anti-Democracy",
    "Authoritarianism", "Populism", "Nativism",
    "Female", "Age", "Education", "Left-right Ideology",
    "Interest in Politics", "Satisfaction with Democracy", "Trust"
  ),
  coef = list(exp(coef(mult_model))),   # relative risk ratios
  p    = list(mult_p),                  ## PATCH: model-based significance stars
  p.auto = FALSE, float = FALSE, header = FALSE, font.size = "scriptsize",
  add.lines = list(
    c("n",                     nrow(df_mnl)),
    c("Pseudo R2 (McFadden)",  round(mult_fit$pseudoR, 3)),
    c("AIC",                   round(mult_fit$AIC, 2)),
    c("BIC",                   round(mult_fit$BIC, 2)),
    c("logLik",                round(mult_fit$logLik, 2))
  ),
  omit.stat = "aic",
  notes = paste("Cells are relative-risk ratios (exp(beta)).",
                "Stars from model-based (non-robust) p-values."),
  out = file.path(output_dir, "table_D1_multinomial.html")
)


## ============================================================================
## 9. AVERAGE MARGINAL EFFECTS (multinomial)
## ============================================================================

predictor_labels <- labels_main
predictor_order  <- c(
  "Unconstrained Government", "Political Rights Curtailment",
  "State Neutral Opposition", "Closed Society", "Anti-democracy",
  "Authoritarianism", "Populism", "Nativism"
)

ame_plot_for <- function(model, groups, file) {
  avg_slopes(model) %>%
    filter(group %in% groups, term %in% names(predictor_labels)) %>%
    mutate(predictor = factor(predictor_labels[term], levels = predictor_order)) %>%
    ggplot(aes(y = reorder(predictor, estimate), x = estimate)) +
    geom_errorbar(aes(xmin = conf.low, xmax = conf.high), width = 0.4, colour = "black") +
    geom_point(size = 0.9, colour = "black") +
    facet_wrap(~ group) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_y_discrete(limits = rev(predictor_order)) +
    coord_cartesian(xlim = c(-0.15, 0.15)) +
    labs(x = "AME (95% CI)", y = "") +
    theme_bw() +
    theme(
      text        = element_text(family = "Times New Roman", size = 12),
      axis.text.y = element_text(face = "italic")
    ) -> p
  ggsave(file, p, width = 10, height = 6, dpi = 300)
  p
}

## 9.1 Figure 2 (main) - Smer, SNS, Republika, Opposition
ame_plot_for(mult_model, c("Smer", "SNS", "Republika", "Opposition"),
             file.path(output_dir, "fig_2_ame_main.png"))

## 9.2 Figure D1 - Hlas only
ame_plot_for(mult_model, "Hlas",
             file.path(output_dir, "fig_D1_ame_hlas.png"))


## ============================================================================
## 10. PREDICTED PROBABILITIES, MULTINOMIAL (Appendix Figures D2-D3)
## ============================================================================

plot_outcomes <- c("Smer", "Opposition", "Republika", "SNS")

mnl_pp <- function(focal, xlab, ymax, file) {
  ggemmeans(mult_model, terms = paste0(focal, " [all]")) %>%
    filter(response.level %in% plot_outcomes) %>%
    ggplot(aes(x = x, y = predicted)) +
    facet_wrap(~ response.level) +
    geom_line(color = "black", linetype = "dashed") +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
    labs(x = xlab, y = expression(italic("Predicted Probability"))) +
    coord_cartesian(ylim = c(0, ymax)) +
    theme_bw() +
    theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10),
          axis.text.y = element_text(face = "italic")) -> p
  ggsave(file, p, width = 8, height = 6, dpi = 300)
  p
}

mnl_pp("authoritarianism", expression(italic("Authoritarianism Attitudes")), 0.5,
       file.path(output_dir, "fig_D2_pp_authoritarianism.png"))
mnl_pp("anti_immigration", expression(italic("Nativist Attitudes")), 0.7,
       file.path(output_dir, "fig_D3_pp_nativism.png"))


## ============================================================================
## 11. MULTINOMIAL ROBUSTNESS: NATIONAL PURITY (Appendix Figure D4)
## ============================================================================

df_mnl_robust <- df_mnl %>%
  filter(between(q38, 1, 5)) %>%
  mutate(national_purity = 6 - q38)

set.seed(2023)
mult_model_robust <- multinom(
  vote_outcome_multinom ~ government_constraints + political_rights +
    state_neutral + closed_society + democracy + authoritarianism +
    pca_populism1 + national_purity + anti_immigration +
    female + age + edu + left_right + interest + satisfaction_democracy + trust,
  weights = weight_FINAL, data = df_mnl_robust, model = TRUE, maxit = 1000
)
if (!identical(mult_model_robust$convergence, 0L))
  warning("mult_model_robust did not converge (convergence code ",
          mult_model_robust$convergence, ").")

predictor_labels_r <- c(predictor_labels, national_purity = "National Purity")
predictor_order_r  <- c(predictor_order, "National Purity")

ame_df_robust <- avg_slopes(mult_model_robust, wts = "weight_FINAL") %>%   
  filter(group %in% c("Smer", "SNS", "Republika", "Opposition"),
         term %in% names(predictor_labels_r)) %>%
  mutate(predictor = factor(predictor_labels_r[term], levels = predictor_order_r))

ame_plot_robust <- ame_df_robust %>%
  ggplot(aes(y = reorder(predictor, estimate), x = estimate)) +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high), width = 0.4, colour = "black") +
  geom_point(size = 0.9, colour = "black") +
  facet_wrap(~ group) +
  scale_x_continuous(limits = c(-0.15, 0.15)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_y_discrete(limits = rev(predictor_order_r)) +
  labs(x = "AME (95% CI)", y = "") +
  theme_bw()

ggsave(file.path(output_dir, "fig_D4_ame_national_purity.png"),
       ame_plot_robust, width = 10, height = 6, dpi = 300)


## ============================================================================
## 12. SESSION INFO (reproducibility)
## ============================================================================
writeLines(capture.output(sessionInfo()),
           file.path(output_dir, "sessionInfo.txt"))
