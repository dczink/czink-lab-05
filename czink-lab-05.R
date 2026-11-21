
# install.packages("dplyr")
# install.packages("tidyverse")
# install.packages("ggplot2")
# install.packages("glmmTMB")
# install.packages("DHARMa")
# install.packages("ggeffects")
# install.packages("emmeans")
# install.packages("broom.mixed")
# install.packages("gt")
# install.packages("AER")

library(dplyr)
library(tidyverse)
library(ggplot2)
library(glmmTMB)
library(DHARMa)
library(ggeffects)
library(emmeans)
library(broom.mixed)
library(gt)
library(AER)

# Data loading, used data function because simply calling didnt work
data("NMES1988")

# select variables and assign reference groups for factor variables
data.clean <- NMES1988 |>
  select(visits, chronic, health, gender, school, insurance) |>
  drop_na() |>
  mutate(
    gender = fct_relevel(factor(gender), "male"),
    health = fct_relevel(factor(health), "average"),
    insurance = fct_relevel(factor(insurance), "no")
  )

n <- nrow(data.clean)

## Explore the data

outcome.summary <- data.clean |>
  summarise(
    n      = n(),
    mean   = mean(visits),
    var    = var(visits),
    p_zero = mean(visits == 0),
    max    = max(visits)
  )

outcome.summary

descriptive.plot <- ggplot(data.clean, aes(x = visits)) +
  geom_bar() +
  labs(x = "Number of Visits", y = "Number of People")

descriptive.plot

## Why not a linear model?

model.gauss <- glmmTMB(
  visits ~ chronic + health + gender + school + insurance,
  data = data.clean,
  family = gaussian(link = "identity")
)

#simulation 
set.seed(9041)
res.gauss <- simulateResiduals(model.gauss, n = 500, plot = FALSE)
plot(res.gauss)

## Nested Poisson Models

null <- glmmTMB(
  visits ~ 1, 
  data = data.clean, 
  family = poisson(link = "log"))

predict1 <- glmmTMB(
  visits ~ chronic, 
  data = data.clean, 
  family = poisson(link = "log"))

poisson <- glmmTMB(
  visits ~ chronic + health + gender + school + insurance,
  data = data.clean,
  family = poisson(link = "log"))

summary(null)
summary(predict1)
summary(poisson)

## Nested Model Comparisons

### Compare the models

tests <- anova(null, predict1, poisson)

tests <- tests |>
  select(Chisq, `Chi Df`, `Pr(>Chisq)`) |>
  rename(
    Nested.df = `Chi Df`,
    Nested.Chisq.stat = Chisq,
    Nested.p.value = `Pr(>Chisq)`)

### Overall model performance

glance1 <- broom.mixed::glance(null)
glance2 <- broom.mixed::glance(predict1)
glance3 <- broom.mixed::glance(poisson)

models <- c(
  "Null",
  "Predict1: Chronic",
  "Full: All Variables"
  )

glances <- rbind(
  glance1,
  glance2,
  glance3)

### Combine nested tests and model performance

comparison.table <- cbind(
  Model = models,
  glances,
  tests)

comparison.table |>
  select(
    Model,
    AIC,
    BIC,
    logLik,
    deviance,
    df.residual,
    Nested.df,
    Nested.Chisq.stat,
    Nested.p.value
  ) |>
  gt() |>
  cols_label(
    Model = "Model",
    AIC = "AIC",
    BIC = "BIC",
    logLik = "Log-likelihood",
    deviance = "Deviance",
    df.residual = "Residual.df",
    Nested.df = "Nested df",
    Nested.Chisq.stat = "Nested χ²",
    Nested.p.value = "Nested p"
  ) |>
  fmt_number(
    columns = c(
      AIC,
      BIC,
      logLik,
      deviance,
      Nested.Chisq.stat
    ),
    decimals = 2
  ) |>
  fmt_number(
    columns = Nested.p.value,
    decimals = 3
  ) |>
  fmt(
    columns = Nested.p.value,
    fns = function(x) {
      ifelse(
        is.na(x),
        "—",
        ifelse(x < .001, "< .001", sprintf("%.3f", x))
      )
    }
  ) |>
  sub_missing(
    columns = c(Nested.df, Nested.Chisq.stat),
    missing_text = "—"
  )

## Select a provisional model

selected.row <- comparison.table |>
  filter(Nested.p.value < .05) |>
  slice_tail(n = 1)

selected.name <- selected.row$Model

selected.model <- switch(
  selected.name,
  "Null" = null,
  "Predict1: Chronic" = predict1,
  "Full: All Variables" = poisson)

## Checking the Poisson

set.seed(9041)

res.pois <- simulateResiduals(poisson, n = 500, plot = FALSE)

plot(res.pois)

disp.pois <- testDispersion(res.pois)
disp.pois


zi.pois <- testZeroInflation(res.pois)
zi.pois

## The negative binomial

# Negative binomial (NB2) model
neg.binom <- glmmTMB(
  visits ~ chronic + health + gender + school + insurance,
  data   = data.clean,
  family = nbinom2(link = "log")
)

summary(neg.binom)

# Does the NB fit better than the Poisson?
lrt.nb<-anova(poisson, neg.binom)
lrt.nb

### Checking the NB

set.seed(9041)

res.nb <- simulateResiduals(neg.binom, n = 500, plot = FALSE)

plot(res.nb)

disp.nb <- testDispersion(res.nb)
disp.nb

zi.nb <- testZeroInflation(res.nb)
zi.nb

### Visualize predicted counts

pred.nb <- predict_response(neg.binom, terms = "chronic")

plot(pred.nb) +
  labs(
    title = "Predicted Visits by Number of Chronic Conditions",
    x     = "Chronic Conditions", y = "Visits"
  )


## Do the zeros need their own model?

# Zero-inflated NB with a constant zero-inflation probability
m.zinb <- glmmTMB(
  visits ~ chronic + health + gender + school + insurance,
  ziformula = ~ 1,
  data = data.clean, family = nbinom2()
)

# Zero-inflated NB with predictors in the zero part
m.zinb.full <- glmmTMB(
  visits ~ chronic + health + gender + school + insurance,
  ziformula = ~ chronic + health + gender + school + insurance,
  data = data.clean, family = nbinom2()
)

summary(m.zinb.full) 

# Hurdle NB
m.hurd <- glmmTMB(
  visits ~ chronic + health + gender + school + insurance,
  ziformula = ~ chronic + health + gender + school + insurance,
  data = data.clean, family = truncated_nbinom2()
)

# Compare nested models
lrt.zero <- anova(neg.binom, m.zinb.full)

# compare hurdle with AIC/BIC 
AIC(neg.binom, m.zinb, m.zinb.full, m.hurd)
BIC(neg.binom, m.zinb, m.zinb.full, m.hurd)

# Then re-check the chosen model with DHARMa
set.seed(9041)
res.zinb <- simulateResiduals(m.zinb.full, n = 500, plot = FALSE)
zi.zi <- testZeroInflation(res.zinb)

## Adding an interaction: Number of Chronic conditions X insurance coverage

### Step 1: Fit the model with the product term

interact <- glmmTMB(
  visits ~ school + health + gender + chronic * insurance,
  data = data.clean, family = nbinom2()
)

summary(interact)

### Step 2: Test whether it's needed

# LRT: does adding the phd x fem interaction improve fit over the additive NB?
lrt.int <- anova(m.zinb.full, interact)
lrt.int

### Step 3: Plot predicted counts with ggeffects

# Predicted mean article counts from the interaction model, on the count scale.
pred.int <- predict_response(interact, terms = c("chronic", "insurance"))

plot(pred.int) +
  scale_colour_manual(values = c(no = "orange", yes = "blue")) +
  scale_fill_manual(values = c(no = "orange", yes = "blue")) +
  labs(
    title  = "Predicted Visits by Chronic Conditions and Insurance",
    x      = "Chronic Conditions", y = "Visits",
    colour = NULL, fill = NULL
  )

### Step 4a: Simple slopes with emmeans

slopes <- emtrends(
  interact, 
  specs = ~ insurance, 
  var = "chronic")

pairs(slopes)

#test simple slopes against zero
slopes.sum <- summary(slopes, infer = TRUE)
slopes.sum

### Step 4b: Simple effects at chosen values of prestige

# Pick the probing values: low (-1 SD), average (mean) and high (+1 SD) prestige.
chronic.mean <- mean(data.clean$chronic, na.rm = TRUE)
chronic.sd <- sd(data.clean$chronic, na.rm = TRUE)
chronic.values <- chronic.mean + c(-1, 0, 1) * chronic.sd

# Estimated marginal means (EMMs) for each gender at each prestige value.
insurance.emmeans <- emmeans(
  interact,
  specs = ~ insurance | chronic,
  at = list(chronic = chronic.values)
)

# Estimated mean article counts in each gender x prestige cell.
summary(insurance.emmeans, type = "response")

# Gender differences at each prestige value, on the log scale.
pairs(insurance.emmeans)

### Step 5: Report the simple effects in a table

# Build a data frame of the gender simple effects, then format it as a table.
contrast.df <- summary(pairs(insurance.emmeans)) |>
  as.data.frame() |>
  mutate(
    chronic.level = c("-1 SD", "Mean", "+1 SD")[match(round(chronic, 6), round(chronic.values, 6))]
  ) |>
  select(chronic.level, contrast, estimate, SE, z.ratio, p.value)

contrast.df |>
  gt() |>
  cols_label(
    chronic.level = "Chronic Conditions", contrast = "Contrast",
    estimate = "Log Contrast", SE = "SE", z.ratio = "z", p.value = "p"
  ) |>
  fmt_number(columns = c(estimate, SE, z.ratio), decimals = 3) |>
  fmt(columns = p.value,
      fns = \(x) scales::pvalue(x, accuracy = 0.001, add_p = FALSE)) |>
  tab_header(
    title    = "Simple Effects of Insurance Across Chronic Conditions",
    subtitle = "Full Zero-Inflated NB"
  ) |>
  tab_source_note(
    "Log contrasts = No − Yes. z-tests and p-values are computed on the log scale."
  )

------------------------------------------------------------------------
  
  ## 11. Results write-up
  
  Article counts were overdispersed (*M* = 1.69, *s*² = 3.71; 30% zeros). Simulation-based residual checks on a Poisson model showed significant overdispersion (dispersion ratio = 1.99, *p* \< .001) and excess zeros (observed/expected = 1.44, *p* \< .001). A negative binomial model fit significantly better, χ²(1) = 180.20, *p* \< .001, and showed no remaining evidence of overdispersion (ratio = 1.05, *p* = .556) or excess zeros (ratio = 0.99, *p* = .772), so no separate model for zeros was needed.

Holding other predictors constant, women published fewer articles than men (*b* = −0.22, 95% CI \[−0.36, −0.07\], *p* = .003). Each additional child aged 5 or younger was associated with fewer articles (*b* = −0.18, 95% CI \[−0.28, −0.07\], *p* \< .001), and each additional mentor article with more articles (*b* = 0.029, 95% CI \[0.023, 0.035\], *p* \< .001). Marital status (*b* = 0.15, 95% CI \[−0.01, 0.31\], *p* = .067) and PhD prestige (*b* = 0.02, 95% CI \[−0.06, 0.09\], *p* = .672) were not significantly associated with article counts.

Adding a PhD prestige × gender interaction did not significantly improve fit, χ²(1) = 3.73, *p* = .053. Simple slopes showed that prestige was not significantly associated with article counts for men (*b* = −0.04, 95% CI \[−0.14, 0.05\], *p* = .363) or women (*b* = 0.10, 95% CI \[−0.01, 0.20\], *p* = .084). Simple effects, however, showed that the gender gap narrowed as prestige increased. Men were expected to publish more articles than women at low prestige (−1 *SD*; *b* = 0.36, *p* \< .001) and at average prestige (*b* = 0.22, *p* = .002), but the difference was not significant at high prestige (+1 *SD*; *b* = 0.09, *p* = .354).

------------------------------------------------------------------------
  
  ## Reference
  
  Long, J. S. (1990). The origins of sex differences in science. *Social Forces, 68*(4), 1297–1316.
