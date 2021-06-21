library(tidyverse)
library(brms)
library(tidybayes)
theme_set(theme_tidybayes())

# data entry ----

regenron <- tribble(
  ~"group", ~"intervention", ~"death", ~"n",
  "seropos", "regenron", 411, 2636,
  "seropos", "control", 383, 2636,
  "seroneg", "regenron", 396, 1633,
  "seroneg", "control", 451, 1520,
  "unknown", "regenron", 137, 570,
  "unknown", "control", 192, 790
) %>% mutate(across(where(is.character), as.factor))

# prior sampling for prior prediction
mr_prior <- brm(death | trials(n) ~ -1 + group + intervention + intervention * group,
  family = binomial(), data = regenron,
  prior = prior(normal(0, 1), class = "b"),
  sample_prior = "only"
)

mr <- brm(death | trials(n) ~ -1 + group + intervention + intervention * group, 
          family = binomial(), 
          data = regenron, 
          prior = prior(normal(0, 1), class = "b"))

regenron %>% # original data
  modelr::data_grid(group, intervention) %>% # generate new data
  mutate(n = 1) %>% # n is the number of trials
  add_fitted_draws(mr) %>% # compute draws from the linear predictor for model `mr`, replace with `mr_prior` for prior prediction
  mutate(.value = brms::inv_logit_scaled(.value)) %>% # inverse logit transformation to probability scale
  compare_levels(.value, by = intervention, fun = `-`) %>% # calculate absolute risk reduction per group
  ggplot(aes(x = .value, y = group, fill = after_stat(ifelse(x > 0, "over", "under")))) +
  stat_halfeye() +
  scale_fill_manual(values = c("over" = "pink", "under" = "skyblue")) +
  theme_tidybayes() +
  theme(legend.position = "none") +
  scale_y_discrete(name = "Group") +
  scale_x_continuous(name = "Absolute risk difference (28d mortality)") +
  labs(title = "REGEN-COV", subtitle = "Bayesian logistic regression, minimally informative sceptical prior", caption = "@load_dependent")



mr <- brm(death | trials(n) ~ -1 + group + intervention + intervention * group, family = binomial(), data = regenron, prior = prior(normal(0, 1), class = "b"))