library(tidyverse)
library(rsample)
library(mirt)
source("R/helpers.R")

# generate some generic item response data
data <-
	simdata(a = rep(1, 10), d = rnorm(10), itemtype = "2PL", N = 100) %>%
	as_tibble()

# fit each model to the full data
models <-
	tibble(
		dim = c(1, 1, 2),
		itemtype = c("Rasch", "2PL", "2PL"),
		method = c("EM", "EM", "EM")
	) %>%
	mutate(
		model = pmap(., ~ mirt(data, ..1, ..2, method = ..3, technical = list(NCYCLES = 500)))
	) %>%
	mutate(nestpars = model %>% map_dbl(~ .@Model$nestpars)) %>%
	arrange(nestpars) %>%
	mutate(
		in_sample_log_lik = model %>% map_dbl(~ .@Fit$logLik),
		anova_p = map2_dbl(model, lag(model), get_p_from_anova),
		aic = model %>% map_dbl(~ .@Fit$AIC),
		bic = model %>% map_dbl(~ .@Fit$BIC)
	)

# setup folds based on persons
splits <-
	vfold_cv(data, v = 3, repeats = 2) %>%
	select(starts_with("id"), splits)

# run cv by fitting each model to the train data and assess on test data
cv <-
	crossing(
		models %>% select(dim:method),
		splits %>% select(-splits)
	) %>%
	mutate(
		in_model =
			pmap(
				list(dim, itemtype, method, map2(id, id2, get_split, splits)),
				~ mirt(training(..4), ..1, ..2, method = ..3, technical = list(NCYCLES = 500)) # NCYCLES
			),

		out_data = map2(id, id2, ~ testing(get_split(.x, .y, splits))),

		log_lik =
			pmap_dbl(
				list(method, in_model, out_data),
				~ mirt_calc_log_lik(..1, ..2, ..3, draws = 100000) %>% sum()
			)
	)

# see what we got
cv %>%
	group_by(dim, itemtype, method) %>%
	summarize(out_log_lik = mean(log_lik)) %>%
	arrange(desc(out_log_lik))
