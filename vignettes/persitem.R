library(tidyverse)
library(mirt)
source("R/helpers.R")

data <-
	simdata(a = rep(1, 10), d = rnorm(10), itemtype = "2PL", N = 100) %>%
	as_tibble()

# this does stratified-by-person masking
mask_item_response_data <- function(data, percent_mask){
	notna <- which(!is.na(as.numeric(data)))
	makena <- sample(notna, round(percent_mask * length(notna)))
	data[makena] <- NA
	data
}

# now let's make train and test data
train <- data %>% apply(1, mask_item_response_data, 0.3) %>% t()
test <- data
test[!is.na(data) & !is.na(train)] <- NA

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

cv <-
	models %>%
	select(dim:method) %>%
	mutate(
		in_model = pmap(., ~ mirt(train, ..1, ..2, method = ..3, technical = list(NCYCLES = 500))),
		fscores = map2(in_model, method, ~ fscores(.x, QMC = .y == "QMCEM", rotate = "none")),
		p = map2(in_model, fscores, get_p),
		eval = p %>% map(evaluate_p, train, test)
	)

cv %>%
	unnest_wider(eval)
