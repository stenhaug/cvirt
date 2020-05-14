
get_p_from_anova <- function(mod1, mod2){
	if("logical" %in% c(class(mod1), class(mod2))) {return(NA_real_)} # catches missing
	anova(mod1, mod2)$p[2]
}

get_split <- function(id, id2, splits){
	splits$splits[[which(splits$id == id & splits$id2 == id2)]]
}

mirt_calc_log_lik <- function(method, model, data, draws){
	if(method == "EM"){
		return(mirt_em_calc_log_lik_marg(model, data))
	}

	if(method == "QMCEM"){
		return(mirt_mci_calc_log_lik_marg(model, data, draws))
	}
}

mirt_em_calc_log_lik_marg <- function(model, data){
	mirt:::Estep.mirt(
		pars = model@ParObjects$pars,
		tabdata = make_fulldata(data),
		freq = rep(1, nrow(data)),
		CUSTOM.IND = model@Internals$CUSTOM.IND,
		Theta = model@Model$Theta,
		prior = model@Internals$Prior[[1]],
		itemloc = model@Model$itemloc,
		full = FALSE,
		Etable = TRUE
	)$expected %>%
		log()
}

mirt_mci_calc_log_lik_marg <- function(model, data, draws){
	dim <- ncol(model@Fit$F)

	mirt:::Estep.mirt(
		pars = model@ParObjects$pars,
		tabdata = make_fulldata(data),
		freq = rep(1, nrow(data)),
		CUSTOM.IND = model@Internals$CUSTOM.IND,
		Theta = MASS::mvrnorm(n = draws, mu = rep(0, dim), Sigma = diag(dim)),
		prior = rep(1 / draws, draws),
		itemloc = model@Model$itemloc,
		full = FALSE,
		Etable = TRUE
	)$expected %>%
		log()
}

make_fulldata <- function(data){
	data <- as.matrix(data)
	wrong <- 1 - data
	right <- data
	colnames(wrong) <- glue::glue("Item.{1:ncol(data)}_1")
	colnames(right) <- glue::glue("Item.{1:ncol(data)}_2")
	out <- cbind(wrong, right)[, order(c(seq(ncol(wrong)), seq(ncol(right))))]

	# I don't totally get this but I think we just make NA = 0
	# So that they are disregarded from the likelihood calculation
	out[is.na(out)] <- 0
	out
}

evaluate_p <- function(p, train, test){
	list(
		train_acc = mean((p > 0.5) == train, na.rm = TRUE),
		test_acc = mean((p > 0.5) == test, na.rm = TRUE),
		train_rmse = sum((train - p)^2, na.rm = TRUE) / sum(!is.na(train)),
		test_rmse = sum((test - p)^2, na.rm = TRUE) / sum(!is.na(test))
	)
}

get_p <- function(model, fscores){
	n_items <- length(model@Data$K)

	1:n_items %>%
		map(~ probtrace(extract.item(model, .), fscores)[ , 2]) %>%
		do.call(cbind, .)
}
