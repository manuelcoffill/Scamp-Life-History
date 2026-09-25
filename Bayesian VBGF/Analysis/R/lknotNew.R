#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
#			Scamp
#		
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Authors: Zach Siders
#Incept: 06/08/2026
###-----------------------------------------------------
#		Initialization
###-----------------------------------------------------
	library(cmdstanr)
	library(truncnorm)
	library(latex2exp)
	library(bayestestR)
	library(truncnorm)

library(qs2)


path <- "C:\\Users\\manue\\OneDrive\\Desktop\\Siders VBGF - v2"
setwd(path)
###-----------------------------------------------------
#		Functions
###-----------------------------------------------------
	vbgm_optim <- function(theta, L0, ages, len){
		linf <- exp(theta[1])
		k <- exp(theta[2])
		sigma <- exp(theta[3])

		l_est <- linf-(linf-L0)*exp(-k*(ages))

		NLL = -1*sum(dnorm(len, l_est, sigma, log=TRUE))

		return(NLL)
	}
	col2rgbA<-function(color,transparency)
	{
		rgb(t(col2rgb(color))/255,alpha=transparency)
	}

	poly.den.partial <- function(x, adj=2, col.poly='black', prop=1, from=0,  to=1, xpart=c(0,1), ypart=c(0,1), alpha=NULL, full=TRUE, median=TRUE){
		d <- density(x, adj=adj, from=from, to=to)
		d$y2 <- d$y/max(d$y) * prop
		d$x2 <- scales::rescale(d$x, from=c(from,to), to=xpart)
		d$y3 <- scales::rescale(d$y2, from=c(0,1), to=ypart)
		d$x.fin <- d$x2
		d$y.fin <- d$y3
		if(!is.null(alpha)){
			d$x.fin[d$y2<alpha] <- NA
			d$y.fin[d$y2<alpha] <- NA
			d$x.fin <- na.omit(d$x.fin)
			d$y.fin <- na.omit(d$y.fin)
		}
		polygon(c(d$x.fin, rev(d$x.fin)), 
		        c(rep(min(ypart),length(d$x.fin)), rev(d$y.fin)), 
		        col=col.poly, border=FALSE, xpd=NA)
		if(full) lines(d$x2, d$y3, xpd=NA)
		if(median){
			med.find <- which.min(abs(median(x) - d$x))
			segments(x0 = d$x2[med.find],
			        x1 = d$x2[med.find],
			        y0 = min(abs(ypart)),
			        y1 = d$y3[med.find], xpd=NA)
		}
	}
	poly.den.partial.flip <- function(x, adj=2, col.poly='grey80', prop=1, from=0,  to=1, xpart=c(0,1), ypart=c(0,1), alpha=NULL, full=TRUE, median=TRUE){
		d <- density(x, adj=adj, from=from, to=to)
		d$y2 <- d$y/max(d$y) * prop
		d$x2 <- scales::rescale(d$x, from=c(from,to), to=xpart)
		d$y3 <- scales::rescale(d$y2, from=c(0,1), to=ypart)
		d$x.fin <- d$x2
		d$y.fin <- d$y3
		if(!is.null(alpha)){
			d$x.fin[d$y2<alpha] <- NA
			d$y.fin[d$y2<alpha] <- NA
			d$x.fin <- na.omit(d$x.fin)
			d$y.fin <- na.omit(d$y.fin)
		}

		polygon(x = c(rep(ypart[which.min(abs(ypart))],length(d$x.fin)), 
		              rev(d$y.fin)),
		        y = c(d$x.fin, rev(d$x.fin)), 
		        col=col.poly, border=FALSE, xpd=NA)
		if(full) lines(d$y3, d$x2, xpd=NA)
		if(median){
			med.find <- which.min(abs(median(x) - d$x))
			segments(x0 = min(abs(ypart)),
			        x1 = d$y3[med.find],
			        y0 = d$x2[med.find],
			        y1 = d$x2[med.find], xpd=NA)
		}
	}
	lines.den.partial.flip <- function(x, adj=2, prop=1, from=0, to=1, xpart=c(0,1), ypart=c(0,1),...){
		d <- density(x, adj=adj, from=from, to=to)
		d$y2 <- d$y/max(d$y) * prop
		d$x2 <- scales::rescale(d$x, from=c(from,to), to=xpart)
		d$y3 <- scales::rescale(d$y2, from=c(0,1), to=ypart)

		lines(x = d$y3, y = d$x2, ...)
		
	}
	lines.den.partial <- function(x, adj=2, from=0, to=1,prop=1,xpart=c(0,1), ypart=c(0,1), ...){
		d <- density(x, adj=adj, from=from,to=to)
		d$y2 <- d$y/max(d$y) * prop
		d$x2 <- scales::rescale(d$x, from=c(from,to), to=xpart)
		d$y2 <- scales::rescale(d$y2, from=c(0,1), to=ypart)
		
		lines(d$x2, d$y2, ...)
	}
	fig.lab <- function(i, xscale=0.05, yscale, cex=1.4){
		text(x = par('usr')[1] + abs(diff(par('usr')[1:2]))*xscale,
		     y = par('usr')[3] + abs(diff(par('usr')[3:4]))*yscale,
		     ifelse(is.numeric(i),LETTERS[i],i), xpd=NA, cex=cex)
	}
	q_paste <- function(x, alpha=0.1){
		v <- quantile(x, probs=c(0.5,(alpha/2),(1-alpha/2))) 
		v <- formatC(v,digits=3, format='f')
		v<- gsub(" ","", v)
		paste0(v[1]," (",v[2],"-",v[3],")")
	}
	q2_save <- function(model){
		try(model$draws(),silent=TRUE)
		try(model$sampler_diagnostics(), silent = TRUE)
		qs2::qs_save(model, file = paste0("./STAN_Outputs/",deparse(substitute(model)),".qs2"))
	}	
###-----------------------------------------------------
#		Data Read In
###-----------------------------------------------------
	scamp <- read.csv("./Data/CSV/Scamp.csv")
	#drop NAs
	scamp <- na.omit(scamp)
###-----------------------------------------------------
#		VBGF Optim
###-----------------------------------------------------
	vbgm.fit.tl <- optim(log(c(900,0.1,10)), 
	                     vbgm_optim, 
	                     ages=scamp$ConsensusFracAge, 
	                     len=scamp$ForkLengthMM, 
	                     L0 = 30)
	vbgm.par.tl <- exp(vbgm.fit.tl$par)
###-----------------------------------------------------
#		Priors
###-----------------------------------------------------
	pars_mu = c(vbgm.par.tl, 30)
	pars_mu[1] <- 900
	pars_sd = abs(pars_mu*0.8)

	labs <- c(expression(italic(L)[infinity]),
	          expression(italic(k)),
	          expression(sigma),
	          expression(italic(L)[0]))
	par(mfrow=c(2,2),mar=c(2.5,2.5,2,1))
	for(i in 1:4){
		d <- density(rtruncnorm(1000, a=0, mean=pars_mu[i], sd=pars_sd[i]), from = 0)
		d$y <- d$y/max(d$y)
		plot(d$x, d$y, type='l', xlab = "", ylab = "",
		     main = labs[i], las = 1)
		abline(v = pars_mu[i], lwd=2, lty = 2)
	}
###-----------------------------------------------------
#		Data (multireader - state space)
###-----------------------------------------------------
	lknot_dat <- list(N = nrow(scamp),
	                  Nread = 3,
	                  age = as.matrix(scamp[,grep('R[[:digit:]]{1}FracAge',colnames(scamp))]),
	                  l = scamp$ForkLengthMM,
	                  age_cons = scamp$ConsensusFracAge,
	                  nseq = length(seq(0,max(scamp[,grep('FracAge',colnames(scamp))]),by=0.1)),
	                  seq_ages = seq(0,max(scamp[,grep('FracAge',colnames(scamp))]),by=0.1),
	                  Linf_prior = c(pars_mu[1],pars_sd[1]),
	                  k_prior = c(pars_mu[2],pars_sd[2]),
	                  Lknot_prior = c(pars_mu[4],pars_sd[4]),
	                  sigma_prior = c(pars_mu[3],pars_sd[3])) 
	seq_ages <- lknot_dat$seq_ages

	dir.create("./STAN_Outputs")
###-----------------------------------------------------
#		Old version (Generalized T)
###-----------------------------------------------------
	lknot_init <- function(chain_id){
		v <- matrix(rnorm(8,mean=pars_mu, sd=abs(pars_mu*0.1)),nrow=4,ncol=2)
		list(Linf = v[1,],
			k = v[2,],
			Lknot = v[4,],
			sigma_proc = rexp(2,10),
			sigma_obs = v[3,1],
			true_age = rowMeans(lknot_dat$age))
	}

	vbgf_lknot <- cmdstan_model("./Analysis/STAN/base_lknot_multi_ss_cons.stan",
	                            compile_model_methods = TRUE)

	vbgf_lknot_fit <- vbgf_lknot$sample(data = lknot_dat,
	                                   chains=4,
	                                   iter_warmup=3000,
	                                   iter_sampling=1000,
	                                   refresh=1000,
	                                   init = lknot_init,
	                                   show_messages = TRUE,
	                                   show_exceptions = FALSE)
	q2_save(vbgf_lknot_fit)
	# vbgf_lknot_fit$save_output_files(dir = "./Outputs/",
	#                                  basename = "genT_freelknot_ln")
###-----------------------------------------------------
#		Old version (Generalized T, RR)
###-----------------------------------------------------

	lknot_nu_init <- function(chain_id){
		v <- matrix(rnorm(8,mean=pars_mu, sd=abs(pars_mu*0.1)),nrow=4,ncol=2)
		list(Linf = v[1,],
			k = v[2,],
			Lknot = v[4,],
			cv_proc = v[3,]/100*0.5,
			sigma_obs = v[3,1],
			nu_proc = rlnorm(2,log(50),0.1),
			true_age = rowMeans(lknot_dat$age))
	}

	vbgf_lknot_RR <- cmdstan_model("./Analysis/STAN/base_lknot_multi_ss_cons_RR.stan",
	                            compile_model_methods = TRUE)

	vbgf_lknot_RR_fit <- vbgf_lknot_RR$sample(data = lknot_dat,
	                                   chains=4,
	                                   iter_warmup=3000,
	                                   iter_sampling=1000,
	                                   refresh=1000,
	                                   init = lknot_nu_init,
	                                   show_messages = TRUE,
	                                   show_exceptions = FALSE)
	vbgf_lknot_RR_fit$summary(c("Linf","k","Lknot",'nu_proc',"cv_proc","sigma_obs",'tknot'))
	q2_save(vbgf_lknot_RR_fit)
	# vbgf_lknot_RR_fit$save_output_files(dir = "./Outputs/",
	#                                  basename = "genT_freelknot_RR")
###-----------------------------------------------------
#		Alternative (Individual Reader version)
###-----------------------------------------------------

	lknot_init_alt <- function(chain_id){
		v <- matrix(rnorm(8,mean=pars_mu, sd=abs(pars_mu*0.1)),nrow=4,ncol=2)
		list(Linf = v[1,],
			k = v[2,],
			Lknot = v[4,],
			sigma_proc = rexp(2,10),
			sigma_obs = rep(v[3,1],lknot_dat$Nread),
			true_age = rowMeans(lknot_dat$age))
	}
	vbgf_lknot_alt <- cmdstan_model("./Analysis/STAN/base_lknot_multi_ss_cons_alt2.stan",
	                                compile_model_methods = TRUE)

	vbgf_lknot_fit_alt <- vbgf_lknot_alt$sample(data = lknot_dat,
	                                   chains=4,
	                                   iter_warmup=3000,
	                                   iter_sampling=1000,
	                                   refresh=1000,
	                                   init = lknot_init_alt,
	                                   show_messages = TRUE,
	                                   show_exceptions = FALSE)
	q2_save(vbgf_lknot_fit_alt)
	# vbgf_lknot_fit_alt$save_output_files(dir = "./Outputs/",
	#                                  basename = "iR_freelknot_ln")
###-----------------------------------------------------
#		Alternative version (IRV, RR)
###-----------------------------------------------------

	lknot_nualt_init <- function(chain_id){
		v <- matrix(rnorm(8,mean=pars_mu, sd=abs(pars_mu*0.1)),nrow=4,ncol=2)
		list(Linf = v[1,],
			k = v[2,],
			Lknot = v[4,],
			cv_proc = v[3,]/100*0.5,
			sigma_obs = rep(v[3,1],lknot_dat$Nread),
			nu_proc = rlnorm(2,log(50),0.1),
			true_age = rowMeans(lknot_dat$age))
	}

	vbgf_lknot_RRalt <- cmdstan_model("./Analysis/STAN/base_lknot_multi_ss_cons_RRalt.stan",
	                            compile_model_methods = TRUE)

	vbgf_lknot_RRalt_fit <- vbgf_lknot_RRalt$sample(data = lknot_dat,
	                                   chains=4,
	                                   iter_warmup=3000,
	                                   iter_sampling=1000,
	                                   refresh=1000,
	                                   init = lknot_nualt_init,
	                                   show_messages = TRUE,
	                                   show_exceptions = FALSE)
	vbgf_lknot_RRalt_fit$summary(c("Linf","k","Lknot",'nu_proc',"cv_proc","sigma_obs",'tknot'))
	q2_save(vbgf_lknot_RRalt_fit)
	# vbgf_lknot_RRalt_fit$save_output_files(dir = "./Outputs/",
	#                                  basename = "iR_freelknot_RR")
###-----------------------------------------------------
#		Known Lknot version
###-----------------------------------------------------
	lknot_dat_kwnL0 <- lknot_dat
	lknot_dat_kwnL0$Lknot_prior <- c(2.5,0.25)

	vbgf_kwnlknot_fit <- vbgf_lknot$sample(data = lknot_dat_kwnL0,
	                                   chains=4,
	                                   iter_warmup=3000,
	                                   iter_sampling=1000,
	                                   refresh=1000,
	                                   init = lknot_init,
	                                   show_messages = TRUE,
	                                   show_exceptions = FALSE)
	q2_save(vbgf_kwnlknot_fit)
	# vbgf_kwnlknot_fit$save_output_files(dir = "./Outputs/",
	#                                  basename = "genT_kwnlknot_ln")
	vbgf_kwnlknot_fit_alt <- vbgf_lknot_alt$sample(data = lknot_dat_kwnL0,
	                                   chains=4,
	                                   iter_warmup=3000,
	                                   iter_sampling=1000,
	                                   refresh=1000,
	                                   init = lknot_init_alt,
	                                   show_messages = TRUE,
	                                   show_exceptions = FALSE)
	q2_save(vbgf_kwnlknot_fit_alt)
	# vbgf_kwnlknot_fit_alt$save_output_files(dir = "./Outputs/",
	#                                  basename = "iR_kwnlknot_ln")
	vbgf_kwnlknot_fit_RR <- vbgf_lknot_RR$sample(data = lknot_dat_kwnL0,
	                                   chains=4,
	                                   iter_warmup=3000,
	                                   iter_sampling=1000,
	                                   refresh=1000,
	                                   init = lknot_nu_init,
	                                   show_messages = TRUE,
	                                   show_exceptions = FALSE)
	q2_save(vbgf_kwnlknot_fit_RR)
	# vbgf_kwnlknot_fit_RR$save_output_files(dir = "./Outputs/",
	#                                  basename = "genT_kwnlknot_RR")
	vbgf_kwnlknot_fit_RRalt <- vbgf_lknot_RRalt$sample(data = lknot_dat_kwnL0,
	                                   chains=4,
	                                   iter_warmup=3000,
	                                   iter_sampling=1000,
	                                   refresh=1000,
	                                   init = lknot_nualt_init,
	                                   show_messages = TRUE,
	                                   show_exceptions = FALSE)
	q2_save(vbgf_kwnlknot_fit_RRalt)
	# vbgf_kwnlknot_fit_RRalt$save_output_files(dir = "./Outputs/",
	#                                  basename = "iR_kwnlknot_RR")
###-----------------------------------------------------
#		Q2 LOAD
###-----------------------------------------------------
	vbgf_lknot_fit <- qs2::qs_read("./STAN_Outputs/vbgf_lknot_fit")
	vbgf_lknot_RR_fit <- qs2::qs_read("./STAN_Outputs/vbgf_lknot_RR_fit")
	vbgf_lknot_fit_alt <- qs2::qs_read("./STAN_Outputs/vbgf_lknot_fit_alt")
	vbgf_lknot_RRalt_fit <- qs2::qs_read("./STAN_Outputs/vbgf_lknot_RRalt_fit")
	vbgf_kwnlknot_fit <- qs2::qs_read("./STAN_Outputs/vbgf_kwnlknot_fit")
	vbgf_kwnlknot_fit_alt <- qs2::qs_read("./STAN_Outputs/vbgf_kwnlknot_fit_alt")
	vbgf_kwnlknot_fit_alt <- qs2::qs_read("./STAN_Outputs/vbgf_kwnlknot_fit_alt")
	vbgf_kwnlknot_fit_RRalt <- qs2::qs_read("./STAN_Outputs/vbgf_kwnlknot_fit_RRalt")
###-----------------------------------------------------
#		(CSV WAY) LOAD INTO IF NEED BE
###-----------------------------------------------------
	vbgf_lknot_fit <- as_cmdstan_fit(list.files("./Outputs/",
	                                 pattern = "genT_freelknot_ln",
	                                 full.names = TRUE))
	vbgf_lknot_RR_fit <- as_cmdstan_fit(list.files("./Outputs/",
	                                    pattern = "genT_freelknot_RR",
	                                    full.names = TRUE))
	vbgf_lknot_fit_alt <- as_cmdstan_fit(list.files("./Outputs/",
	                                     pattern = "iR_freelknot_ln",
	                                     full.names = TRUE))
	vbgf_lknot_RRalt_fit <- as_cmdstan_fit(list.files("./Outputs/",
	                                       pattern = "iR_freelknot_RR",
	                                       full.names = TRUE))
	vbgf_kwnlknot_fit <- as_cmdstan_fit(list.files("./Outputs/",
	                                    pattern = "genT_kwnlknot_ln",
	                                    full.names = TRUE))
	vbgf_kwnlknot_fit_alt <- as_cmdstan_fit(list.files("./Outputs/",
	                                        pattern = "iR_kwnlknot_ln",
	                                        full.names = TRUE))
	vbgf_kwnlknot_fit_RR <- as_cmdstan_fit(list.files("./Outputs/",
	                                       pattern = "genT_kwnlknot_RR",
	                                       full.names = TRUE))
	vbgf_kwnlknot_fit_RRalt <- as_cmdstan_fit(list.files("./Outputs/",
	                                          pattern = "iR_kwnlknot_RR",
	                                          full.names = TRUE))
###-----------------------------------------------------
#		STAN extract (old version)
###-----------------------------------------------------
	pars <- c("Linf","k","Lknot","sigma_proc","sigma_obs",'tknot')
	ext_par <- as.data.frame(vbgf_lknot_fit$draws(pars, format="matrix"))
	ext_Lage <- as.data.frame(vbgf_lknot_fit$draws("Lage", format="matrix"))
	
	q_Lage <- apply(ext_Lage, 2, quantile, probs=c(0.025,0.5,0.975))
	write.csv(describe_posterior(ext_par, centrality='all'),
	          file = './Tables/desc_post_GenT_freeLzero_LN.csv',
	          row.names=FALSE)
	write.csv(vbgf_lknot_fit$summary(pars),file="./Tables/post_summ_GenT_freeLzero_LN.csv", row.names=FALSE)
###-----------------------------------------------------
#		STAN extract (old version, robust regression)
###-----------------------------------------------------
	pars2 <- c("Linf","k","Lknot","cv_proc","sigma_obs",'tknot','nu_proc')
	ext_rr_par <- as.data.frame(vbgf_lknot_RR_fit$draws(pars2, format="matrix"))
	ext_rr_Lage <- as.data.frame(vbgf_lknot_RR_fit$draws("Lage", format="matrix"))
	
	q_rr_Lage <- apply(ext_rr_Lage, 2, quantile, probs=c(0.025,0.5,0.975))
	write.csv(describe_posterior(ext_rr_par, centrality='all'),
	          file = './Tables/desc_post_GenT_freeLzero_RR.csv',
	          row.names=FALSE)
	write.csv(vbgf_lknot_RR_fit$summary(pars2),file="./Tables/post_summ_GenT_freeLzero_RR.csv", row.names=FALSE)
###-----------------------------------------------------
#		STAN extract (old version, known L0)
###-----------------------------------------------------
	pars <- c("Linf","k","Lknot","sigma_proc","sigma_obs",'tknot')
	ext_kwn_par <- as.data.frame(vbgf_kwnlknot_fit$draws(pars, format="matrix"))
	ext_kwn_Lage <- as.data.frame(vbgf_kwnlknot_fit$draws("Lage", format="matrix"))
	
	q_kwn_Lage <- apply(ext_kwn_Lage, 2, quantile, probs=c(0.025,0.5,0.975))
	write.csv(describe_posterior(ext_kwn_par, centrality='all'),
	          file = './Tables/desc_post_GenT_knownLzero_LN.csv',
	          row.names=FALSE)
	write.csv(vbgf_kwnlknot_fit$summary(pars),file="./Tables/post_summ_GenT_knownLzero_LN.csv", row.names=FALSE)
###-----------------------------------------------------
#		STAN extract (RR version, known L0)
###-----------------------------------------------------
	pars2 <- c("Linf","k","Lknot","cv_proc","sigma_obs",'tknot','nu_proc')
	ext_kwnRR_par <- as.data.frame(vbgf_kwnlknot_fit_RR$draws(pars2, format="matrix"))
	ext_kwnRR_Lage <- as.data.frame(vbgf_kwnlknot_fit_RR$draws("Lage", format="matrix"))
	
	q_kwnRR_Lage <- apply(ext_kwnRR_Lage, 2, quantile, probs=c(0.025,0.5,0.975))
	write.csv(describe_posterior(ext_kwnRR_par, centrality='all'),
	          file = './Tables/desc_post_GenT_knownLzero_RR.csv',
	          row.names=FALSE)
	write.csv(vbgf_kwnlknot_fit_RR$summary(pars2),file="./Tables/post_summ_GenT_knownLzero_RR.csv", row.names=FALSE)
###-----------------------------------------------------
#		STAN extract (new version)
###-----------------------------------------------------
	ext_alt_par <- as.data.frame(vbgf_lknot_fit_alt$draws(pars, format="matrix"))
	ext_alt_Lage <- as.data.frame(vbgf_lknot_fit_alt$draws("Lage", format="matrix"))
	
	q_alt_Lage <- apply(ext_alt_Lage, 2, quantile, probs=c(0.025,0.5,0.975))
	write.csv(describe_posterior(ext_alt_par, centrality='all'),
	          file = './Tables/desc_post_IndR_freeLzero_LN.csv',
	          row.names=FALSE)
	write.csv(vbgf_lknot_fit_alt$summary(pars),file="./Tables/post_summ_IndR_freeLzero_LN.csv", row.names=FALSE)
###-----------------------------------------------------
#		STAN extract (new version, RR)
###-----------------------------------------------------
	ext_altRR_par <- as.data.frame(vbgf_lknot_RRalt_fit$draws(pars2, format="matrix"))
	ext_altRR_Lage <- as.data.frame(vbgf_lknot_RRalt_fit$draws("Lage", format="matrix"))
	
	q_altRR_Lage <- apply(ext_altRR_Lage, 2, quantile, probs=c(0.025,0.5,0.975))
	write.csv(describe_posterior(ext_altRR_par, centrality='all'),
	          file = './Tables/desc_post_IndR_freeLzero_RR.csv',
	          row.names=FALSE)
	write.csv(vbgf_lknot_RRalt_fit$summary(pars2),file="./Tables/post_summ_IndR_freeLzero_RR.csv", row.names=FALSE)
###-----------------------------------------------------
#		STAN extract (new version, known L0)
###-----------------------------------------------------
	ext_kwnalt_par <- as.data.frame(vbgf_kwnlknot_fit_alt$draws(pars, format="matrix"))
	ext_kwnalt_Lage <- as.data.frame(vbgf_kwnlknot_fit_alt$draws("Lage", format="matrix"))
	
	q_kwnalt_Lage <- apply(ext_kwnalt_Lage, 2, quantile, probs=c(0.025,0.5,0.975))
	write.csv(describe_posterior(ext_kwnalt_par, centrality='all'),
	          file = './Tables/desc_post_IndR_knownLzero_LN.csv',
	          row.names=FALSE)
	write.csv(vbgf_kwnlknot_fit_alt$summary(pars),file="./Tables/post_summ_IndR_knownLzero_LN.csv", row.names=FALSE)
###-----------------------------------------------------
#		STAN extract (RRalt version, known L0)
###-----------------------------------------------------

	ext_kwnRRalt_par <- as.data.frame(vbgf_kwnlknot_fit_RRalt$draws(pars2, format="matrix"))
	ext_kwnRRalt_Lage <- as.data.frame(vbgf_kwnlknot_fit_RRalt$draws("Lage", format="matrix"))
	
	q_kwnRRalt_Lage <- apply(ext_kwnRRalt_Lage, 2, quantile, probs=c(0.025,0.5,0.975))
	write.csv(describe_posterior(ext_kwnRRalt_par, centrality='all'),
	          file = './Tables/desc_post_IndR_knownLzero_RR.csv',
	          row.names=FALSE)
	write.csv(vbgf_kwnlknot_fit_RRalt$summary(pars2),file="./Tables/post_summ_IndR_knownLzero_RR.csv", row.names=FALSE)
###-----------------------------------------------------
#		STAN extract - compare consensus with true age
###-----------------------------------------------------
	ext_Tage <- as.data.frame(vbgf_lknot_fit$draws("true_age", format="matrix"))
	q.Tage <- apply(ext_Tage, 2, quantile, probs=c(0.025,0.5,0.975))

	ext_Tage_alt <- as.data.frame(vbgf_lknot_fit_alt$draws("true_age", format="matrix"))
	q.Tage_alt <- apply(ext_Tage_alt, 2, quantile, probs=c(0.025,0.5,0.975))

	par(mfrow=c(1,2), mgp=c(2.5,1,0),
	    mar=c(4,4,1.5,1))
	#old version
	plot(lknot_dat$age_cons, q.Tage[2,],type='n',
	     las = 1, xlab = "Consensus Age (y)", ylab = expression(paste('Latent Age (',hat(italic(t))[i],")")), xpd=NA)
	segments(x0 =lknot_dat$age_cons, x1=lknot_dat$age_cons,
	         y0 = q.Tage[1,], y1 = q.Tage[3,])
	points(lknot_dat$age_cons, q.Tage[2,],pch=21,col='white',bg='black', cex=1.2)
	abline(a=0,b=1,col='red')
	title('Generalized T version')

	#new version
	plot(lknot_dat$age_cons, q.Tage_alt2[2,],type='n',xlab = "Consensus Age (y)", ylab = expression(paste('Latent Age (',hat(italic(t))[i],")")), xpd=NA, las=1)
	segments(x0 =lknot_dat$age_cons, x1=lknot_dat$age_cons,
	         y0 = q.Tage_alt2[1,], y1 = q.Tage_alt2[3,])
	points(lknot_dat$age_cons, q.Tage_alt2[2,],pch=21,col='white',bg='black', cex=1.2)
	abline(a=0,b=1,col='red')
	title('Individual Reader version')
###-----------------------------------------------------
#		STAN (loo)
###-----------------------------------------------------
	old_free_read <- vbgf_lknot_fit$loo(variables="log_like_read")
	new_free_read <- vbgf_lknot_fit_alt$loo(variables="log_like_read")
	old_kwn_read <- vbgf_kwnlknot_fit$loo(variables="log_like_read")
	new_kwn_read <- vbgf_kwnlknot_fit_alt$loo(variables="log_like_read")
	rr_free_read <- vbgf_lknot_RR_fit$loo(variables="log_like_read")
	rr_kwn_read <- vbgf_kwnlknot_fit_RR$loo(variables="log_like_read")
	rralt_free_read <- vbgf_lknot_RRalt_fit$loo(variables="log_like_read")
	rralt_kwn_read <- vbgf_kwnlknot_fit_RRalt$loo(variables="log_like_read")
	old_free_cons <- vbgf_lknot_fit$loo(variables="log_like_cons")
	new_free_cons <- vbgf_lknot_fit_alt$loo(variables="log_like_cons")
	old_kwn_cons <- vbgf_kwnlknot_fit$loo(variables="log_like_cons")
	new_kwn_cons <- vbgf_kwnlknot_fit_alt$loo(variables="log_like_cons")
	rr_free_cons <- vbgf_lknot_RR_fit$loo(variables="log_like_cons")
	rr_kwn_cons <- vbgf_kwnlknot_fit_RR$loo(variables="log_like_cons")
	rralt_free_cons <- vbgf_lknot_RRalt_fit$loo(variables="log_like_cons")
	rralt_kwn_cons <- vbgf_kwnlknot_fit_RRalt$loo(variables="log_like_cons")

	loo_read <- list(old_free = old_free_read,
	                 new_free = new_free_read,
	                 old_kwn = old_kwn_read,
	                 new_kwn = new_kwn_read,
	                 rr_free = rr_free_read,
	                 rr_kwn = rr_kwn_read,
	                 rrnew_free = rralt_free_read,
	                 rrnew_kwn = rralt_kwn_read)
	write.csv(loo::loo_compare(loo_read),file="./Tables/loo_compare_reads.csv")
	loo_cons <- list(old_free = old_free_cons,
	                 new_free = new_free_cons,
	                 old_kwn = old_kwn_cons,
	                 new_kwn = new_kwn_cons,
	                 rr_free = rr_free_cons,
	                 rr_kwn = rr_kwn_cons,
	                 rrnew_free = rralt_free_cons,
	                 rrnew_kwn = rralt_kwn_cons)
	write.csv(loo::loo_compare(loo_cons),file="./Tables/loo_compare_cons.csv")
###-----------------------------------------------------
#		Figures
###-----------------------------------------------------
	pdf("./Figures/vbgf_lknot_LN.pdf",width=9.546296,height=6.166667)
		par(mar=c(3,5,3,3.75),mfrow=c(1,1))
		read.col <- rcartocolor::carto_pal(4,'SunsetDark')[-1]
		#VBGM
			plot(seq_ages, rep(1,length(seq_ages)), 
			     xlim=c(-0.1,max(seq_ages)+5), ylim=c(0,950),xaxs='i',
			     type='n', las=1, xaxs='i', xaxt='n',
			     xlab="", ylab="", yaxs='i', bty='l')
			# polygon(x = c(seq_ages, rev(seq_ages)),
			#         y = c(q_Lage[1,], rev(q_Lage[3,])),
			#         col='gray80', border='gray60')
			axis(1)
		#Length distribution
			poly.den.partial.flip(lknot_dat$l, 
			                      from=min(lknot_dat$l), 
			                      to=max(lknot_dat$l),
			                      col.poly='grey90', 
			                      ypart=c(max(seq_ages)+0.05,
			                              max(seq_ages)+3.5),
			                      xpart=range(lknot_dat$l), adj=0.5,
			                      full=TRUE, median=FALSE)
		#Age Distribution
			poly.den.partial(lknot_dat$age_cons, from = min(lknot_dat$age_cons), to = max(lknot_dat$age_cons),
			                 col.poly='grey90', ypart=c(par('usr')[4], par('usr')[4]+75), 
			                 xpart=range(lknot_dat$age_cons), full=TRUE, median=FALSE, adj=0.5)
		#Size at birth
			poly.den.partial.flip(ext_alt_par$"Lknot[1]", 
			                      from=min(ext_alt_par$"Lknot[1]"), 
			                      to=max(ext_alt_par$"Lknot[1]"),
			                      col.poly=col2rgbA('gray50',0.7), 
			                      ypart=c(0,3.5), xpart=range(ext_alt_par$"Lknot[1]"), 
			                      alpha=0.1, full=TRUE, median=TRUE)
		#Median fits
			polygon(x = c(seq_ages, rev(seq_ages)),
			        y = c(q_kwnalt_Lage[1,1:length(seq_ages)],
			              q_kwnalt_Lage[3,length(seq_ages):1]),
			        col = col2rgbA('dodgerblue1',0.2),
			        border = NA)
			polygon(x = c(seq_ages, rev(seq_ages)),
			        y = c(q_alt_Lage[1,1:length(seq_ages)],
			              q_alt_Lage[3,length(seq_ages):1]),
			        col = col2rgbA('gray',0.2),
			        border = NA)
			# lines(seq_ages, q_Lage[2,1:length(seq_ages)], col='black')
			# lines(seq_ages, q_Lage[2,(length(seq_ages)+1):ncol(q_Lage)], col='black')
			lines(seq_ages, q_alt_Lage[2,1:length(seq_ages)], col='black',lwd = 3)
			lines(seq_ages, q_alt_Lage[1,1:length(seq_ages)], col='black',lty=2)
			lines(seq_ages, q_alt_Lage[3,1:length(seq_ages)], col='black',lty=2)
			# lines(seq_ages, q_alt_Lage[2,(length(seq_ages)+1):ncol(q_Lage)], col='gray50', lty=3,lwd = 3)

			# lines(seq_ages, q_kwn_Lage[2,1:length(seq_ages)], col='red')
			# lines(seq_ages, q_kwn_Lage[2,(length(seq_ages)+1):ncol(q_Lage)], col='red', lty=3)
			lines(seq_ages, q_kwnalt_Lage[2,1:length(seq_ages)], col='dodgerblue3',lwd = 3)
			# lines(seq_ages, q_kwnalt_Lage[2,(length(seq_ages)+1):ncol(q_Lage)], col='dodgerblue1', lty=3,lwd = 3)
			
			lines(seq_ages, q_kwnalt_Lage[3,1:length(seq_ages)], col='dodgerblue3',lty=2)
			lines(seq_ages, q_kwnalt_Lage[1,1:length(seq_ages)], col='dodgerblue3', lty=2)
		#Linf
			poly.den.partial.flip(ext_alt_par$"Linf[1]", from=min(ext_alt_par$"Linf[1]"), 
			                      to=max(ext_alt_par$"Linf[1]"),
		                      col.poly=col2rgbA('gray50',0.7), 
		                      ypart=c(max(seq_ages)+0.05, max(seq_ages)+4),
		                      xpart=range(ext_alt_par$"Linf[1]"), 
		                      alpha=0.1, full=TRUE, median=TRUE)
			poly.den.partial.flip(ext_kwnalt_par$"Linf[1]", from=min(ext_kwnalt_par$"Linf[1]"), 
			                      to=max(ext_kwnalt_par$"Linf[1]"),
		                      col.poly=col2rgbA('dodgerblue3',0.7), 
		                      ypart=c(max(seq_ages)+0.05, max(seq_ages)+4),
		                      xpart=range(ext_kwnalt_par$"Linf[1]"), 
		                      alpha=0.1, full=TRUE, median=TRUE)
			# points(y=gen_nester/10,
			#        x=jitter(rep(max(seq_ages)+0.75, length(gen_nester))),
			#        pch=18, col='palegreen4', xpd=NA)
			Linf.prior <- rtruncnorm(1000, a=0, mean=lknot_dat$Linf_prior[1], sd = lknot_dat$Linf_prior[2])
			lines.den.partial.flip(Linf.prior, from = min(lknot_dat$l), 
			                       to = max(lknot_dat$l), 
			                       ypart = c(max(seq_ages)+0.75, max(seq_ages)+4), 
			                       xpart = range(lknot_dat$l),
			                       col='black', lty=1, lwd=2)
		#k
			k.r <- range(c(ext_alt_par$"k[1]",ext_kwnalt_par$"k[1]"))
			poly.den.partial(ext_alt_par$"k[1]", from=min(k.r), to = max(k.r),
			                 col.poly = col2rgbA('gray50',0.7), ypart=c(100,180),
			                 xpart = c(10,25), alpha=0.1, full=TRUE, median=TRUE)
			poly.den.partial(ext_kwnalt_par$"k[1]", from=min(k.r), to = max(k.r),
			                 col.poly = col2rgbA('dodgerblue3',0.7), ypart=c(100,180),
			                 xpart = c(10,25), alpha=0.1, full=TRUE, median=TRUE)
			k.sim <- rtruncnorm(10000, a = 0, b = Inf, lknot_dat$k_prior[1], lknot_dat$k_prior[2])
			lines.den.partial(k.sim, from=min(k.r), to = max(k.r),
			                 ypart=c(100,180), xpart = c(10,25),
			                 col = 'black', lwd=2, lty=1)
			k.axis <- pretty(k.r)
			axis(1, at=seq(10,25, length.out=length(k.axis)), 
			     labels = NA, pos=100, tck=-0.01)
			text(x = seq(10,25, length.out=length(k.axis)), 
			     y = rep(80, length(k.axis)),
			     k.axis, adj=c(1,0.5), srt=45, cex=0.8)
		#sigma
			sigma.r <- range(c(ext_alt_par$"sigma_proc[1]",ext_kwnalt_par$"sigma_proc[1]"))
			poly.den.partial(ext_alt_par$"sigma_proc[1]", from=min(sigma.r), to = max(sigma.r),
			                 col.poly = col2rgbA('gray50',0.7), ypart=c(260,340),
			                 xpart = c(10,25), alpha=0.1, full=TRUE, median=TRUE)
			poly.den.partial(ext_kwnalt_par$"sigma_proc[1]", from=min(sigma.r), to = max(sigma.r),
			                 col.poly = col2rgbA('dodgerblue3',0.7), ypart=c(260,340),
			                 xpart = c(10,25), alpha=0.1, full=TRUE, median=TRUE)
			sigma.sim <- rtruncnorm(10000, a = 0, b = Inf, lknot_dat$sigma_prior[1], lknot_dat$sigma_prior[2])
			lines.den.partial(sigma.sim, from=min(sigma.r), to = max(sigma.r),
			                 ypart=c(260,340), xpart = c(10,25),
			                 col = 'black', lwd=2, lty=1)
			sigma.axis <- pretty(sigma.r)
			axis(1, at=seq(10,25, length.out=length(sigma.axis)), 
			     labels = NA, pos=260, tck=-0.01)
			text(x = seq(10,25, length.out=length(sigma.axis)), 
			     y = rep(240, length(sigma.axis)),
			     sigma.axis, adj=c(1,0.5), srt=45, cex=0.8)
		#tknot
			poly.den.partial(ext_alt_par$'tknot[1]', from=min(ext_alt_par$'tknot[1]'), to = max(ext_alt_par$'tknot[1]'),
			                 col.poly = col2rgbA('gray50',0.7), ypart=c(340,420),
			                 xpart = c(30,40), alpha=0.1, full=TRUE, median=TRUE)
			t0.axis <- pretty(range(ext_alt_par$'tknot[1]'))
			axis(1, at=seq(30,40, length.out=length(t0.axis)), 
			     labels = NA, pos=340, tck=-0.01)
			text(x = seq(30,40, length.out=length(t0.axis)), 
			     y = rep(320, length(t0.axis)),
			     t0.axis, adj=c(1,0.5), srt=45, cex=0.8)
		#sigma obs
			sigmao.r <- range(c(ext_alt_par$"sigma_obs[1]",ext_alt_par$"sigma_obs[2]",ext_alt_par$"sigma_obs[2]"))
			poly.den.partial(ext_alt_par$"sigma_obs[1]", from=min(sigmao.r), to = max(sigmao.r),
			                 col.poly = col2rgbA(read.col[1],0.7), ypart=c(100,180),
			                 xpart = c(30,45), alpha=0.1, full=TRUE, median=TRUE)
			poly.den.partial(ext_alt_par$"sigma_obs[2]", from=min(sigmao.r), to = max(sigmao.r),
			                 col.poly = col2rgbA(read.col[2],0.7), ypart=c(100,180),
			                 xpart = c(30,45), alpha=0.1, full=TRUE, median=TRUE)
			poly.den.partial(ext_alt_par$"sigma_obs[3]", from=min(sigmao.r), to = max(sigmao.r),
			                 col.poly = col2rgbA(read.col[3],0.7), ypart=c(100,180),
			                 xpart = c(30,45), alpha=0.1, full=TRUE, median=TRUE)
			
			sigmao.sim <- rtruncnorm(10000, a = 0, b = Inf, lknot_dat$sigma_prior[1], lknot_dat$sigma_prior[2])
			lines.den.partial(sigmao.sim, from=min(sigmao.r), to = max(sigmao.r),
			                 ypart=c(100,180), xpart = c(30,45),
			                 col = 'black', lwd=2, lty=1)
			sigmao.axis <- pretty(sigmao.r)
			axis(1, at=seq(30,45, length.out=length(sigmao.axis)), 
			     labels = NA, pos=100, tck=-0.01)
			text(x = seq(30,45, length.out=length(sigmao.axis)), 
			     y = rep(80, length(sigmao.axis)),
			     sigmao.axis, adj=c(1,0.5), srt=45, cex=0.8)
		#points
			for(i in 1:nrow(lknot_dat$age)){
				segments(x0 = min(lknot_dat$age[i,]),
				         x1 = max(lknot_dat$age[i,]),
				        y0 = lknot_dat$l[i],y1 = lknot_dat$l[i],
				        col='gray')
				points(lknot_dat$age[i,], rep(lknot_dat$l[i],lknot_dat$Nread), pch=21,
				       col = 'white', bg = read.col)
			}
		# axis(2, las=1)
		mtext("Fork Length (mm)", side=2, line=3)
		# axis(1)
		mtext("Age (yr.)", side=1, line=2)

		text(x=c(44.21253,27.52791,33.65369),
		     y=c(143.5754,143.5754,143.5754),
		     c(expression(sigma["obs,R1"]),
		       expression(sigma["obs,R2"]),
		       expression(sigma["obs,R3"])),
		     cex = 1.2, xpd=NA, adj=c(0,0.5))
		text(x=c(3.5,9.568459,23.162650,9.568459,23.162650,30.12294,46,46),
		     y=c(146.7847,141.4791,141.4791,306.8801,306.8801,398.9479,743.3436,663.6855),
		     labels=c(expression(italic(L)[0][',free']),expression(italic(k)['free']),expression(italic(k)['known']),expression(italic(sigma)['free']),expression(italic(sigma)['known']),expression(italic(t)[0][',free']),expression(italic(L)[infinity][',free']),expression(italic(L)[infinity][',known'])), cex=1.2, xpd=NA, adj = c(0,0.5))
		legend('topleft',legend=c('R1','R2','R3','Free L0','Known L0'), lty = c(NA,NA,NA,1,1), pch = c(rep(16,3),NA,NA), lwd = 3, col = c(read.col,'black','dodgerblue3'), cex = 0.8, inset = c(0.025,0.025), xpd=NA, pt.cex = 1.2)
	dev.off()
	
	
	
	
	
	
	
###-----------------------------------------------------
#		Figures (RR)
###-----------------------------------------------------
	## modification to show the 90% CRI for length-at-age rather than 95%
	q_altRR_Lagez <- apply(ext_altRR_Lage, 2, quantile, probs=c(0.05,0.5,0.95))
	q_kwnRRalt_Lagez <- apply(ext_kwnRRalt_Lage, 2, quantile, probs=c(0.05,0.5,0.95))
	
	
	
	#	pdf("./Figures/vbgf_lknot_RR.pdf",width=9.546296,height=6.166667)
	
	png(
	  filename = "./Figures/vbgf_lknot_RR_90CRI.png",
	  width = 10,
	  height = 8,
	  units = "in",
	  res = 600
	)
	
	
	
	
	
	
	
	
	
	
	
	# pdf("./Figures/vbgf_lknot_RR.pdf",width=9.546296,height=6.166667)
		
	
	
	par(mar=c(3,5,3,3.75),mfrow=c(1,1))
		read.col <- rcartocolor::carto_pal(4,'SunsetDark')[-1]
		#VBGM
			plot(seq_ages, rep(1,length(seq_ages)), 
			     xlim=c(-0.1,max(seq_ages)+5), ylim=c(0,950),xaxs='i',
			     type='n', las=1, xaxs='i', xaxt='n',
			     xlab="", ylab="", yaxs='i', bty='l')
			# polygon(x = c(seq_ages, rev(seq_ages)),
			#         y = c(q_Lage[1,], rev(q_Lage[3,])),
			#         col='gray80', border='gray60')
			axis(1)
		#Length distribution
			poly.den.partial.flip(lknot_dat$l, 
			                      from=min(lknot_dat$l), 
			                      to=max(lknot_dat$l),
			                      col.poly='grey90', 
			                      ypart=c(max(seq_ages)+0.05,
			                              max(seq_ages)+3.5),
			                      xpart=range(lknot_dat$l), adj=0.5,
			                      full=TRUE, median=FALSE)
		#Age Distribution
			poly.den.partial(lknot_dat$age_cons, from = min(lknot_dat$age_cons), to = max(lknot_dat$age_cons),
			                 col.poly='grey90', ypart=c(par('usr')[4], par('usr')[4]+75), 
			                 xpart=range(lknot_dat$age_cons), full=TRUE, median=FALSE, adj=0.5)
		#Size at birth
			poly.den.partial.flip(ext_altRR_par$"Lknot[1]", 
			                      from=min(ext_altRR_par$"Lknot[1]"), 
			                      to=max(ext_altRR_par$"Lknot[1]"),
			                      col.poly=col2rgbA('gray50',0.7), 
			                      ypart=c(0,3.5), xpart=range(ext_altRR_par$"Lknot[1]"), 
			                      alpha=0.1, full=TRUE, median=TRUE)
		#Median fits
			polygon(x = c(seq_ages, rev(seq_ages)),
			        y = c(q_kwnRRalt_Lagez[1,1:length(seq_ages)],
			              q_kwnRRalt_Lagez[3,length(seq_ages):1]),
			        col = col2rgbA('dodgerblue1',0.2),
			        border = NA)
			polygon(x = c(seq_ages, rev(seq_ages)),
			        y = c(q_altRR_Lagez[1,1:length(seq_ages)],
			              q_altRR_Lagez[3,length(seq_ages):1]),
			        col = col2rgbA('gray',0.2),
			        border = NA)
			# lines(seq_ages, q_Lage[2,1:length(seq_ages)], col='black')
			# lines(seq_ages, q_Lage[2,(length(seq_ages)+1):ncol(q_Lage)], col='black')
			lines(seq_ages, q_altRR_Lage[2,1:length(seq_ages)], col='black',lwd = 3)
			lines(seq_ages, q_altRR_Lagez[1,1:length(seq_ages)], col='black',lty=2)
			lines(seq_ages, q_altRR_Lagez[3,1:length(seq_ages)], col='black',lty=2)
			# lines(seq_ages, q_alt_Lage[2,(length(seq_ages)+1):ncol(q_Lage)], col='gray50', lty=3,lwd = 3)

			# lines(seq_ages, q_kwn_Lage[2,1:length(seq_ages)], col='red')
			# lines(seq_ages, q_kwn_Lage[2,(length(seq_ages)+1):ncol(q_Lage)], col='red', lty=3)
			lines(seq_ages, q_kwnRRalt_Lage[2,1:length(seq_ages)], col='dodgerblue3',lwd = 3)
			# lines(seq_ages, q_kwnRRalt_Lage[2,(length(seq_ages)+1):ncol(q_Lage)], col='dodgerblue1', lty=3,lwd = 3)
			
			lines(seq_ages, q_kwnRRalt_Lagez[3,1:length(seq_ages)], col='dodgerblue3',lty=2)
			lines(seq_ages, q_kwnRRalt_Lagez[1,1:length(seq_ages)], col='dodgerblue3', lty=2)
		#Linf
			poly.den.partial.flip(ext_altRR_par$"Linf[1]", from=min(ext_altRR_par$"Linf[1]"), 
			                      to=max(ext_altRR_par$"Linf[1]"),
		                      col.poly=col2rgbA('gray50',0.7), 
		                      ypart=c(max(seq_ages)+0.05, max(seq_ages)+4),
		                      xpart=range(ext_altRR_par$"Linf[1]"), 
		                      alpha=0.1, full=TRUE, median=TRUE)
			poly.den.partial.flip(ext_kwnRRalt_par$"Linf[1]", from=min(ext_kwnRRalt_par$"Linf[1]"), 
			                      to=max(ext_kwnRRalt_par$"Linf[1]"),
		                      col.poly=col2rgbA('dodgerblue3',0.7), 
		                      ypart=c(max(seq_ages)+0.05, max(seq_ages)+4),
		                      xpart=range(ext_kwnRRalt_par$"Linf[1]"), 
		                      alpha=0.1, full=TRUE, median=TRUE)
			# points(y=gen_nester/10,
			#        x=jitter(rep(max(seq_ages)+0.75, length(gen_nester))),
			#        pch=18, col='palegreen4', xpd=NA)
			Linf.prior <- rtruncnorm(1000, a=0, mean=lknot_dat$Linf_prior[1], sd = lknot_dat$Linf_prior[2])
			lines.den.partial.flip(Linf.prior, from = min(lknot_dat$l), 
			                       to = max(lknot_dat$l), 
			                       ypart = c(max(seq_ages)+0.75, max(seq_ages)+4), 
			                       xpart = range(lknot_dat$l),
			                       col='black', lty=1, lwd=2)
		#k
			k.r <- range(c(ext_altRR_par$"k[1]",ext_kwnRRalt_par$"k[1]"))
			poly.den.partial(ext_altRR_par$"k[1]", from=min(k.r), to = max(k.r),
			                 col.poly = col2rgbA('gray50',0.7), ypart=c(100,180),
			                 xpart = c(10,25), alpha=0.1, full=TRUE, median=TRUE)
			poly.den.partial(ext_kwnRRalt_par$"k[1]", from=min(k.r), to = max(k.r),
			                 col.poly = col2rgbA('dodgerblue3',0.7), ypart=c(100,180),
			                 xpart = c(10,25), alpha=0.1, full=TRUE, median=TRUE)
			k.sim <- rtruncnorm(10000, a = 0, b = Inf, lknot_dat$k_prior[1], lknot_dat$k_prior[2])
			lines.den.partial(k.sim, from=min(k.r), to = max(k.r),
			                 ypart=c(100,180), xpart = c(10,25),
			                 col = 'black', lwd=2, lty=1)
			k.axis <- pretty(k.r)
			axis(1, at=seq(10,25, length.out=length(k.axis)), 
			     labels = NA, pos=100, tck=-0.01)
			text(x = seq(10,25, length.out=length(k.axis)), 
			     y = rep(80, length(k.axis)),
			     k.axis, adj=c(1,0.5), srt=45, cex=0.8)
		#sigma
			sigma.r <- range(c(ext_altRR_par$"cv_proc[1]",ext_kwnRRalt_par$"cv_proc[1]"))
			poly.den.partial(ext_altRR_par$"cv_proc[1]", from=min(sigma.r), to = max(sigma.r),
			                 col.poly = col2rgbA('gray50',0.7), ypart=c(260,340),
			                 xpart = c(10,25), alpha=0.1, full=TRUE, median=TRUE)
			poly.den.partial(ext_kwnRRalt_par$"cv_proc[1]", from=min(sigma.r), to = max(sigma.r),
			                 col.poly = col2rgbA('dodgerblue3',0.7), ypart=c(260,340),
			                 xpart = c(10,25), alpha=0.1, full=TRUE, median=TRUE)
			sigma.sim <- rtruncnorm(10000, a = 0, b = Inf, lknot_dat$sigma_prior[1], lknot_dat$sigma_prior[2])
			lines.den.partial(sigma.sim, from=min(sigma.r), to = max(sigma.r),
			                 ypart=c(260,340), xpart = c(10,25),
			                 col = 'black', lwd=2, lty=1)
			sigma.axis <- pretty(sigma.r*100)
			axis(1, at=seq(10,25, length.out=length(sigma.axis)), 
			     labels = NA, pos=260, tck=-0.01)
			text(x = seq(10,25, length.out=length(sigma.axis)), 
			     y = rep(240, length(sigma.axis)),
			     sigma.axis, adj=c(1,0.5), srt=45, cex=0.8)
		#tknot
			poly.den.partial(ext_altRR_par$'tknot[1]', from=min(ext_altRR_par$'tknot[1]'), to = max(ext_altRR_par$'tknot[1]'),
			                 col.poly = col2rgbA('gray50',0.7), ypart=c(340,420),
			                 xpart = c(30,40), alpha=0.1, full=TRUE, median=TRUE)
			t0.axis <- pretty(range(ext_altRR_par$'tknot[1]'))
			axis(1, at=seq(30,40, length.out=length(t0.axis)), 
			     labels = NA, pos=340, tck=-0.01)
			text(x = seq(30,40, length.out=length(t0.axis)), 
			     y = rep(320, length(t0.axis)),
			     t0.axis, adj=c(1,0.5), srt=45, cex=0.8)
		#sigma obs
			sigmao.r <- range(c(ext_altRR_par$"sigma_obs[1]",ext_altRR_par$"sigma_obs[2]",ext_altRR_par$"sigma_obs[2]"))
			poly.den.partial(ext_altRR_par$"sigma_obs[1]", from=min(sigmao.r), to = max(sigmao.r),
			                 col.poly = col2rgbA(read.col[1],0.7), ypart=c(100,180),
			                 xpart = c(30,45), alpha=0.1, full=TRUE, median=TRUE)
			poly.den.partial(ext_altRR_par$"sigma_obs[2]", from=min(sigmao.r), to = max(sigmao.r),
			                 col.poly = col2rgbA(read.col[2],0.7), ypart=c(100,180),
			                 xpart = c(30,45), alpha=0.1, full=TRUE, median=TRUE)
			poly.den.partial(ext_altRR_par$"sigma_obs[3]", from=min(sigmao.r), to = max(sigmao.r),
			                 col.poly = col2rgbA(read.col[3],0.7), ypart=c(100,180),
			                 xpart = c(30,45), alpha=0.1, full=TRUE, median=TRUE)
			
			sigmao.sim <- rtruncnorm(10000, a = 0, b = Inf, lknot_dat$sigma_prior[1], lknot_dat$sigma_prior[2])
			lines.den.partial(sigmao.sim, from=min(sigmao.r), to = max(sigmao.r),
			                 ypart=c(100,180), xpart = c(30,45),
			                 col = 'black', lwd=2, lty=1)
			sigmao.axis <- pretty(sigmao.r)
			axis(1, at=seq(30,45, length.out=length(sigmao.axis)), 
			     labels = NA, pos=100, tck=-0.01)
			text(x = seq(30,45, length.out=length(sigmao.axis)), 
			     y = rep(80, length(sigmao.axis)),
			     sigmao.axis, adj=c(1,0.5), srt=45, cex=0.8)
		#points
			for(i in 1:nrow(lknot_dat$age)){
				segments(x0 = min(lknot_dat$age[i,]),
				         x1 = max(lknot_dat$age[i,]),
				        y0 = lknot_dat$l[i],y1 = lknot_dat$l[i],
				        col='gray')
				points(lknot_dat$age[i,], rep(lknot_dat$l[i],lknot_dat$Nread), pch=21,
				       col = 'white', bg = read.col)
			}
		# axis(2, las=1)
		mtext("Fork length (mm)", side=2, line=3)
		# axis(1)
		mtext("Age (y)", side=1, line=2)

		text(x=c(44.21253,27.52791,33.65369),
		     y=c(143.5754,143.5754,143.5754),
		     c(expression(sigma["obs,R1"]),
		       expression(sigma["obs,R2"]),
		       expression(sigma["obs,R3"])),
		     cex = 1.2, xpd=NA, adj=c(0,0.5))
		text(x=c(3.5,9.568459,23.162650,9.568459,23.162650,30.12294,46,46),
		     y=c(146.7847,141.4791,141.4791,306.8801,306.8801,398.9479,743.3436,663.6855),
		     labels=c(expression(italic(L)[0][',free']),expression(italic(k)['free']),expression(italic(k)['known']),expression('CV'['free']),expression('CV'['known']),expression(italic(t)[0][',free']),expression(italic(L)[infinity][',free']),expression(italic(L)[infinity][',known'])), cex=1.2, xpd=NA, adj = c(0,0.5))
		legend('topleft',legend=c('R1','R2','R3','Free L0','Known L0'), lty = c(NA,NA,NA,1,1), pch = c(rep(16,3),NA,NA), lwd = 3, col = c(read.col,'black','dodgerblue3'), cex = 0.8, inset = c(0.025,0.025), xpd=NA, pt.cex = 1.2)
	dev.off()
	
	
	
	
	
	
	################################################################################
	## add on by MCR
	## saving LAA preds
	## UPDATE now saves 90% CRI rahter than 95
	
	agepreds <- data.frame(
	  "Age" = seq_ages,
	  "Known low CRI" = q_kwnRRalt_Lagez[1,],
	  "Known med" = q_kwnRRalt_Lagez[2,],
	  "Known high CRI" = q_kwnRRalt_Lagez[3,],
	  "Free low CRI" = q_altRR_Lagez[1,],
	  "Free med CRI" = q_altRR_Lagez[2,],
	  "Free high CRI" = q_altRR_Lagez[3,]
	)
	
	agepreds
	nrow(agepreds)
	
	# Load package
	library(openxlsx)
	# Save dataframe as an Excel file
	write.xlsx(agepreds, "LengthAtAgePreds90CRI.xlsx", overwrite = TRUE)
	
	
	
	###########################################################
	####save posterior draws n=4000
	#### full distribution distributions (n=4000)
	
	# free L0
	postparamsfree <- ext_altRR_par
	head(postparamsfree)
	nrow(postparamsfree)
	
	#known L0
	postparamsknown <- ext_kwnRRalt_par 
	head(postparamsknown)
	
	## save these
	library(openxlsx)
	
	write.xlsx(postparamsfree, "PosteriorsFreeL0.xlsx", overwrite = TRUE)
	write.xlsx(postparamsknown, "PosteriorsKnownL0.xlsx", overwrite = TRUE)
	
	
	
	
	
	
	
	
	