library(ggplot2)
library(viridis)
library(readxl)
library(here)
library(dplyr)
library(tidyverse)
library(MuMIn)
library(viridis)

set.seed(123)

dat <- read_excel(here("Scamp staging.xlsx"))
nrow(dat)
table(dat$HistoSex)

#####dat################################################################################
################################################################################
# male and female df
dat <- dat[dat$HistoSex == "F" | dat$HistoSex == "M",]
nrow(dat)
summary(dat$FL)


# add a binary column for maturity
dat <- dat %>%
  mutate(
    Male = case_when(
      HistoSex == "M" ~ 1,
      HistoSex == "F" ~ 0,
      TRUE ~ NA_real_
    )
  )

unique(dat$HistoSex)
table(dat$Male, useNA = "ifany")

## check lengths
ggplot(dat, aes(FL)) +
  geom_histogram()


### ready for transition analysis :)



################################################################################
################################################################################

############################
# 3. FIT GLMs (MULTIPLE LINKS)
############################
dat$Length <- dat$FL
nrow(dat)
table(dat$Male)


m_logit   <- glm(Male ~ Length, family = binomial(link = "logit"),   data = dat)
m_probit  <- glm(Male ~ Length, family = binomial(link = "probit"),  data = dat)
m_cloglog <- glm(Male ~ Length, family = binomial(link = "cloglog"), data = dat)
m_cauchit <- glm(Male ~ Length, family = binomial(link = "cauchit"), data = dat)

model_set <- list(
  logit   = m_logit,
  probit  = m_probit,
  cloglog = m_cloglog,
  cauchit = m_cauchit
)


############################
# 4. MODEL SELECTION (AIC)
############################
print(model.sel(model_set))
# save this!



############################
# 5. L50 ESTIMATION FUNCTION
############################
# Handles link-specific math correctly

get_L50 <- function(model, link) {
  coefs <- coef(model)
  b0 <- coefs[1]
  b1 <- coefs[2]
  
  if (link %in% c("logit", "probit", "cauchit")) {
    L50 <- -b0 / b1
  } else if (link == "cloglog") {
    L50 <- (log(-log(1 - 0.5)) - b0) / b1
  } else {
    stop("Unknown link")
  }
  
  return(L50)
}

# Point estimates
L50_estimates <- mapply(
  get_L50,
  model = model_set,
  link = names(model_set)
)

print(L50_estimates)





## quick model summary list
model_summary_full <- map_df(names(model_set), function(link) {
  
  model <- model_set[[link]]
  
  tibble(
    Model = link,
    Intercept = coef(model)[1],
    Slope = coef(model)[2],
    AIC = AIC(model),
    Deviance = deviance(model),
    L50 = get_L50(model, link)
  )
})
model_summary_full


############################
# 6. DELTA METHOD UNCERTAINTY
############################

### removed. bootstrapping is more robust, although slower


############################
# 7. BOOTSTRAP UNCERTAINTY (10,000 DRAWS)
############################
# This resamples individuals WITH replacement,
# refits models, and recalculates L50 each time

n_boot <- 5000

bootstrap_L50 <- function(data, link) {
  
  # Resample rows
  boot_dat <- data %>%
    slice_sample(n = nrow(data), replace = TRUE)
  
  # Fit model
  model <- glm(Male ~ Length,
               family = binomial(link = link),
               data = boot_dat)
  
  # Extract L50
  return(get_L50(model, link))
}

# Run bootstrap for each link
boot_results <- lapply(names(model_set), function(link) {
  
  cat("Running bootstrap for:", link, "\n")
  
  replicate(n_boot, bootstrap_L50(dat, link))
})

names(boot_results) <- names(model_set)


############################
# 8. SUMMARIZE BOOTSTRAP RESULTS
############################
boot_summary <- map_df(names(boot_results), function(link) {
  
  vals <- boot_results[[link]]
  
  tibble(
    Model = link,
    L50_median = median(vals, na.rm = TRUE),
    L50_mean = mean(vals, na.rm = TRUE),
    L50_sd   = sd(vals, na.rm = TRUE),
    Lower_2.5 = quantile(vals, 0.025, na.rm = TRUE),
    Upper_97.5 = quantile(vals, 0.975, na.rm = TRUE)
  )
})

print(boot_summary)
# save this!

############################
# 9. GENERATE PREDICTIONS
############################
newdat <- tibble(
  Length = seq(0,850, length.out = 850)
)

preds <- map_df(names(model_set), function(name) {
  model <- model_set[[name]]
  
  newdat %>%
    mutate(
      fit = predict(model, newdata = newdat, type = "response"),
      Model = name
    )
})


############################
# 10. PLOT OGIVES
############################
ggplot(dat, aes(Length, Male)) +
  geom_jitter(height = 0.05, alpha = 0.5) +
  geom_line(data = preds, aes(y = fit, color = Model), linewidth = 1.2) +
  scale_color_viridis(discrete=T, option="D") +
  theme_bw() +
  labs(
    x = "Fork length (mm)",
    y = "Probability Male"
  ) +
  coord_cartesian(xlim = c(0, 850))



############################
# 11. PLOT WITH BOOTSTRAP L50 WITH MEAN
############################
ggplot(dat, aes(Length, Male)) +
  geom_jitter(height = 0.03, alpha = 0.5, size=2) +
  geom_line(data = preds, aes(y = fit, color = Model), linewidth = 1.5) +
  scale_color_viridis(discrete=T, option="D",
                      labels = tools::toTitleCase,
                      name = "Link") +
  geom_vline(
    data = boot_summary,
    aes(xintercept = L50_mean, color = Model),
    linetype = "dashed",
    linewidth=1
  ) +
  theme_bw() +
  labs(
    x = "Fork length (mm)",
    y = "Probability Male",
  ) +
  coord_cartesian(xlim = c(0, 850)) +
  theme(
    legend.position = c(0.2, 0.8),       # Position inside plot
    legend.background = element_rect(    # Add background for readability
      fill = "white", 
      color = "black"
    ),
    #legend.title = element_text(size = 10),
    #legend.text = element_text(size = 9)
  )



############################
# 12. OPTIONAL: VISUALIZE BOOTSTRAP DISTRIBUTIONS
############################
boot_long <- map_df(names(boot_results), function(link) {
  tibble(
    Model = link,
    L50 = boot_results[[link]]
  )
})

bootmean <- mean(boot_long$L50)
bootmed <- median(boot_long$L50)



ggplot(boot_long, aes(L50, fill = Model)) +
  geom_density(alpha = 0.4) +
  scale_fill_viridis(discrete=T, option="D") +
  geom_vline(
    xintercept = bootmed, color = "red",
    linetype = "dashed",
    linewidth=1
  ) +
  geom_vline(
    xintercept = bootmean, color = "black",
    linetype = "dashed",
    linewidth=1
  ) +
  #scale_fill_brewer(palette="Dark2") +
  theme_bw() +
  labs(
    title = "Bootstrap Distributions of L50",
    x = "L50 (Length at 50% Maturity)"
  )

ggsave("ScampLengthBasedTransitionL50Bootstrap062526.png", width = 10, height = 8, dpi = 1000)




############################
# 13. BOOTSTRAP PREDICTION UNCERTAINTY (RIBBONS)
############################


###### NOTE THAT THIS PLOTS MEAN VALUES. NOT MEDIAN LIKE THE LENGTH BASED MATURITY MODEL!!!



# Function to get predictions from one bootstrap draw
bootstrap_preds <- function(data, link, newdata) {
  
  # Resample
  boot_dat <- data %>%
    slice_sample(n = nrow(data), replace = TRUE)
  
  # Fit model
  model <- glm(Male ~ Length,
               family = binomial(link = link),
               data = boot_dat)
  
  # Predict
  predict(model, newdata = newdata, type = "response")
}

# Reduce number slightly for speed (you can increase if needed)
n_boot_pred <- 5000

# Generate bootstrap predictions for each model
boot_pred_list <- lapply(names(model_set), function(link) {
  
  cat("Bootstrapping predictions for:", link, "\n")
  
  preds_mat <- replicate(
    n_boot_pred,
    bootstrap_preds(dat, link, newdat)
  )
  
  # Summarize across bootstrap draws
  tibble(
    Length = newdat$Length,
    Model = link,
    fit = apply(preds_mat, 1, median, na.rm = TRUE),
    lower = apply(preds_mat, 1, quantile, 0.025, na.rm = TRUE),
    upper = apply(preds_mat, 1, quantile, 0.975, na.rm = TRUE)
  )
})

boot_preds <- bind_rows(boot_pred_list)








lengthtrans <- ggplot(dat, aes(Length, Male)) +
  geom_jitter(height = 0.03, alpha = 0.5, size=3) +
  
  geom_ribbon(
    data = boot_preds,
    aes(x = Length, ymin = lower, ymax = upper, fill = Model),
    alpha = 0.2,
    inherit.aes = FALSE
  ) +
  
  geom_line(data = preds, aes(y = fit, color = Model), linewidth = 1.5) +

  
  geom_vline(
    data = boot_summary,
    aes(xintercept = L50_mean, color = Model),
    linetype = "dashed",
    linewidth=1
  ) +
  scale_color_viridis(discrete=T, option="D",
                      labels = tools::toTitleCase,
                      name = "Link") +
  scale_fill_viridis(discrete=T, option="D",
                     labels = tools::toTitleCase,
                     name = "Link") +
  #scale_color_brewer(palette="Dark2") +
  #scale_fill_brewer(palette="Dark2") +
  theme_bw() +
  labs(
    x = "Fork length (mm)",
    y = "Probability male",
  ) +
  coord_cartesian(xlim = c(0, 850)) +
  theme(
    text = element_text(size = 14),
    legend.position = c(0.2, 0.8),       # Position inside plot
    legend.background = element_rect(    # Add background for readability
      fill = "white", 
      color = "black"
    ),
    #legend.title = element_text(size = 10),
    #legend.text = element_text(size = 9)
  )

lengthtrans



ggsave("ScampLengthBasedTransitionMedian062526.png", width = 10, height = 8, dpi = 1000)






#########
# merging maturity and transition plots togehter
########## REVISIT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
library(ggpubr)
library(grid)

p55 <- p5 + xlab(NULL)

final <- ggarrange(p55, lengthtrans, ncol = 1, nrow = 2,
                   label.x = c(0.125, 0.125),  # the x positions of the labels
                   label.y = c(0.98, 0.98),
                   font.label = list(size = 20),
                   labels = c("A", "B","C","D"))
final

#ggsave("SCPlengthbasedmerged.png", width = 10, height = 8, dpi = 1000)







































############# RESUME here!!! 
####need to change stuff below to the transition parameters and correct link1!!!



###############################################################################
####################### comparing to SEDAR68 DW
###############################################################
# ADD EXTERNAL MODEL: SEDAR68DW
###############################################################

# Old model parameters
###############################################################
# EXTERNAL MODEL: SEDAR68DW (PROBIT LINK)
###############################################################

a <- -7.9
b <- 2.17e-2

old_preds <- newdat %>%
  mutate(
    fit = pnorm(a + b * Length),   # probit link
    Model = "SEDAR68DW"
  )

# add to prediction table
preds <- bind_rows(preds, old_preds)

###############################################################
# L50 FOR PROBIT MODEL
###############################################################

# in probit, 50% occurs at linear predictor = 0
L50_old <- -a / b

boot_summary <- bind_rows(
  boot_summary,
  tibble(
    Model = "SEDAR68DW",
    L50_median = L50_old
  )
)
###############################################################
# OPTIONAL: improve legend labels for plotting
###############################################################

label_fun <- function(x) dplyr::recode(x,
                                       logit = "Logit",
                                       probit = "Probit",
                                       cloglog = "Cloglog",
                                       cauchit = "Cauchit",
                                       SEDAR68DW = "SEDAR68DW"
)

###############################################################
# USE IN PLOTS (replace your scale_* calls)
###############################################################

# scale_color_viridis(
#   discrete = TRUE,
#   option = "D",
#   labels = label_fun,
#   name = "Model"
# ) +
# scale_fill_viridis(
#   discrete = TRUE,
#   option = "D",
#   labels = label_fun,
#   name = "Model"
# )
###############################################################
############################
# 11. PLOT WITH BOOTSTRAP L50 + EXTERNAL MODEL
############################

ggplot(dat, aes(Length, Male)) +
  
  # raw data
  geom_jitter(height = 0.03, alpha = 0.5, size = 2) +
  
  # model curves (including SEDAR68DW via preds)
  geom_line(
    data = preds,
    aes(y = fit, color = Model),
    linewidth = 1.5
  ) +
  
  # L50 vertical lines (bootstrap + external model)
  geom_vline(
    data = boot_summary,
    aes(xintercept = L50_median, color = Model),
    linetype = "dashed",
    linewidth = 1
  ) +
  
  # color scale + labels
  scale_color_viridis(
    discrete = TRUE,
    option = "D",
    name = "Link",
    labels = function(x) dplyr::recode(x,
                                       logit = "Logit",
                                       probit = "Probit",
                                       cloglog = "Cloglog",
                                       cauchit = "Cauchit",
                                       SEDAR68DW = "SEDAR68DW"
    )
  ) +
  
  # theme + styling
  theme_bw() +
  labs(
    x = "Fork length (mm)",
    y = "Probability Male"
  ) +
  
  coord_cartesian(xlim = c(0, 850)) +
  
  theme(
    legend.position = c(0.8, 0.2),
    legend.background = element_rect(
      fill = "white",
      color = "black"
    )
  )

############################
# SAVE FIGURE
############################

ggsave(
  "ScampLengthBasedMaturity.png",
  width = 10,
  height = 8,
  dpi = 1000
)











###############################################################
# END OF SCRIPT
###############################################################

