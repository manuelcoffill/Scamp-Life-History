###############################################################################
# Compare SEDAR68 and Bayesian Growth Models
###############################################################################

library(ggplot2)
library(readxl)
library(dplyr)

###############################################################################
# 1. SEDAR68DW VBGF (convert cm -> mm)
###############################################################################

Linf <- 70.22
k <- 0.134
t0 <- -1.762

age <- seq(0, 34, by = 0.1)

length_vbgf <- Linf * (1 - exp(-k * (age - t0)))

df_vbgf <- data.frame(
  age = age,
  length = length_vbgf * 10,   # cm -> mm
  lower = NA,
  upper = NA,
  method = "SEDAR68DW"
)


###############################################################################
# 2. SEDAR68OA Assessment Growth (convert cm -> mm)
###############################################################################

L_Amin <- 24.6953
L_Amax <- 77.2886
K_assess <- 0.0727
cv_assess <- 0.1298

length_assess <- ifelse(
  age < 1,
  L_Amin * age,
  L_Amax - (L_Amax - L_Amin) * exp(-K_assess * (age - 1))
)

sigma_assess <- sqrt(log(1 + cv_assess^2))

upper_assess <- ifelse(
  age < 1,
  NA,
  length_assess * exp(1.96 * sigma_assess)
)

lower_assess <- ifelse(
  age < 1,
  NA,
  length_assess * exp(-1.96 * sigma_assess)
)


df_assess <- data.frame(
  age = age,
  length = length_assess * 10,  # cm -> mm
  lower = lower_assess * 10,    # cm -> mm
  upper = upper_assess * 10,    # cm -> mm
  method = "SEDAR68OA"
)


###############################################################################
# 3. Bayesian Growth Models
###############################################################################

path <- "C:\\Users\\manue\\OneDrive\\Desktop\\Siders VBGF - v2"
setwd(path)

dat <- read_excel("LengthAtAgePreds90CRI.xlsx")


# Known t0 model
df_known <- data.frame(
  age = dat$Age,
  length = dat$Known.med,
  lower = dat$Known.low.CRI,
  upper = dat$Known.high.CRI,
  method = "Bayesian Known L0"
)


# Free t0 model
df_free <- data.frame(
  age = dat$Age,
  length = dat$Free.med.CRI,
  lower = dat$Free.low.CRI,
  upper = dat$Free.high.CRI,
  method = "Bayesian Free L0"
)


###############################################################################
# 4. Combine all models
###############################################################################

growth_all <- bind_rows(
  df_vbgf,
  df_assess,
  df_known,
  df_free
)


###############################################################################
# 5. Plot
###############################################################################

ggplot(growth_all,
       aes(x = age,
           y = length,
           color = method,
           fill = method)) +
  
  # uncertainty ribbons
  geom_ribbon(aes(ymin = lower,
                  ymax = upper),
              alpha = 0.20,
              na.rm = TRUE,
              color = NA) +
  
  # growth curves
  geom_line(linewidth = 1.5) +
  
  labs(
    x = "Age (y)",
    y = "Fork length (mm)",
    color = "Growth model",
    fill = "Growth model"
  ) +
  
  scale_color_manual(
    values = c(
      "SEDAR68DW" = "#D55E00",       # vermillion
      "SEDAR68OA" = "#0072B2",       # blue
      "Bayesian Known L0" = "#009E73", # green
      "Bayesian Free L0" = "#CC79A7"  # purple-pink
    ),
    labels = c(
      "SEDAR68DW" = "SEDAR68DW",
      "SEDAR68OA" = "SEDAR68OA",
      "Bayesian Known L0" = expression("Bayesian Known " * italic(L)[0]),
      "Bayesian Free L0" = expression("Bayesian Free " * italic(L)[0])
    )
  ) +
  
  scale_fill_manual(
    values = c(
      "SEDAR68DW" = "#D55E00",       # vermillion
      "SEDAR68OA" = "#0072B2",       # blue
      "Bayesian Known L0" = "#009E73", # green
      "Bayesian Free L0" = "#CC79A7"  # purple-pink
    ),
    labels = c(
      "SEDAR68DW" = "SEDAR68DW",
      "SEDAR68OA" = "SEDAR68OA",
      "Bayesian Known L0" = expression("Bayesian Known " * italic(L)[0]),
      "Bayesian Free L0" = expression("Bayesian Free " * italic(L)[0])
    )
  ) +
  
  scale_x_continuous(
    breaks = seq(0, 40, by = 5)
  ) +
  
  scale_y_continuous(
    breaks = seq(0, 1200, by = 200)
  ) +
  
  theme_bw() +
  
  theme(
    text = element_text(size = 14),
    legend.position = c(0.78, 0.20),
    legend.background = element_rect(
      fill = "white",
      color = "black"
    )
  )

# ggsave("ScampGrowthModelsComparison90CRI.png", width = 10, height = 8, dpi = 1000)


