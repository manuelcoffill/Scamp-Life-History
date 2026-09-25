# Load packages
library(ggplot2)
library(viridis)
library(readxl)
library(here)
library(dplyr)
library(tidyverse)
library(MuMIn)

# Read data
repro <- read_excel(here("Scamp staging.xlsx"))
age <- read_excel(here("Scamp ages for growth analysis w metadata.xlsx"))

# Check that IDs are stored as character strings
repro <- repro %>%
  mutate(`Cassette ID...9` = as.character(`Cassette ID...9`))

age <- age %>%
  mutate(`Fish ID` = as.character(`Fish ID`))

# Keep only the columns needed from the age dataset
age_lookup <- age %>%
  select(`Fish ID`, FINALAGE)

# Join age data onto reproductive data
repro_age <- repro %>%
  left_join(
    age_lookup,
    by = c("Cassette ID...9" = "Fish ID")
  )


# Rename the imported age column if desired
repro_age <- repro_age %>%
  rename(Age = FINALAGE)

# Check results
head(repro_age)


# How many fish matched?
sum(!is.na(repro_age$Age))
# only 74 :(



# Which reproductive IDs did not find a match?
unmatched <- repro_age %>%
  filter(is.na(Age)) %>%
  select(`Cassette ID...9`)

unmatched
print(unmatched, n=100)




# Optional: save the merged dataset
library(writexl)

# Save merged dataset as an Excel file
# write_xlsx(
#   repro_age,
#   path = "Scamp_reproductive_data_with_ages_062526.xlsx"
# )

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## some summary stats
dat <- repro_age %>%
  filter(!is.na(Age))

nrow(dat)

# sex breakdown
table(dat$HistoSex)

# are there any immature females?
fem <- dat[dat$HistoSex == "F",]
table(fem$ReproState)


# female summary
nrow(fem)
hist(fem$FL)
summary(fem$FL)

hist(fem$Age)
summary(fem$Age)


# male summary
mal <- dat[dat$HistoSex == "M",]
nrow(mal)

hist(mal$FL)
summary(mal$FL)

hist(mal$Age)
summary(mal$Age)



################################################################################
### comparing male to female
library(ggplot2)
library(dplyr)
library(viridis)

# Ensure HistoSex is a factor
dat <- dat %>%
  mutate(HistoSex = as.factor(HistoSex))

# Calculate sample sizes and positions for labels
n_labels <- dat %>%
  group_by(HistoSex) %>%
  summarise(
    n = n(),
    y_pos = 25,
    .groups = "drop"
  )

# Plot
ggplot(dat, aes(x = HistoSex, y = Age, fill = HistoSex)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 3) +
  geom_jitter(width = 0.15, alpha = 0.5, size=3) +
  geom_text(
    data = n_labels,
    aes(x = HistoSex, y = y_pos, label = paste0("n = ", n)),
    inherit.aes = FALSE
  ) +
  scale_fill_viridis_d(end=0.8) +
  coord_cartesian(ylim = c(0, 25)) +
  labs(
    x = "Sex",
    y = "Age (y)",
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 12),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )

ggsave("ScampSexSpecificAges062526.png", width = 10, height = 8, dpi = 1000)



ggplot(dat, aes(x = HistoSex, y = Age, fill = HistoSex)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.08, fill = "white", outlier.size = 3,outlier.alpha = 0.4) +
  geom_jitter(width = 0.1, alpha = 0.4, size=3) +
  geom_text(
    data = n_labels,
    aes(x = HistoSex, y = y_pos, label = paste0("n = ", n)),
    inherit.aes = FALSE
  ) +
  scale_fill_viridis_d(end=0.8) +
  labs(
    x = "Sex",
    y = "Age (y)"
  ) +
  theme_bw() +
  theme(
    legend.position = "none"
  )


#ggsave("ScampSexSpecificAges.png", width = 10, height = 8, dpi = 1000)




ggplot(dat, aes(x = Age, fill = HistoSex)) +
  geom_histogram(binwidth = 1, position = "identity", col="black") +
  scale_fill_viridis_d(end=0.8) +
  coord_cartesian(xlim = c(0, 25)) +
  labs(
    y = "Count",
    x = "Age (y)",
    fill="Sex"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 12),
    legend.position = c(0.1,0.85),
    legend.background = element_rect(fill = "white", color = "black")
  )

ggsave("ScampSexSpecificAgesHist062526.png", width = 10, height = 8, dpi = 1000)




###############################################################################
shapiro.test(dat$Age[dat$HistoSex == "F"])
shapiro.test(dat$Age[dat$HistoSex == "M"])
# neither is normally distributed
# The normality assumption for a t-test is technically violated.

t.test(Age ~ HistoSex, data = dat)
# Interpretation
#There is a highly significant difference in mean age between sexes
#Males are older on average than females
#Estimated difference ??? 6.3 years (males older)

#Important point
#Even though data are not normal, the t-test is still valid here because:

#  sample sizes are moderate (43 and 25)
#Welch correction is used
#t-test is robust to moderate non-normality



wilcox.test(Age ~ HistoSex, data = dat)

#Interpretation
#Confirms a significant difference in distributions
#Same conclusion as t-test
#Slightly smaller p-value difference is normal
#Warning message

#"cannot compute exact p-value with ties"

#This is expected because:
  
#  your ages are integers
#many tied values exist

#?????? This is NOT a problem and does not invalidate the test.











#4. Overall conclusion (what all of this means together)

#All tests agree:
  
#  ??? Females are significantly younger than males
#Mean females: ~10.4 years
#Mean males: ~16.7 years
#Difference: ~6 years
#Strong statistical evidence:
 # t-test: p ??? 5e-07
#Wilcoxon: p ??? 2e-06




