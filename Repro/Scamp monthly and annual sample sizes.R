library(ggplot2)
library(viridis)
library(readxl)
library(here)
library(dplyr)
library(tidyr)

dat <- read_excel(here("Scamp staging.xlsx"))
nrow(dat)



# yearly and monthly sample breakdown
table(dat$Month)
table(dat$Year)

data <- dat %>% group_by(Year, Month) %>% tally

data %>% expand(Year, Month, full_seq(n, 0))
data <- as.data.frame(data)
data$Year <- as.numeric(data$Year)
data$Month <- as.numeric(data$Month)


## df with all possible combinations of years and months
Year <- seq(2022,2025,1)
Month <- seq(01,12,01)
dd <- expand.grid(Year=Year, Month=Month)
dd <- as.data.frame(dd)
dd


# Fill in values based on multiple conditions of year and month
ddd <- dd %>%
  left_join(data, by = c("Year", "Month"))
ddd$n[is.na(ddd$n)] <- 0
ddd




################################################################################
# Give extreme colors:
ggplot(ddd, aes(Month, Year, fill= n)) + 
  geom_tile() +
  scale_fill_gradient(low="white", high="slateblue") +
  geom_text(aes(label = n), size = 4) +
  scale_x_continuous(breaks = seq(1,12,1)) +
  theme_bw() +
  guides(fill = "none") +
  theme(
    axis.title = element_text(size = 15),
    axis.text  = element_text(size = 15)
  )

# Color Brewer palette
ggplot(ddd, aes(Month, Year, fill= n)) + 
  geom_tile() +
  scale_fill_viridis(discrete=FALSE, option="G") +
  geom_text(aes(label = ifelse(n > 0, n, "")), 
            color = "white", 
            size = 3) +
  scale_x_continuous(breaks = seq(1,12,1)) +
  theme_bw() +
  guides(fill = "none") +
  theme(
    axis.title = element_text(size = 15),
    axis.text  = element_text(size = 15)
  )

# Color Brewer palette
ggplot(ddd, aes(Month, Year, fill= n)) + 
  geom_tile() +
  scale_fill_distiller(palette = "Purples", direction = 1) +
  geom_text(aes(label = n), size = 4) +
  scale_x_continuous(breaks = seq(1,12,1)) +
  theme_bw() +
  guides(fill = "none") +
  theme(
    axis.title = element_text(size = 15),
    axis.text  = element_text(size = 15)
  )
# play with palettes Blues and Purples





###############################
# adding yearly and monthly total

# ---- Ensure consistent types ----
ddd <- ddd %>%
  mutate(Year = as.character(Year))

# ---- Yearly totals (new column: Month = 13) ----
year_totals <- ddd %>%
  group_by(Year) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  mutate(Month = 13)

# ---- Monthly totals (new row: Year = "Total") ----
month_totals <- ddd %>%
  group_by(Month) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  mutate(Year = "Total")

# ---- Combine all ----
data_totals <- bind_rows(ddd, year_totals, month_totals)

# ---- Order Year so "Total" is last ----
data_totals$Year <- factor(
  data_totals$Year,
  levels = c(sort(unique(ddd$Year)), "Total")
)

# ---- Plot ----
ggplot(data_totals, aes(x = Month, y = Year, fill = n)) +
  geom_tile(color = "white") +
  geom_text(aes(label = n), size = 3) +
  scale_x_continuous(
    breaks = c(1:12, 13),
    labels = c(1:12, "Total")
  ) +
  labs(x = "Month", y = "Year", fill = "Sample size") +
  theme_minimal()


# Color Brewer palette
ggplot(data_totals, aes(Month, Year, fill= n)) + 
  geom_tile() +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  geom_text(aes(label = n), size = 4) +
  scale_x_continuous(
    breaks = c(1:12, 13),
    labels = c(1:12, "Total")
  ) +
  theme_bw() +
  guides(fill = "none") +
  theme(
    axis.title = element_text(size = 15),
    axis.text  = element_text(size = 15)
  )
# play with palettes Blues and Purples

ggsave("ScampStagedSampleSizesFinal.png", width = 10, height = 8, dpi = 1000)





#####################################
### consider adding the other things that i did on the age analysis
# e.g., sources
