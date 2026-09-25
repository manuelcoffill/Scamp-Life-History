library(ggplot2)
library(viridis)
library(readxl)
library(here)
library(dplyr)

dat <- read_excel(here("Scamp staging.xlsx"))

################################################################################
################################################################################

## quick usmmary of lengths by sex
trans <- dat[dat$HistoSex == "T",]
males <- dat[dat$HistoSex == "M",]

summary(trans$FL)
summary(males$FL)



############################## monthly phase plot
# female df
females <- dat[dat$HistoSex == "F",]
nrow(females)
summary(females$FL)

table(females$ReproPhase)

# Proportions of phases
prop.table(table(females$ReproPhase))

# remove immature fish
table(females$ReproState)
matfem <- females[!females$ReproState == "Immature",]
table(matfem$ReproState) #should be all mature
nrow(matfem)


# remove fish missing phase assignments
is.na(matfem$ReproPhase)
table(matfem$ReproPhase) # remove unassigned, if any
matfem <- matfem[!is.na(matfem$ReproPhase),]
nrow(matfem)

# monthly data for plot
q3dat <- count(matfem, ReproPhase, Month)
q3dat

# monthly sample sizes
q3monthn <- count(matfem, Month)
q3monthn

# factor and sort so it makes sense biologically
table(q3dat$ReproPhase)
q3dat$ReproPhase <- factor(q3dat$ReproPhase,
                            levels = c("Regenerating", "Early developing", 
                                       "Late developing","Spawning A",
                                       "Spawning S","Regressing"))

names(q3dat)[names(q3dat) == "ReproPhase"] <- "Phase"

# proportion
ggplot(q3dat, aes(fill=Phase, y=n, x=Month)) + 
  geom_bar(position="fill", stat="identity", col="black") +
  scale_fill_viridis(discrete = T, option = "D") +
  #scale_fill_brewer(palette = "Dark2") +
  xlab("Month") + ylab("Proportion") +
  scale_x_continuous(breaks = seq(1, 12, by = 1)) +
  geom_text(aes(Month, 1.05, label = n, fill = NULL), data = q3monthn) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    text = element_text(size = 12)
  ) +
  guides(
    fill = guide_legend(nrow = 1)
  )

#ggsave("SCPFemaleMonthlyPhases062526n75.png", width = 10, height = 8, dpi = 1000)


#### consider binning all spawning subphases


################################################################################
################################################################################
################################################################################
####### sex-specific length frequency
table(dat$HistoSex)

# filter unsexed fish
sexed <- dat[!dat$HistoSex == "??",]
nrow(sexed)
table(sexed$HistoSex)

# current male proportion
pct_male <- nrow(sexed[sexed$HistoSex == "M",]) / nrow(sexed) * 100
label_text <- paste0(round(pct_male), "% male")
label_text


# hist
library(grid)

ggplot(sexed, aes(x=FL, fill=HistoSex)) +
  geom_histogram(col="black", binwidth=50, boundary=250,
                 position = "stack") +
  scale_fill_viridis(discrete=T,option = "E") +
  #scale_fill_brewer(palette="Dark2") +
  #scale_fill_manual(values = c("#0072B2", "#E69F00", "#009E73")) +
  scale_x_continuous(breaks = seq(300,800,100)) +
  coord_cartesian(xlim = c(250, 850)) +
  theme_bw() +
  labs(
    x = "Fork length (mm)",
    y = "Count"
  ) +
  theme(
    text = element_text(size = 12),
    legend.position = c(0.8, 0.8),       # Position inside plot
    legend.background = element_rect(    # Add background for readability
      fill = "white", 
      color = "black"
    ),
    #legend.title = element_text(size = 10),
    #legend.text = element_text(size = 9)
  ) +
  labs(fill = "Sex",
       y="Frequency")

#ggsave("SCPLengthFreqStacked062526.png", width = 10, height = 8, dpi = 1000)



# consider including male proportion in plot
#annotate(
 #   "text",
  #  x = -Inf, y = Inf,
   # label = "33% male",
    #hjust = -0.2, vjust = 1.2,
    #size = 6
  #)




################################################################################
################################################################################
####### sex-specific length frequency (mirrored)

#### not worth it if there are transitionals

# define breaks to match your histogram exactly
breaks <- seq(250, 850, by = 50)

sexed2 <- sexed %>%
  mutate(
    # assign bins based on same boundary/binwidth logic
    bin_id = findInterval(FL, breaks, rightmost.closed = TRUE),
    bin_mid = breaks[bin_id] + 25
  ) %>%
  filter(bin_id > 0 & bin_id < length(breaks)) %>%
  count(bin_mid, HistoSex) %>%
  mutate(
    # males below (negative), females above (positive)
    count_plot = ifelse(HistoSex == "M", -n, n)
  )

# recompute male proportion for annotation
pct_male <- nrow(sexed[sexed$HistoSex == "M",]) / nrow(sexed) * 100
label_text <- paste0(round(pct_male), "% male")



ggplot(sexed2, aes(x = bin_mid, y = count_plot, fill = HistoSex)) +
  geom_col(color = "black", width = 50) +
  scale_y_continuous(labels = abs, breaks = c(-10,-5,0,5,10,15,20)) +
  scale_x_continuous(breaks = seq(300,800,100)) +
  geom_hline(yintercept = 0, linetype = "solid") +
  scale_fill_viridis(discrete = TRUE, option = "E") +
  coord_cartesian(xlim = c(250, 850), ylim = c(-10,20)) +
  theme_bw() +
  labs(
    x = "Fork length (mm)",
    y = "Frequency",
    fill = "Sex"
  ) +
  theme(
    legend.position = c(0.8, 0.8),
    legend.background = element_rect(fill = "white", color = "black")
  ) +
  annotate(
    "text",
    x = 300, y = -6,
    label = label_text,
    size = 5
  )

#ggsave("SCPLengthFreqMirrored060826.png", width = 10, height = 8, dpi = 1000)


