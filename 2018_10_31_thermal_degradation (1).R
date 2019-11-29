# Script to do initial processing of fluorogenic substrate proxy thermogenic degradation experiments

library(tidyverse) # Does all kinds of stuff with respect to data processing and visualization
library(lubridate) # handles dates and times

# Read in the data

d <- read_csv("data/2018_10_31_RFU - Sheet1.csv")

d <- d %>%
  mutate(time.posix = hms(time)) %>%
  group_by(sample, substrate, replicate, std, temperature) %>%
  mutate(elapsed.time = as.numeric(time - min(time))/ 3600) # hist(d$elapsed.time, breaks = 10)
  
# Visualize the raw data
p_raw <- ggplot(d, aes(x=elapsed.time, y=RFU, colour = substrate, shape = as.factor(replicate))) + 
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  expand_limits(ymin = 0) + 
  facet_grid(substrate~temperature, scales = "free_y")
print(p_raw)

# Calculate rates

mutate(d, time.posix = hms(time))
