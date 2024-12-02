#setwd("U:/My Documents/AU_NCRR/1_iPSYCH design/Illustrations")

# ==============================================================================
# Title: Lexisplot illustrating iPSYCH cohorts and follow-up (sFigure 1)
# Author: TW
# Date: Nov 2024
# Description: Plot lexis grid and layers using lexis_grid() and lexis_polygon()
#===============================================================================

# Loaded packages --------------------------------------------------------------
library(LexisPlotR)
library(ggplot2)

# empty lexis diagram with grid
lexis <- lexis_grid(
  1980,
  2025,
  0,
  42,
  delta = 5,
  lwd = 0.3,
  force_equal = TRUE) + 
  theme(panel.grid.major = element_line(color = "gray", linetype = "dashed"))
lexis
          
# Study cohorts:
# born 1981-2008, followed from age 1+, followed to end of 2021
# born 1981-2008, followed from age 1+, followed to end of 2015
# born 1981-2005, followed from age 1+, followed to end of 2012
          
d1994 <- as.Date("1994-01-01")
d2012 <- as.Date("2012-12-31")
d2015 <- as.Date("2015-12-31")
d2021 <- as.Date("2021-12-31")

# adding 3 polygons
lexis1 <- lexis_polygon(
  lg = lexis,
  x = c("1982-05-01", "2010-01-01", "2021-12-31", "2021-12-31"),
  y = c(1, 1, 13, 40.67),
  fill = "grey",
  alpha = 0.9
)

lexis2 <- lexis_polygon(
  lg = lexis1,
  x = c("1982-05-01", "2010-01-01", "2015-12-31", "2015-12-31"),
  y = c(1, 1, 7, 34.67),
  fill = "darkgrey",
  alpha = 0.9
)
lexis2

lexis3 <- lexis_polygon(
  lg = lexis2,
  x = c("1982-05-01", "2007-01-01", "2012-12-31", "2012-12-31"),
  y = c(1, 1, 7, 31.67),
  fill = "lightgrey",
  alpha = 0.6
) +
  geom_vline(xintercept = as.numeric(d1994), linetype = "dashed") +
  geom_vline(xintercept = as.numeric(d2012), linetype = "dashed") +
  geom_vline(xintercept = as.numeric(d2015)) +
  geom_vline(xintercept = as.numeric(d2021)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_hline(yintercept = 10, linetype = "dashed") +
  labs(y = "Age in years", x = "Calendar year") +
  ggtitle("iPSYCH cohorts and follow-up time")
lexis3

# adding text in plot area
lexis4 <- lexis3 +
  geom_text(
    data = data.frame(x = as.Date("2008-01-01"), y = 5, label = "iPSYCH2012"),
    aes(x = x, y = y, label = label),
    color = "white",
    size = 3.5,
    angle = 45,
    vjust = 0.5,
    hjust = 0.5
  ) +
  geom_text(
    data = data.frame(x = as.Date("2012-06-01"), y = 5, label = "iPSYCH2015"),
    aes(x = x, y = y, label = label),
    color = "white",
    size = 3.5,
    angle = 45,
    vjust = 0.5,
    hjust = 0.5
  ) +
  geom_text(
    data = data.frame(x = as.Date("2018-09-01"), y = 14, label = "   Extended\nfollow-up"),
    aes(x = x, y = y, label = label),
    color = "white",
    size = 3.5,
    angle = 45,
    vjust = 0.5,
    hjust = 0.5
  )
lexis4

ggsave("sfig1_lexis_iPSYCH.png", width = 7,height = 6,dpi = 1200)
          





