# libraries
library(KODAMA)
library(ggplot2)
# sampling
set.seed(1996)
N <- 300
data <- KODAMA::spirals(n = c(N/3,N/3,N/3), sd = c(0.3,0.25,0.15))
truth <- rep(1:3, each = N/3)
y <- scale(data)
data_df <- data.frame(y1 = y[,1], y2 = y[,2], s0 = as.factor(truth))
# plotting
ggplot(data_df, aes(x=y1, y = y2, color = s0)) + geom_point() + theme_bw() +
  xlab(" ") + ylab(" ") +
  labs(title = " ", caption = "Example of the spirals data.") +
  theme(panel.border = element_blank(),
        #panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.caption = element_text(hjust=0.5),
        text = element_text(size=14))
