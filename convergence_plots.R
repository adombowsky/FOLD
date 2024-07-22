# packages
require(cowplot)
require(mclust)
require(RColorBrewer)
# sourcing function
source("extra_R/convergence_line.R")
l1 <- convergence_line(n=100,seed=1,S=20000)
l2 <- convergence_line(n=250,seed=1,S=20000)
l3 <- convergence_line(n=500,seed=1,S=20000)
l4 <- convergence_line(n=1000,seed=1,S=20000)
risk_legend <- get_legend(l1$plot)
plot_grid(l1$plot + theme(legend.position = "none"),
          l2$plot + theme(legend.position = "none"),
          l3$plot + theme(legend.position = "none"),
          l4$plot + theme(legend.position = "none"),
          risk_legend,
          nrow=1,
          rel_widths = c(1,1,1,1,0.5),
          labels=c("(a)", "(b)", "(c)", "(d)"))
# comparing oracles
adjustedRandIndex(l1$c.fold, l1$c.oracle)
adjustedRandIndex(l2$c.fold, l2$c.oracle)
adjustedRandIndex(l3$c.fold, l3$c.oracle)
adjustedRandIndex(l4$c.fold, l4$c.oracle)
table(l4$c.fold, l4$c.oracle)
# computing number of clusters
table(l1$c.fold, l1$c.oracle)
table(l2$c.fold, l2$c.oracle)
table(l3$c.fold, l3$c.oracle)
table(l4$c.fold, l4$c.oracle)
# plotting for n=1000
plot.df <- data.frame(x1 = l4$data[,1],
                      x2 = l4$data[,2],
                      FOLD = factor(l4$c.fold),
                      oracle=factor(l4$c.oracle),
                      binder = factor(l4$c.B),
                      binderoracle = factor(l4$c.B.oracle))
fold_plot <- ggplot(plot.df,aes(x=x1,y=x2,col=FOLD)) + geom_point(size = 1.5) +
  theme_bw() +
  scale_color_brewer(palette = "Dark2") +
  xlab(" ") + ylab(" ") +
  labs(title="FOLD") +
  theme(panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        text = element_text(size=15))
oracle_plot <- ggplot(plot.df,aes(x=x1,y=x2,col=oracle)) + geom_point(size = 1.5) +
  theme_bw() +
  xlab(" ") + ylab(" ") +
  scale_color_brewer(palette = "Dark2") +
  labs(title="Oracle FOLD") +
  theme(panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        text = element_text(size=15))
binder_plot <- ggplot(plot.df,aes(x=x1,y=x2,col=binder)) + geom_point(size = 1.5) +
  theme_bw() +
  xlab(" ") + ylab(" ") +
  scale_color_brewer(palette = "Dark2") +
  labs(title="Binder") +
  theme(panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        text = element_text(size=15))
binderoracle_plot <- ggplot(plot.df,aes(x=x1,y=x2,col=binderoracle)) + geom_point(size = 1.5) +
  theme_bw() +
  xlab(" ") + ylab(" ") +
  scale_color_brewer(palette = "Dark2") +
  labs(title="Oracle Binder") +
  theme(panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        text = element_text(size=15))
plot_grid(fold_plot,oracle_plot,binder_plot,binderoracle_plot,
          labels = c("(a)", "(b)", "(c)", "(d)"))
# number of clusters
table(l4$c.fold)
table(l4$c.oracle)
table(l4$c.B)
table(l4$c.B.oracle)
