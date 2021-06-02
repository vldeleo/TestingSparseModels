library(ggplot2)
library(ggridges)
library(reshape2)

# how well do different models do at predicting unobserved values
predR2many <- read.csv("many_PredR2.csv")
predR2few <- read.csv("few_PredR2.csv")
predR2many$sparse <- "dense"
predR2few$sparse <- "sparse"
predR2 <- rbind(predR2many, predR2few)
predR2$X <- as.character(predR2$X) # turns out I can't use an integer as an ID column
predR2 <- melt(predR2, id = c("X", "sparse"), variable.name = "model")
png(filename = "predR2comparison_6.2.21.png", width = 600, height = 400)
ggplot(predR2, aes(x = model, y = value, color = sparse, fill = sparse)) +  geom_boxplot(width=0.5,lwd=1) + ggtitle("R2 of Predicted vs Actual OOS") + theme_bw()
dev.off()

# how well do different models do at finding Betas
BetaR2many <- read.csv("many_BetaR2.csv")
BetaR2few <- read.csv("few_BetaR2.csv")
BetaR2many$sparse <- "dense"
BetaR2few$sparse <- "sparse"
BetaR2 <- rbind(BetaR2many, BetaR2few)
BetaR2$X <- as.character(BetaR2$X) # turns out I can't use an integer as an ID column
BetaR2 <- melt(BetaR2, id = c("X", "sparse"), variable.name = "model")
png(filename = "BetaR2comparison_6.2.21.png", width = 600, height = 400)
ggplot(BetaR2, aes(x = model, y = value, color = sparse, fill = sparse)) +  geom_boxplot(width=0.5,lwd=1) + ggtitle("R2 of Model vs Actual Beta") + theme_bw()
dev.off()



# how about distributions of R2?
# and better yet, how distributions of R2 change from few to many predictors for each model

ggplot(predR2, aes(x = Variable, y = model, color = sparse, point_color = sparse, fill = sparse)) +
  geom_density_ridges(
    jittered_points = TRUE, scale = .95, rel_min_height = .01,
    point_shape = "|", point_size = 3, size = 0.25,
    position = position_points_jitter(height = 0)
  ) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0), name = "r^2") +
  scale_fill_manual(values = c("darkblue", "darkgold"), labels = c("sparse", "dense")) +
  scale_color_manual(values = c("darkblue", "darkgold"), guide = "none") +
  #scale_discrete_manual("point_color", values = c("#D55E00", "#0072B2"), guide = "none") +
  coord_cartesian(clip = "off") +
  guides(fill = guide_legend(
    override.aes = list(
      fill = c("darkblue", "darkgold"),
      color = NA, point_color = NA)
  )
  ) +
  ggtitle("r^2 Predictions of OOS") +
  theme_ridges(center = TRUE)

png(filename = "R2dens_5.28.21.png", width = 600, height = 400)
density(predR2tog, predR2tog)
dev.off()
