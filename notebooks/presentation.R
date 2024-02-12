library("tidyverse")
library("xkcd")


## Bias illustration

expand_grid(
  N = 100L,
  Prev = c(0.1, 0.3, 0.5),
  RefSe = seq(1,0.75,by=-0.01),
  RefSp = seq(0.99,0.95,by=-0.01),
  TrueSe = 0.85,
  TrueSp = 0.99,
) |>
  mutate(TP = TrueSe*RefSe*Prev + (1-TrueSp)*(1-RefSp)*(1-Prev)) |>
  mutate(FP = TrueSe*(1-RefSe)*Prev + (1-TrueSp)*(RefSp)*(1-Prev)) |>
  mutate(FN = (1-TrueSe)*(RefSe)*Prev + TrueSp*(1-RefSp)*(1-Prev)) |>
  mutate(TN = (1-TrueSe)*(1-RefSe)*Prev + (TrueSp)*(RefSp)*(1-Prev)) |>
  mutate(Total = TP+FP+FN+TN, NewSe = TP/(TP+FN), NewSp = TN/(TN+FP)) |>
  pivot_longer(NewSe:NewSp) ->
  data


theme_set(theme_light())

data |>
  mutate(`True Prevalence` = factor(str_c(Prev*100,"%")), RefSp = str_c("Reference Sp=", RefSp*100, "%")) |>
  mutate(name = case_match(name, "NewSe" ~ "New Test Se", "NewSp" ~ "New Test Sp")) |>
  ggplot(aes(x=RefSe*100, y=value*100, col=`True Prevalence`)) +
  geom_hline(aes(x=NULL, y=NULL, col=NULL, yintercept=yint), tibble(name=c("New Test Se", "New Test Sp"), yint=c(0.85,0.99)*100), lty="dashed") +
  geom_line() +
  facet_grid(name ~ RefSp, scales="free_y") +
  theme(legend.position = "bottom") +
  xlab("Reference Se (%)") + ylab("Estimate for New Test")


data |>
  mutate(`True Prevalence` = factor(str_c(Prev*100,"%")), RefSp = str_c("Reference Sp=", RefSp*100, "%")) |>
  mutate(name = case_match(name, "NewSe" ~ "New Test Se", "NewSp" ~ "New Test Sp")) |>
  mutate(bias = value - if_else(name=="New Test Se", TrueSe, TrueSp)) |>
  ggplot(aes(x=RefSe*100, y=bias*100, col=`True Prevalence`)) +
  geom_hline(yintercept=0, lty="dashed") +
  geom_line() +
  facet_grid(name ~ RefSp, scales="free_y") +
  theme(legend.position = "bottom") +
  xlab("Reference Se (%)") + ylab("Bias for New Test (% points)")
ggsave("notebooks/bias.pdf", width=8, height=4)


## ROC

tibble(
  mean_n = 0,
  mean_p = 3,
  sd_n = 1,
  sd_p = 1
) |>
  expand_grid(cutoff = seq(-10,10,by=0.01)) |>
  mutate(Se = pnorm(cutoff, mean_p, sd_p, lower.tail=FALSE)) |>
  mutate(Sp = pnorm(cutoff, mean_n, sd_n, lower.tail=TRUE)) |>
  ggplot(aes((1-Sp)*100, Se*100)) +
  geom_abline(slope=1, intercept=0, lty="dashed") +
  geom_line() +
  ylab("Sensitivity (%)") + xlab("1 - Specificity (%)")
ggsave("notebooks/ROC.pdf", width=5, height=4)

  Youden = 1.6-1, Se = seq(0.6, 1, by=0.01)) |>
  mutate(Sp = Youden+1-Se) |>
  ggplot(aes(x=1-Sp, y=Se)) +
  geom_line() +
  xlim(0,1) + ylim(0,1)


## Fun illustration
datascaled <- data.frame(x=c(-3,3),y=c(-30,30))
p <- ggplot(data=datascaled, aes(x=x,y=y)) + geom_point()
xrange <- range(datascaled$x)
yrange <- range(datascaled$y)
ratioxy <- diff(xrange) / diff(yrange)

mapping <- aes(x=x,
  y=y,
  scale=scale,
  ratioxy=ratioxy,
  angleofspine = angleofspine,
  anglerighthumerus = anglerighthumerus,
  anglelefthumerus = anglelefthumerus,
  anglerightradius = anglerightradius,
  angleleftradius = angleleftradius,
  anglerightleg =  anglerightleg,
  angleleftleg = angleleftleg,
  angleofneck = angleofneck,
  color = color )

dataman <- data.frame( x= c(-1,0,1), y=c(-10,0,10),
  scale = c(10,7,5),
  ratioxy = ratioxy,
  angleofspine =  seq(- pi / 2, -pi/2 + pi/8, l=3) ,
  anglerighthumerus = -pi/6,
  anglelefthumerus = pi + pi/6,
  anglerightradius = 0,
  angleleftradius = runif(3,- pi/4, pi/4),
  angleleftleg = 3*pi/2  + pi / 12 ,
  anglerightleg = 3*pi/2  - pi / 12,
  angleofneck = runif(3, min = 3 * pi / 2 - pi/10 , max = 3 * pi / 2 + pi/10),
  color=c("A","B","C"))

p + xkcdman(mapping,dataman)


## Potential bias

expand_grid(
  Prevalence = c(0.1, 0.3, 0.5),
  GoldSe = seq()
)
