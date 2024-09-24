# metaMest

Random-Effect Meta-Analysis with Robust Between-Study Variance via M-estimator. The package implement the estimation procedure of random-effect meta-analysis described in Hanada and Sugimoto (2024)


# Sample code of Hanada and Sugimoto (2024)

```{r}

# install.packages("openxlsx")
# install.packages("metafor")
# install.packages("rma.exact")
# devtools::install_github("metaMest")
library(openxlsx)
library(metafor)
library(rma.exact)
library(metaMest)

set.seed(1234)

dat01 <- data.frame(
  study = c("Prahalad(2006)", "Hinks(2006)", "Lindner(2007)"),
  logOR = c(-0.13, -0.24, -0.20),
  CIL95 = c(-0.33, -0.41, -0.47),
  CIU95 = c(0.08, -0.06, 0.07),
  SD = c(0.10459184, 0.08928571, 0.13775510)
)


yk <- dat01$logOR
sig2k <- dat01$SD^2
K <- length(yk)
alpha <- 0.05

dl <- rma.uni(yi=yk, vi=sig2k, method="DL")
sj <- rma.uni(yi=yk, vi=sig2k, method="SJ")

# DL
dl.orci <- c(exp(dl$ci.lb), exp(dl$ci.ub))
# SJ
sj.orci <- c(exp(sj$ci.lb), exp(sj$ci.ub))


# HKSJ-DL
wk_tau2dl <- sapply(1:K, function(k){1/(sig2k[k]+dl$tau2)})
v2_hksj <- sum(wk_tau2dl*(yk-as.numeric(dl$b))^2) / ((K-1)*sum(wk_tau2dl))
cildl_hksj <- dl$b + qt(alpha/2, df=K-1) * sqrt(v2_hksj)
ciudl_hksj <- dl$b + qt(1-alpha/2, df=K-1) * sqrt(v2_hksj)
ci_hksj_dl <- c(cildl_hksj, ciudl_hksj)

# HKSJ-SJ
wk_tau2sj <- sapply(1:K, function(k){1/(sig2k[k]+sj$tau2)})
v2_hksj <- sum(wk_tau2dl*(yk-as.numeric(sj$b))^2) / ((K-1)*sum(wk_tau2sj))
cilsj_hksj <- sj$b + qt(alpha/2, df=K-1) * sqrt(v2_hksj)
ciusj_hksj <- sj$b + qt(1-alpha/2, df=K-1) * sqrt(v2_hksj)
ci_hksj_sj <- c(cilsj_hksj, ciusj_hksj)


# DL with M-estimator
y <- 1:1000/1000
pval_dl <- ptheta_with_M_estimator(qval=y, yi=yk, vi=sig2k, theta=0, tau2=dl$tau2, method="DL")
alpha <- 0.05
yU <- y[which.min(abs(pval_dl$qy-(1-alpha/2)))]
cim_dl <- as.numeric(dl$b) + c(-1,1)*yU

# SJ with M-estimator
pval_sj <- ptheta_with_M_estimator(qval=y, yi=yk, vi=sig2k, theta=0, tau2=sj$tau2, method="SJ")
yU <- y[which.min(abs(pval_sj$qy-(1-alpha/2)))]
cim_sj <- as.numeric(sj$b) + c(-1,1)*yU

# HKSJ-DL with M-estimator
y <- c(-100:100/1000+dl$ci.lb, -100:100/1000+dl$ci.ub)
pval_hksj_dl <- ptheta_with_M_estimator(qval=y, yi=yk, vi=sig2k, theta=0, tau2=dl$tau2, method="HKSJ-DL")
cim_hksj_dl <- c(y[which.min(abs(pval_hksj_dl$qy - (alpha/2)))], y[which.min(abs(pval_hksj_dl$qy - (1-alpha/2)))])

# HKSJ-SJ with M-estimator
y <- c(-100:100/1000+sj$ci.lb, -100:100/1000+sj$ci.ub)
pval_hksj_sj <- ptheta_with_M_estimator(qval=y, yi=yk, vi=sig2k, theta=0, tau2=sj$tau2, method="HKSJ-SJ")
cim_hksj_sj <- c(y[which.min(abs(pval_hksj_sj$qy - (alpha/2)))], y[which.min(abs(pval_hksj_sj$qy - (1-alpha/2)))])

mic <- rma.exact(yi=yk, vi=sig2k)
mic

beta <- 0.0005
tau2bound <- as.numeric( c(0, quantile(rand_tau2dl_tau2(n=10000, K=K, tau2=pval_dl$tau2m, sig2k=sig2k), probs=c(1-beta)) ) )
mmic <- rma.exact(yi=yk, vi=sig2k, tau2.bounds=tau2bound)
mmic


result <- data.frame(Method=c("DL", "SJ", "HKSJ-DL", "HKSJ-SJ", "mDL", "mSJ", "HKSJ-mDL", "HKSJ-mSJ", "Mi", "mMi"),
                     round(
                       rbind(c(dl$b, dl$ci.lb, dl$ci.ub),
                             c(sj$b, sj$ci.lb, sj$ci.ub),
                             c(dl$b, cim_dl),
                             c(sj$b, cim_sj),
                             c(dl$b, ci_hksj_dl),
                             c(sj$b, ci_hksj_sj),
                             c(dl$b, cim_hksj_dl),
                             c(sj$b, cim_hksj_sj),
                             c((mic[2]+mic[1])/2, mic[1], mic[2]),
                             c((mmic[2]+mmic[1])/2, mmic[1], mmic[2])
                       ),3
                     )
                     )
colnames(result) <- c("Method", "Estimate", "CIL", "CIU")
result


### forest plot
# install.packages("forestplot")
# install.packages("dplyr")
library(forestplot)
library(dplyr)

pdat <- data.frame(Method = c("Prahalad (2006)", "Hinks (2006)", "Lindner (2007)"),
                   Estimate = yk,
                   CIL = round(yk - qnorm(1-alpha/2) * sqrt(sig2k), 3),
                   CIU = round(yk + qnorm(1-alpha/2) * sqrt(sig2k), 3)
)

dummy <- data.frame(Method = "Methods", Estimate = 0, CIL = 0, CIU = 0)

pdat <- rbind(pdat, dummy, result)
pdat$Method <- paste("   ", pdat$Method, "   ", sep="")
pdat$clab <- paste(pdat$Estimate, " (", pdat$CIL, ",", pdat$CIU, ")", sep="")
pdat[pdat$Method=="   Methods   ",] <- data.frame("Methods", NA, NA, NA, NA)

pdat %>% 
  forestplot(
    mean = Estimate,
    lower = CIL,
    upper = CIU,
    labeltext = c(Method, clab),
    lwd.ci = 3,
    boxsize = 0.25,
    ci.vertices.height = 0.2,
    vertices = TRUE,
    xticks = -5:1/10,
    xticks.digits = 1,
    xlim = c(-0.5, 0.1)
  ) %>%
  fp_set_style(
    box = "black",
    line = "black",
    summary = "black",
    align = "lccccc",
    hrz_lines = "black",
    txt_gp = fpTxtGp(ticks = gpar(cex = 1))
  ) %>%
  fp_add_header(
    Method = " Study",
    clab = "log OR (95%CI)"
  ) %>%
  fp_decorate_graph(graph.pos = 2) %>%
  fp_set_zebra_style("#EFEFEF") -> p1
p1
```


Hanada, K., & Sugimoto, T. (2024). Random-Effect Meta-Analysis with Robust Between-Study Variance (Submitted).


