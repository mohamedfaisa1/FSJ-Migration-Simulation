library(reticulate)
library(ggplot2)
library(cowplot)
np <- import("numpy")
theme_set(theme_cowplot())
mytheme = list(
  theme_classic()+
    theme(panel.background = element_blank(),
          strip.background = element_rect(colour=NA, fill=NA),
          panel.border = element_rect(fill = NA, color = "black"),
          legend.title = element_blank(),legend.position="bottom", 
          strip.text = element_text(face="bold", size=9),
          axis.text=element_text(face="bold"),
          axis.title = element_text(face="bold"),
          plot.title = element_text(face = "bold", hjust = 0.5,size=13))
)

# Fst
AC_n <- np$load("AC_no_mig_0.05.npy")
AC_u <- np$load("AC_unid_0.05.npy")
AC_r <- np$load("AC_recip_0.05.npy")
AC_d <- np$load("AC_dirmig_0.05.npy")

OC_n <- np$load("OC_no_mig_0.05.npy")
OC_u <- np$load("OC_unid_0.05.npy")
OC_r <- np$load("OC_recip_0.05.npy")
OC_d <- np$load("OC_dirmig_0.05.npy")

OA_n <- np$load("OA_no_mig_0.05.npy")
OA_u <- np$load("OA_unid_0.05.npy")
OA_r <- np$load("OA_recip_0.05.npy")
OA_d <- np$load("OA_dirmig_0.05.npy")


Fst <- c(AC_n, AC_u, AC_r, AC_d, OC_n, OC_u, OC_r, OC_d, OA_n, OA_u, OA_r, OA_d)

pair_pops <- c(rep("ABS/CST", 400), rep("ONF/CST", 400), rep("ONF/ABS", 400))

nm.tag <- rep("No Migration", 100)
ud.tag <- rep("Unidirectional Migration", 100)
rc.tag <- rep("Reciprocal Migration", 100)
dm.tag <- rep("Directional Migration", 100)

mgs <- c(nm.tag, ud.tag, rc.tag, dm.tag)
mgs <- c(mgs, mgs, mgs)


fst.df <- cbind.data.frame(x=pair_pops, y=Fst, m=mgs)

conf95 <- function(x){
  return(((1.96/sqrt(100))*sd(x)))
}

get_se <- function(v, p, m){
  res <- tapply(v, p, 
                function(x) tapply(x, m, function(w) conf95(w)))
  return(res)
}


get_mu <- function(v, p, m){
  res <- tapply(v, p, 
                function(x) tapply(x, m, function(w) mean(w)))
  return(res)
}


get_pvals <- function(v, p, m){
  res <- tapply(v, p, 
                function(x) tapply(x, m, function(w) return(w)))
}

gstat_vectorize <- function(w){
  p1 <- c(rep(w[[3]][2], 100), rep(w[[3]][4], 100),
          rep(w[[3]][3], 100), rep(w[[3]][1], 100))
  p2 <- c(rep(w[[2]][2], 100), rep(w[[2]][4], 100),
          rep(w[[2]][3], 100), rep(w[[2]][1], 100))
  p3 <- c(rep(w[[1]][2], 100), rep(w[[1]][4], 100),
          rep(w[[1]][3], 100), rep(w[[1]][1], 100))
  p <- c(p3,p1,p2)
  return(p)
}

fst.df$fm <- gstat_vectorize(get_mu(Fst, pair_pops, mig_tags))
fst.df$fse <- gstat_vectorize(get_se(Fst, pair_pops, mig_tags))


fst.plot <- ggplot(fst.df, aes(x=factor(x), y=y, color=factor(m)))+
   xlab("Pairwise Populations") + 
  ylab("Fst") +
  ggtitle("Pairwise Fst (m = 0.05)") +
  labs(color="Migration Scheme") + mytheme+
  geom_point(position=position_dodge(width=0.5), 
             alpha=0.2)+ 
  stat_summary(position=position_dodge(width=0.5),fun=mean, geom="point", 
               shape=20, size=8, alpha=0.2) +
  geom_errorbar(aes(ymin=fm-fse,
                    ymax=fm+fse),
                position=position_dodge(.5), width=0.4)

fst.plot

ggsave("fst_box.png", fst.plot)
s