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

het_dm_abs <- np$load("hetero_cts_abs_dm.npy")
het_dm_cst <- np$load("hetero_cts_cst_dm.npy")
het_dm_onf <- np$load("hetero_cts_onf_dm.npy")

het_nm_abs <- np$load("hetero_cts_abs_nm.npy")
het_nm_cst <- np$load("hetero_cts_cst_nm.npy")
het_nm_onf <- np$load("hetero_cts_onf_nm.npy")

het_ud_abs <- np$load("hetero_cts_abs_ud.npy")
het_ud_cst <- np$load("hetero_cts_cst_ud.npy")
het_ud_onf <- np$load("hetero_cts_onf_ud.npy")

het_rc_abs <- np$load("hetero_cts_abs_rc.npy")
het_rc_cst <- np$load("hetero_cts_cst_rc.npy")
het_rc_onf <- np$load("hetero_cts_onf_rc.npy")

hom_dm_abs <- np$load("homo_cts_abs_dm.npy")
hom_dm_cst <- np$load("homo_cts_cst_dm.npy")
hom_dm_onf <- np$load("homo_cts_onf_dm.npy")

hom_nm_abs <- np$load("homo_cts_abs_nm.npy")
hom_nm_cst <- np$load("homo_cts_cst_nm.npy")
hom_nm_onf <- np$load("homo_cts_onf_nm.npy")

hom_ud_abs <- np$load("homo_cts_abs_ud.npy")
hom_ud_cst <- np$load("homo_cts_cst_ud.npy")
hom_ud_onf <- np$load("homo_cts_onf_ud.npy")

hom_rc_abs <- np$load("homo_cts_abs_rc.npy")
hom_rc_cst <- np$load("homo_cts_cst_rc.npy")
hom_rc_onf <- np$load("homo_cts_onf_rc.npy")

nm.tag.cst <- rep("No Migration", 50)
nm.tag.abs <- rep("No Migration", 50)
nm.tag.onf <- rep("No Migration", 50)

dm.tag.cst <- rep("Directional Migration", 50)
dm.tag.abs <- rep("Directional Migration", 50)
dm.tag.onf <- rep("Directional Migration", 50)

ud.tag.cst <- rep("Unidirectional Migration", 50)
ud.tag.abs <- rep("Unidirectional Migration", 50)
ud.tag.onf <- rep("Unidirectional Migration", 50)

rc.tag.cst <- rep("Reciprocal Migration", 50)
rc.tag.abs <- rep("Reciprocal Migration", 50)
rc.tag.onf <- rep("Reciprocal Migration", 50)

cst.tag <- rep("CST", 50)
abs.tag <- rep("ABS", 50)
onf.tag <- rep("ONF", 50)

conf950 <- function(x){
  return(((1.96/sqrt(50))*sd(x)))
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
  p <- c(rep(w[[1]][1], 50), rep(w[[2]][1], 50),
         rep(w[[3]][1], 50), rep(w[[1]][2], 50),
         rep(w[[2]][2], 50), rep(w[[3]][2], 50),
         rep(w[[1]][4], 50), rep(w[[2]][4], 50),
         rep(w[[3]][4], 50), rep(w[[1]][3], 50),
         rep(w[[2]][3], 50), rep(w[[3]][3], 50))
  
  return(p)
}

het_vals <- c(het_dm_abs, het_dm_cst, het_dm_onf, het_nm_abs, het_nm_cst,
              het_nm_onf, het_ud_abs, het_ud_cst, het_ud_onf, het_rc_abs, 
              het_rc_cst, het_rc_onf)

hom_vals <- c(hom_dm_abs, hom_dm_cst, hom_dm_onf, hom_nm_abs, hom_nm_cst,
              hom_nm_onf, hom_ud_abs, hom_ud_cst, hom_ud_onf, hom_rc_abs, 
              hom_rc_cst, hom_rc_onf)

pops <- c(abs.tag, cst.tag, onf.tag, abs.tag, cst.tag, onf.tag,
          abs.tag, cst.tag, onf.tag, abs.tag, cst.tag, onf.tag)

migs <- c(dm.tag.abs, dm.tag.cst, dm.tag.onf, nm.tag.abs, nm.tag.cst,
          nm.tag.onf, ud.tag.abs, ud.tag.cst, ud.tag.onf, rc.tag.abs,
          rc.tag.cst, rc.tag.onf)

bm <- c(dm.tag.abs, nm.tag.abs, ud.tag.abs, rc.tag.abs)

het.df <- cbind.data.frame(x=pops, y=het_vals, m=migs)
het.df$htm <- gstat_vectorize(get_mu(het_vals, pops, bm))
het.df$htse <- gstat_vectorize(get_se(het_vals, pops, bm))

hom.df <- cbind.data.frame(x=pops, y=hom_vals, m=migs)
hom.df$hm <- gstat_vectorize(get_mu(hom_vals, pops, bm))
hom.df$hse <- gstat_vectorize(get_se(hom_vals, pops, bm))
het.plot <- ggplot(het.df, aes(x=factor(x), y=y, color=factor(m)))+
  geom_point(position=position_dodge(width=0.5), 
             alpha=0.2)+ xlab("Population") + 
  ylab("Number of Heterozygous Sites") +
  ggtitle("Number Heterozygous Sites per Individual (m = 0.05)") +
  labs(color="Migration Scheme") + mytheme +
  stat_summary(position=position_dodge(width=0.5),fun=mean, geom="point", 
               shape=20, size=8, alpha=0.2) +
  geom_errorbar(aes(ymin=htm-htse,
                    ymax=htm+htse),
                position=position_dodge(.5), width=0.4)

het.plot
ggsave("het_box.png", het.plot)

hom.plot <- ggplot(hom.df, aes(x=factor(x), y=y, color=factor(m)))+
  xlab("Population") + 
  geom_point(position=position_dodge(width=0.5), 
             alpha=0.2)+ xlab("Population") + 
  ylab("Number of Homozygous Sites") +
  ggtitle("Number of Homozygous Sites per Individual (m = 0.05)") +
  labs(color="Migration Scheme") + mytheme+
  stat_summary(position=position_dodge(width=0.5),fun=mean, geom="point", 
               shape=20, size=8, alpha=0.2) +
  geom_errorbar(aes(ymin=hm-hse,
                    ymax=hm+hse),
                position=position_dodge(.5), width=0.4)

hom.plot
ggsave("hom_box.png", hom.plot)
  
