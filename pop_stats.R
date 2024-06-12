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

conf95 <- function(x, d=F){
  if(isTRUE(d)){x <- x/100}
  return(((1.96/sqrt(100))*sd(x)))
}

nm.tag <- rep("No Migration", 100)
ud.tag <- rep("Unidirectional Migration", 100)
rc.tag <- rep("Reciprocal Migration", 100)
dm.tag <- rep("Directional Migration", 100)

mig_tags <- c(nm.tag, ud.tag, rc.tag, dm.tag)
mig.tags <- c(mig_tags, mig_tags, mig_tags)

pi_nm <- np$load("pi_no_mig_0.05.npy")
pi_ud <- np$load("pi_unidirectional_0.05.npy")
pi_rc <- np$load("pi_reciprocal_0.05.npy")
pi_dm <- np$load("pi_directional_0.05.npy")

pi_onf <- c(pi_nm[1,], pi_ud[1,], pi_rc[1,], pi_dm[1,])
pi_abs <- c(pi_nm[2,], pi_ud[2,], pi_rc[2,], pi_dm[2,])
pi_cst <- c(pi_nm[3,], pi_ud[3,], pi_rc[3,], pi_dm[3,])


tD_nm <- np$load("tajimaD_no_mig_0.05.npy")
tD_ud <- np$load("tajimaD_unidirectional_0.05.npy")
tD_rc <- np$load("tajimaD_reciprocal_0.05.npy")
tD_dm <- np$load("tajimaD_directional_0.05.npy")

tD_onf <- c(tD_nm[1,], tD_ud[1,], tD_rc[1,], tD_dm[1,])
tD_abs <- c(tD_nm[2,], tD_ud[2,], tD_rc[2,], tD_dm[2,])
tD_cst <- c(tD_nm[3,], tD_ud[3,], tD_rc[3,], tD_dm[3,])


sfs_onf <- np$load("sfs_onf_0.05.npy")
mu_onf <- apply(sfs_onf, 1, function(x) return(mean(x)/100))
ci_onf <- apply(sfs_onf, 1, function(x) return(conf95(x, T)))
sfs_onf <- apply(sfs_onf, 1, function(x) x/sum(x))
sfs_onf <- c(sfs_onf[,1], sfs_onf[,2], sfs_onf[,3], sfs_onf[,4])

sfs_abs <- np$load("sfs_abs_0.05.npy")
mu_abs <- apply(sfs_abs, 1, function(x) return(mean(x)/100))
ci_abs <- apply(sfs_abs, 1, function(x) return(conf95(x, T)))
sfs_abs <- apply(sfs_abs, 1, function(x) x/sum(x))
sfs_abs <- c(sfs_abs[,1], sfs_abs[,2], sfs_abs[,3], sfs_abs[,4])

sfs_cst <- np$load("sfs_cst_0.05.npy")
mu_cst <- apply(sfs_cst, 1, function(x) return(mean(x)/100))
ci_cst <- apply(sfs_cst, 1, function(x) return(conf95(x, T)))
sfs_cst <- apply(sfs_cst, 1, function(x) x/sum(x))
sfs_cst <- c(sfs_cst[,1], sfs_cst[,2], sfs_cst[,3], sfs_cst[,4])

sfs.df <- cbind.data.frame(onf=sfs_onf, abs=sfs_abs, cst=sfs_cst,
                           mt=mig_tags, x=rep((1:100)/100, 4))

cst.tag <- rep("CST", 400)
abs.tag <- rep("ABS", 400)
onf.tag <- rep("ONF", 400)
pop.tags <- c(onf.tag, abs.tag, cst.tag)


pi_vals <- c(pi_onf, pi_abs, pi_cst) 
tD_vals <- c(tD_onf, tD_abs, tD_cst) 


pop_stats.df <- cbind.data.frame(pv=pi_vals, tv=tD_vals, 
                                 ptag=pop.tags, mtag=mig.tags)

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
  p <- c(p1,p3,p2)
  return(p)
}

pop_stats.df$pm <- gstat_vectorize(get_mu(pi_vals, pop.tags, mig_tags))
pop_stats.df$tm <- gstat_vectorize(get_mu(tD_vals, pop.tags, mig_tags))
pop_stats.df$pse <- gstat_vectorize(get_se(pi_vals, pop.tags, mig_tags))
pop_stats.df$tse <- gstat_vectorize(get_se(tD_vals, pop.tags, mig_tags))

pi.plot <- ggplot(pop_stats.df, aes(x=factor(ptag), y=pv,
                                    color=factor(mtag)))+
  geom_point(position=position_dodge(width=0.5), 
             alpha=0.2)+ xlab("Population") + 
  ylab(bquote(pi)) +
  ggtitle(bquote("Genetic Diversity"~(pi)~"per Population (m = 0.05)")) +
  labs(color="Migration Scheme") + mytheme +
  stat_summary(position=position_dodge(width=0.5),fun=mean, geom="point", 
               shape=20, size=8, alpha=0.2) +
  geom_errorbar(aes(ymin=pm-pse,
                    ymax=pm+pse),
                position=position_dodge(.5), width=0.4)

pi.plot
ggsave("pi_box.png",pi.plot)


tD.plot <- ggplot(pop_stats.df, aes(x=factor(ptag), y=tv,
                                    color=factor(mtag)))+
  geom_point(position=position_dodge(width=0.5), 
             alpha=0.2)+ xlab("Population") + 
  ylab(bquote(D)) +
  ggtitle(bquote("Tajima's D"~(D)~"per Population (m = 0.05)")) +
  labs(color="Migration Scheme") + mytheme +
  stat_summary(position=position_dodge(width=0.5),fun=mean, geom="point", 
               shape=20, size=8, alpha=0.2) +
  geom_errorbar(aes(ymin=tm-tse,
                    ymax=tm+tse),
                position=position_dodge(.5), width=0.4)


tD.plot
ggsave("tD_box.png", tD.plot)

# pid: 1 -> onf, 2 -> abs, 3 -> cst
plot_sfs_mig <- function(df, pid, muv, civ){
  pops <- c("ONF", "ABS", "CST")
  p <- pops[pid]
  dat <- tapply(df[,pid], df[,4], function(x) return(x))
  tle <- paste("Simulated",p,"SFS (m = 0.05)")
  ndf <- cbind.data.frame(dat[[2]], dat[[4]], 
                          dat[[3]], dat[[1]],
                          (1:100)/100)
  colnames(ndf) <- c("nm", "ud", "rc", "dm", "x")
  
  p1 <- ggplot(ndf, aes(x=x, y=nm, fill=x)) +
    geom_bar(stat="identity", position="dodge") + 
    mytheme + ggtitle("No Migration") +
    xlab("Derived Allele Frequency") +
    ylab("Proportion") +
    geom_vline(xintercept=muv[1], colour="red", linewidth=0.8) +
    geom_rect(aes(xmin=muv[1]-civ[1], xmax=muv[1]+civ[1]), 
                  fill="grey90", alpha=0.05,
                  ymin=-Inf, ymax=1) +
    scale_fill_viridis_c() +
    annotate("text", x=0.4, y=(0.9*ndf$nm[1]),
             label=paste0("AF = ",as.character(round(muv[1],3)),"±",
                          as.character(round(civ[1],3))),
             color="red2")
  
    
  p2 <- ggplot(ndf, aes(x=x, y=ud, fill=x)) +
    geom_bar(stat="identity", position="dodge") + 
    mytheme + ggtitle("Unidirectional Migration") +
    xlab("Derived Allele Frequency") +
    ylab("Proportion") +
    geom_vline(xintercept=muv[2], colour="red", linewidth=0.8) +
    geom_rect(aes(xmin=muv[2]-civ[2], xmax=muv[2]+civ[2]), 
                  fill="grey90", alpha=0.05,
                  ymin=-Inf, ymax=1) +
    scale_fill_viridis_c() +
    annotate("text", x=0.4, y=(0.9*ndf$ud[1]),
             label=paste0("AF = ",as.character(round(muv[2],3)),"±",
                          as.character(round(civ[2],3))),
             color="red2")
  
  p3 <- ggplot(ndf, aes(x=x, y=rc, fill=x)) +
    geom_bar(stat="identity", position="dodge") + 
    mytheme + ggtitle("Reciprocal Migration") +
    xlab("Derived Allele Frequency") +
    ylab("Proportion") + 
    geom_vline(xintercept=muv[3], colour="red", linewidth=0.8) +
    geom_rect(aes(xmin=muv[3]-civ[3], xmax=muv[3]+civ[3]), 
              fill="grey90", alpha=0.05,
              ymin=-Inf, ymax=1) +
    scale_fill_viridis_c() +
    annotate("text", x=0.4, y=(0.9*ndf$rc[1]),
             label=paste0("AF = ",as.character(round(muv[3],3)),"±",
                          as.character(round(civ[3],3))),
             color="red2")
  
  p4 <- ggplot(ndf, aes(x=x, y=dm, fill=x)) +
    geom_bar(stat="identity", position="dodge") + 
    mytheme + ggtitle("Directional Migration") +
    xlab("Derived Allele Frequency") +
    ylab("Proportion") + 
    geom_vline(xintercept=muv[4], colour="red", linewidth=0.8) +
    geom_rect(aes(xmin=muv[4]-civ[4], xmax=muv[4]+civ[4]), 
                  fill= "grey90", alpha=0.05,
                  ymin=-Inf, ymax=1) +
    scale_fill_viridis_c() +
    annotate("text", x=0.4, y=(0.9*ndf$dm[1]),
             label=paste0("AF = ",as.character(round(muv[4],3)),"±",
                          as.character(round(civ[4],3))),
             color="red2")
  
  
  fp <- ggpubr::ggarrange(p1,p2,p3,p4, 
                          ncol=2, nrow=2, legend="none") 
  fp <- ggpubr::annotate_figure(fp,
                                top=ggpubr::text_grob(tle, 
                                color="black", 
                                face="bold", 
                                size=18))
  return(fp)
}

o5 <- plot_sfs_mig(sfs.df, 1, mu_onf, ci_onf)
a5 <-  plot_sfs_mig(sfs.df, 2, mu_abs, ci_abs)
c5 <- plot_sfs_mig(sfs.df, 3, mu_cst, ci_cst)



fall <- ggpubr::ggarrange(pi.plot,het.plot,pi1.plot,het1.plot, labels=LETTERS[1:4],
                          ncol=2, nrow=2)

small <- ggpubr::ggarrange(fst.plot,hom.plot,fst1.plot,hom1.plot, labels=LETTERS[1:4],
                           ncol=2, nrow=2)


joint5 <- ggpubr::ggarrange(o5,a5,c5,tD.plot, labels=LETTERS[1:4],
                            ncol=2, nrow=2)

ggsave("j5.png", joint5)

onf.plot <- ggplot(sfs.df, aes(x=x, y=onf, fill=mt)) +
  geom_bar(stat="identity", position="dodge") + 
  mytheme + ggtitle("ONF SFS (m=0.05)") + xlab("Derived Allele Frequency") +
  ylab("Proportion") + scale_color_brewer(palette = "PuOr")

onf.plot 
ggsave("onf_sfs.png", onf.plot)

abs.plot <- ggplot(sfs.df, aes(x=x, y=abs, fill=mt)) +
  geom_bar(stat="identity", position="dodge") + 
  mytheme + ggtitle("ABS SFS (m=0.05)") + xlab("Derived Allele Frequency") +
  ylab("Proportion")

abs.plot 
ggsave("abs_sfs.png", abs.plot)

cst.plot <- ggplot(sfs.df, aes(x=x, y=cst, fill=mt)) +
  geom_bar(stat="identity", position="dodge") +
  mytheme + ggtitle("CST SFS (m=0.05)") + xlab("Derived Allele Frequency") +
  ylab("Proportion")

cst.plot 
ggsave("cst_sfs.png", cst.plot)



help_density <- function(x, sname){
  df <- data.frame(x=unlist(x))
  p <- ggplot(df, aes(x = x)) +
      geom_density(color = "red", fill = "blue", alpha = 0.328) + 
      mytheme + ylab("Density") + xlab("Value") + 
      scale_x_continuous(breaks=median(df$x)[1],
                         labels=round(mean(df$x),7))
  return(p)
}


plot_density <- function(v, p, m, sname, migr=0.05){
  res <- tapply(v, p, 
         function(x) tapply(x, m, function(w) return(w)))
  
  p1 <- help_density(res[[1]][1], sname) 
  p1 <- p1 + ggtitle(paste("ABS", "DM"))
  p2 <- help_density(res[[1]][2], sname) 
  p2 <- p2 + ggtitle(paste("ABS", "NM"))
  p3 <- help_density(res[[1]][3], sname)
  p3 <- p3 + ggtitle(paste("ABS", "RC"))
  p4 <- help_density(res[[1]][4], sname) 
  p4 <- p4 + ggtitle(paste("ABS","UD"))
  
  p5 <- help_density(res[[2]][1], sname) 
  p5 <- p5 + ggtitle(paste("CST", "DM"))
  p6 <- help_density(res[[2]][2], sname) 
  p6 <- p6 + ggtitle(paste("CST", "NM"))
  p7 <- help_density(res[[2]][3], sname)
  p7 <- p7 + ggtitle(paste("CST", "RC"))
  p8 <- help_density(res[[2]][4], sname) 
  p8 <- p8 + ggtitle(paste("CST","UD"))
  
  p9 <- help_density(res[[3]][1], sname) 
  p9 <- p9 + ggtitle(paste("ONF", "DM"))
  p10 <- help_density(res[[3]][2], sname) 
  p10 <- p10 + ggtitle(paste("ONF","NM"))
  p11 <- help_density(res[[3]][3], sname)
  p11 <- p11 + ggtitle(paste("ONF", "RC"))
  p12 <- help_density(res[[3]][4], sname) 
  p12 <- p12 + ggtitle(paste("ONF", "UD"))
  
  
  p <- ggpubr::ggarrange(p10, p12, p11, p9,
                         p2, p4, p3, p1,
                         p6, p8, p7, p5,
                         labels=LETTERS[1:12], 
                         ncol=4, nrow=3, legend="none")
  p <- ggpubr::annotate_figure(p,
                               top=ggpubr::text_grob(paste(sname, 
                                                           "Density Plots",
                                                           "(m = ", migr, ")"), 
                                                     color="black", 
                                                     face="bold", 
                                                     size=18))
  return(p)
}


fst_plot_density <- function(v, p, m, sname="Fst", migr=0.05){
  res <- tapply(v, p, 
                function(x) tapply(x, m, function(w) return(w)))
  
  p1 <- help_density(res[[1]][1], sname) 
  p1 <- p1 + ggtitle(paste("ABS/CST", "DM"))
  p2 <- help_density(res[[1]][2], sname) 
  p2 <- p2 + ggtitle(paste("ABS/CST", "NM"))
  p3 <- help_density(res[[1]][3], sname)
  p3 <- p3 + ggtitle(paste("ABS/CST", "RC"))
  p4 <- help_density(res[[1]][4], sname) 
  p4 <- p4 + ggtitle(paste("ABS","UD"))
  
  p5 <- help_density(res[[2]][1], sname) 
  p5 <- p5 + ggtitle(paste("ONF/ABS", "DM"))
  p6 <- help_density(res[[2]][2], sname) 
  p6 <- p6 + ggtitle(paste("ONF/ABS", "NM"))
  p7 <- help_density(res[[2]][3], sname)
  p7 <- p7 + ggtitle(paste("ONF/ABS", "RC"))
  p8 <- help_density(res[[2]][4], sname) 
  p8 <- p8 + ggtitle(paste("ONF/ABS","UD"))
  
  p9 <- help_density(res[[3]][1], sname) 
  p9 <- p9 + ggtitle(paste("ONF/CST", "DM"))
  p10 <- help_density(res[[3]][2], sname) 
  p10 <- p10 + ggtitle(paste("ONF/CST","NM"))
  p11 <- help_density(res[[3]][3], sname)
  p11 <- p11 + ggtitle(paste("ONF/CST", "RC"))
  p12 <- help_density(res[[3]][4], sname) 
  p12 <- p12 + ggtitle(paste("ONF/CST", "UD"))
  
  
  p <- ggpubr::ggarrange(p10, p12, p11, p9,
                         p2, p4, p3, p1,
                         p6, p8, p7, p5,
                         labels=LETTERS[1:12], 
                         ncol=4, nrow=3, legend="none")
  p <- ggpubr::annotate_figure(p,
                               top=ggpubr::text_grob(paste(sname, 
                                                           "Density Plots",
                                                           "(m = ", migr, ")"), 
                                                     color="black", 
                                                     face="bold", 
                                                     size=18))
  return(p)
}

ggsave("pden5.png", plot_density(pi_vals, pop.tags, 
                                 mig_tags, "Genetic Diversity"))

ggsave("tden5.png", plot_density(tD_vals, pop.tags, 
                                 mig_tags, "Tajima's D"))


do_shapiro <- function(v, p, m){
  res <- tapply(v, p, 
                function(x) tapply(x, m, function(w) shapiro.test(w)$p.value))
  return(res)
}


comp_to_no_mig <- function(v, p, m){
  res <- tapply(v, p, 
                function(x) tapply(x, m, function(w) return(w)))
  
  nmabs <- unlist(res[[1]][2])
  nmcst <- unlist(res[[2]][2])
  nmonf <- unlist(res[[3]][2])
  
  res[[1]][2] <- NULL
  res[[2]][2] <- NULL
  res[[3]][2] <- NULL
  
  op <- lapply(res[[3]], 
               function(x) t.test(unlist(x), nmonf, paired=T)$p.value)
  ap <- lapply(res[[1]], 
               function(x) t.test(unlist(x), nmabs, paired=T)$p.value)
  cp <- lapply(res[[2]], 
               function(x) t.test(unlist(x), nmcst, paired=T)$p.value)
  
  pvs <- list(op, ap, cp)
  
  return(pvs)  
}
