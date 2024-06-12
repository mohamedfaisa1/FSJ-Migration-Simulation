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

empiricalSFS <- read.table(file="~/Downloads/empiricalNeutralSFS.txt", 
                           header=T)
inds_sav <- which(empiricalSFS$Pop == "SAV")
empiricalSFS$Pop[inds_sav] <- "CST"

res_emp <- tapply(empiricalSFS$AF, empiricalSFS$Pop, function(x) return(x))

odf <- data.frame(onf=res_emp$ONF)
adf <- data.frame(abs=res_emp$ABS)
csdf <- data.frame(cst=res_emp$CST)

e_onf.plot <- ggplot(odf, aes(x=onf)) +
  geom_histogram(binwidth=0.065, bins=100, fill="red", 
                 colour="black", aes(y=..count../sum(..count..))) + 
  mytheme + ggtitle("ONF") + 
  xlab("Derived Allele Frequency") +
  ylab("Proportion") +
  geom_vline(xintercept=mean(odf$onf), colour="black", linewidth=0.8) +
  annotate("text", x=0.25, y=0.25,
           label=paste("AF = ",
             as.character(round(mean(odf$onf),3))),
           color="black")

e_onf.plot 

e_abs.plot <- ggplot(adf, aes(x=abs)) +
  geom_histogram(binwidth=0.0955, bins=100, fill="lightblue3", 
                 colour="black", aes(y=..count../sum(..count..))) +
  mytheme + ggtitle("ABS") + 
  xlab("Derived Allele Frequency") +
  ylab("Proportion") +
  geom_vline(xintercept=mean(adf$abs), colour="black", linewidth=0.8) +
  annotate("text", x=0.25, y=0.25,
           label=paste("AF = ",
                       as.character(round(mean(adf$abs),3))),
           color="black")

e_abs.plot 

e_cst.plot <- ggplot(csdf, aes(x=cst)) +
  geom_histogram(binwidth=0.0355, fill="magenta2", 
                 colour="black", aes(y=..count../sum(..count..))) +
  mytheme + ggtitle("CST") +
  xlab("Derived Allele Frequency") +
  ylab("Proportion") +
  geom_vline(xintercept=mean(csdf$cst), colour="black", linewidth=0.8)+
  annotate("text", x=0.35, y=0.14,
           label=paste("AF = ",
                       as.character(round(mean(csdf$cst),3))),
           color="black")

e_cst.plot 

emp_p <-  ggpubr::ggarrange(e_onf.plot, e_abs.plot, e_cst.plot, 
                                 labels=c("A", "B", "C"), 
                                 ncol=1, nrow=3, legend="none") 
emp_p <- ggpubr::annotate_figure(emp_p,
                              top=ggpubr::text_grob("Empirical SFS", 
                                                    color="black", 
                                                    face="bold", 
                                                    size=18))
emp_p

ggsave("emp_sfs.png", emp_p)
