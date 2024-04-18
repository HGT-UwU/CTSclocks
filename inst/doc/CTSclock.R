## ----echo=FALSE---------------------------------------------------------------
knitr::opts_chunk$set(crop = NULL)

## ----message=FALSE------------------------------------------------------------
library(CTSclocks)

data(MurphyGSE88890)
agePred.v = CTSclockAge(beta.m, CTSclock = 'Neu-In', dataType = 'bulk', CTF.m = NULL, tissue = 'brain')

## ----echo=FALSE---------------------------------------------------------------
library(ggplot2)
phenotype.df$DNAmAgePred = as.numeric(agePred.v)
## Make a scatter plot for DNAm age against chronological age
PvalFormat = function(PvalRaw){
  PvalRaw = as.numeric(PvalRaw)
  if (PvalRaw < 1e-300){
    return('P<1e-300')
  }else{
    if (PvalRaw >= 0.01 & PvalRaw <= 1){
      return(paste0('P=',sprintf("%0.3f", PvalRaw)))
    }else if(PvalRaw < 0.01 & PvalRaw >= 1e-300){
      return(paste0('P=', format(PvalRaw, digits = 1, scientific = T)))
    }
  }
}
p = ggplot(phenotype.df, aes(x = Age, y = DNAmAgePred)) + geom_point(shape = 20, color = 'black') + theme_bw() +
    theme(axis.title.x = element_text(vjust = 2), axis.title.y = element_text(vjust = 2),
          axis.text.x = element_text(vjust = 2), axis.text.y = element_text(vjust = 2),
          plot.title = element_text(hjust = 0.5),
          legend.background = element_rect(linetype = "solid", colour = "black",linewidth = 0.2)) + 
  ggtitle("Neu-In_Murphy(n=38)") + xlab("Chronological Age") + ylab("Predicted Age") + xlim(0, 100) + ylim(0, 100)
corValue = cor.test(phenotype.df$Age, phenotype.df$DNAmAgePred)$estimate
pValue = cor.test(phenotype.df$Age, phenotype.df$DNAmAgePred)$p.value
p = p + annotate("text", x = 20, y = 65, label=paste("PCC=", sprintf("%0.2f", corValue), sep = ""), color = 'darkred', size = 5, fontface = 'italic') + annotate("text", x = 20, y = 85, label = PvalFormat(pValue), color = 'darkred', size = 5, fontface = 'italic')
lm.o = lm(phenotype.df$DNAmAgePred ~ phenotype.df$Age)
a.v = summary(lm.o)$coeff[, 1]
p = p + geom_abline(intercept = 0, slope = 1, color = 'black', linetype = 'dashed', lwd = 0.2) + geom_abline(intercept = a.v[1], slope = a.v[2], color = '#650023', linetype = 'dashed')
print(p)

## ----message=FALSE------------------------------------------------------------
library(CTSclocks)

data(PaiGSE112179)
agePred.v = CTSclockAge(beta.m, CTSclock = 'Neu-In', dataType = 'sorted', CTF.m = NULL, tissue = 'brain')

## ----echo=FALSE---------------------------------------------------------------
library(ggplot2)
phenotype.df$DNAmAgePred = as.numeric(agePred.v)
## Make a scatter plot for DNAm age against chronological age
PvalFormat = function(PvalRaw){
  PvalRaw = as.numeric(PvalRaw)
  if (PvalRaw < 1e-300){
    return('P<1e-300')
  }else{
    if (PvalRaw >= 0.01 & PvalRaw <= 1){
      return(paste0('P=',sprintf("%0.3f", PvalRaw)))
    }else if(PvalRaw < 0.01 & PvalRaw >= 1e-300){
      return(paste0('P=', format(PvalRaw, digits = 1, scientific = T)))
    }
  }
}
p = ggplot(phenotype.df, aes(x = Age, y = DNAmAgePred)) + geom_point(shape = 20, color = 'black') + theme_bw() +
    theme(axis.title.x = element_text(vjust = 2), axis.title.y = element_text(vjust = 2),
          axis.text.x = element_text(vjust = 2), axis.text.y = element_text(vjust = 2),
          plot.title = element_text(hjust = 0.5),
          legend.background = element_rect(linetype = "solid", colour = "black",linewidth = 0.2)) + 
  ggtitle("Neu-In_Murphy(n=38)") + xlab("Chronological Age") + ylab("Predicted Age") + xlim(0, 100) + ylim(0, 100)
corValue = cor.test(phenotype.df$Age, phenotype.df$DNAmAgePred)$estimate
pValue = cor.test(phenotype.df$Age, phenotype.df$DNAmAgePred)$p.value
p = p + annotate("text", x = 20, y = 65, label=paste("PCC=", sprintf("%0.2f", corValue), sep = ""), color = 'darkred', size = 5, fontface = 'italic') + annotate("text", x = 20, y = 85, label = PvalFormat(pValue), color = 'darkred', size = 5, fontface = 'italic')
lm.o = lm(phenotype.df$DNAmAgePred ~ phenotype.df$Age)
a.v = summary(lm.o)$coeff[, 1]
p = p + geom_abline(intercept = 0, slope = 1, color = 'black', linetype = 'dashed', lwd = 0.2) + geom_abline(intercept = a.v[1], slope = a.v[2], color = '#650023', linetype = 'dashed')
print(p)

## ----message=FALSE------------------------------------------------------------
library(CTSclocks)

data(MurphyGSE88890)
agePred.v = CTSclockAge(beta.m, CTSclock = 'Neu-Ex', dataType = 'bulk', CTF.m = NULL, tissue = 'brain')

## ----echo=FALSE---------------------------------------------------------------
library(ggplot2)
phenotype.df$DNAmAgePred = as.numeric(agePred.v)
## Make a scatter plot for DNAm age against chronological age
PvalFormat = function(PvalRaw){
  PvalRaw = as.numeric(PvalRaw)
  if (PvalRaw < 1e-300){
    return('P<1e-300')
  }else{
    if (PvalRaw >= 0.01 & PvalRaw <= 1){
      return(paste0('P=',sprintf("%0.3f", PvalRaw)))
    }else if(PvalRaw < 0.01 & PvalRaw >= 1e-300){
      return(paste0('P=', format(PvalRaw, digits = 1, scientific = T)))
    }
  }
}
p = ggplot(phenotype.df, aes(x = Age, y = DNAmAgePred)) + geom_point(shape = 20, color = 'black') + theme_bw() +
    theme(axis.title.x = element_text(vjust = 2), axis.title.y = element_text(vjust = 2),
          axis.text.x = element_text(vjust = 2), axis.text.y = element_text(vjust = 2),
          plot.title = element_text(hjust = 0.5),
          legend.background = element_rect(linetype = "solid", colour = "black",linewidth = 0.2)) + 
  ggtitle("Neu-In_Murphy(n=38)") + xlab("Chronological Age") + ylab("Predicted Age") + xlim(0, 100) + ylim(0, 100)
corValue = cor.test(phenotype.df$Age, phenotype.df$DNAmAgePred)$estimate
pValue = cor.test(phenotype.df$Age, phenotype.df$DNAmAgePred)$p.value
p = p + annotate("text", x = 20, y = 65, label=paste("PCC=", sprintf("%0.2f", corValue), sep = ""), color = 'darkred', size = 5, fontface = 'italic') + annotate("text", x = 20, y = 85, label = PvalFormat(pValue), color = 'darkred', size = 5, fontface = 'italic')
lm.o = lm(phenotype.df$DNAmAgePred ~ phenotype.df$Age)
a.v = summary(lm.o)$coeff[, 1]
p = p + geom_abline(intercept = 0, slope = 1, color = 'black', linetype = 'dashed', lwd = 0.2) + geom_abline(intercept = a.v[1], slope = a.v[2], color = '#650023', linetype = 'dashed')
print(p)

## ----message=FALSE------------------------------------------------------------
library(CTSclocks)

data(PaiGSE112179)
agePred.v = CTSclockAge(beta.m, CTSclock = 'Neu-Ex', dataType = 'sorted', CTF.m = NULL, tissue = 'brain')

## ----echo=FALSE---------------------------------------------------------------
library(ggplot2)
phenotype.df$DNAmAgePred = as.numeric(agePred.v)
## Make a scatter plot for DNAm age against chronological age
PvalFormat = function(PvalRaw){
  PvalRaw = as.numeric(PvalRaw)
  if (PvalRaw < 1e-300){
    return('P<1e-300')
  }else{
    if (PvalRaw >= 0.01 & PvalRaw <= 1){
      return(paste0('P=',sprintf("%0.3f", PvalRaw)))
    }else if(PvalRaw < 0.01 & PvalRaw >= 1e-300){
      return(paste0('P=', format(PvalRaw, digits = 1, scientific = T)))
    }
  }
}
p = ggplot(phenotype.df, aes(x = Age, y = DNAmAgePred)) + geom_point(shape = 20, color = 'black') + theme_bw() +
    theme(axis.title.x = element_text(vjust = 2), axis.title.y = element_text(vjust = 2),
          axis.text.x = element_text(vjust = 2), axis.text.y = element_text(vjust = 2),
          plot.title = element_text(hjust = 0.5),
          legend.background = element_rect(linetype = "solid", colour = "black",linewidth = 0.2)) + 
  ggtitle("Neu-In_Murphy(n=38)") + xlab("Chronological Age") + ylab("Predicted Age") + xlim(0, 100) + ylim(0, 100)
corValue = cor.test(phenotype.df$Age, phenotype.df$DNAmAgePred)$estimate
pValue = cor.test(phenotype.df$Age, phenotype.df$DNAmAgePred)$p.value
p = p + annotate("text", x = 20, y = 65, label=paste("PCC=", sprintf("%0.2f", corValue), sep = ""), color = 'darkred', size = 5, fontface = 'italic') + annotate("text", x = 20, y = 85, label = PvalFormat(pValue), color = 'darkred', size = 5, fontface = 'italic')
lm.o = lm(phenotype.df$DNAmAgePred ~ phenotype.df$Age)
a.v = summary(lm.o)$coeff[, 1]
p = p + geom_abline(intercept = 0, slope = 1, color = 'black', linetype = 'dashed', lwd = 0.2) + geom_abline(intercept = a.v[1], slope = a.v[2], color = '#650023', linetype = 'dashed')
print(p)

## ----message=FALSE------------------------------------------------------------
library(CTSclocks)

data(ExampleData_Liver)
agePred.v = CTSclockAge(Test.m, CTSclock = 'Hep', dataType = 'bulk', CTF.m = NULL, tissue = 'otherTissue')

## ----echo=FALSE, message=FALSE------------------------------------------------
library(ggpubr)
mae <- median(abs(Age - agePred.v))
# Create a data frame for ggplot2/ggpubr
data <- data.frame(ChronologicalAge = Age, PredictedAge = agePred.v)
ggscatter(
    data, x = "ChronologicalAge", y = "PredictedAge",
    add = "reg.line",                   # Remove automatic regression line
    xlab = "Chronological age (years)",
    ylab = "Predicted age (years)",
    color = "black", 
  ) + 
    geom_abline(slope = 1, intercept = 0, linetype="dashed", color="black") + 
    stat_cor(method = "pearson", label.x = 20, label.y = 90) +
    annotate("text",label=paste0("MedAE: ", format(mae, digits=2), " years"),
              x=30, y=85, size=4) +
    theme_minimal()+
    theme(panel.grid.major = element_blank(),                
          panel.grid.minor = element_blank(),                
          panel.background = element_blank(),                
          axis.line = element_line(),
          axis.ticks=element_line())+
    scale_x_continuous(breaks = seq(10, 90, by = 20), limits = c(10, 90)) +
    scale_y_continuous(breaks = seq(10, 90, by = 20), limits = c(10, 90))


## ----message=FALSE------------------------------------------------------------
library(CTSclocks)

data(ExampleData_Liver)
agePred.v = CTSclockAge(Test.m, CTSclock = 'Liver', dataType = 'bulk', CTF.m = NULL, tissue = 'otherTissue')

## ----echo=FALSE, message=FALSE------------------------------------------------
library(ggpubr)
mae <- median(abs(Age - agePred.v))
# Create a data frame for ggplot2/ggpubr
data <- data.frame(ChronologicalAge = Age, PredictedAge = agePred.v)
ggscatter(
    data, x = "ChronologicalAge", y = "PredictedAge",
    add = "reg.line",                   # Remove automatic regression line
    xlab = "Chronological age (years)",
    ylab = "Predicted age (years)",
    color = "black", 
  ) + 
    geom_abline(slope = 1, intercept = 0, linetype="dashed", color="black") + 
    stat_cor(method = "pearson", label.x = 20, label.y = 90) +
    annotate("text",label=paste0("MedAE: ", format(mae, digits=2), " years"),
              x=30, y=85, size=4) +
    theme_minimal()+
    theme(panel.grid.major = element_blank(),                
          panel.grid.minor = element_blank(),                
          panel.background = element_blank(),                
          axis.line = element_line(),
          axis.ticks=element_line())+
    scale_x_continuous(breaks = seq(10, 90, by = 20), limits = c(10, 90)) +
    scale_y_continuous(breaks = seq(10, 90, by = 20), limits = c(10, 90))


