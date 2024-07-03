# set working directory
setwd("C:/master/SommerSemester_24/RPI Project/RPI_Project_Christian/data")

# get the raw data
data <- read.csv("FruitDatabase.shared_Dehiscence_n_FvD_Unq.Sp_tidy.FAU.sortnOrd.csv")
unique(data$John.Dehiscence.n.FvsD)

# subset data focal columns (species and dehiscence)
dat <- data[, c("Genus_Species_Accepted", "John.Dehiscence.n.FvsD")]

# load necessary library to transform dehiscence column
library(tidyr)
library(dplyr)

# separate the traits into two columns
dat1 <- dat %>%
  separate(John.Dehiscence.n.FvsD, into = c("dehiscence","type"), sep = "_")


# set the first column as rownames
rownames(dat1) <- dat1$Genus_Species_Accepted

# remove the first column
df <- dat1[, -1]
table(df)

# phylogeny
#install.packages("devtools")
library("devtools")
#Sys.setenv(GITHUB_PAT = "...") # use your github_pat
#install_github("jinyizju/V.PhyloMaker2")
library("V.PhyloMaker2")
library(phytools)


# get phylogeny
phy <- GBOTB.extended.LCVP

#######################################
# ... feeding the tree with your data
# Prepare species_data dataframe
# 
species_data <- data.frame(species = data$Genus_Species_Accepted, 
                           genus = data$Genus_Accepted, family = data$Family, 
                           order = data$Order)

# Omit rows with NA values
# we don't have to do that, there aren't any NAs
# sp.df <- na.omit(species_data)

# feed the phylogeny with your data
# phy2 <- phylo.maker(sp.list = species_data, tree = GBOTB.extended.LCVP, 
#                    scenarios = "S3")
# scenario 3 doesn't work, only scenario 1 works for me
phy2 <- phylo.maker(sp.list = species_data, tree = GBOTB.extended.LCVP, 
                    scenarios = "S1")

saveRDS(phy2, file = "phylogeny.rds")
phy2 <- phylogeny


#prune the tree and data 
library(geiger)
treedat <- treedata(phy2$scenario.1, df, sort=FALSE, warnings=TRUE)

# data after merging | 8713 species left
df1 <- treedat$data

x <- df1[,1] # dehiscence
names(x) <- rownames(df1)
y <- df1[,2]
names(y) <- rownames(df1)


tree <- phy2$scenario.1

###############################
# ancestral state reconstruction (phylo tree)
windows()
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(x[tree$tip.label],c("dehiscent","indehiscent")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
legend(fill=c("blue","red", "darkviolet", "darkgreen"),legend=c("dehiscent","indehiscent", "dry", "fleshy"),
       x=35,y=1750)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(y[tree$tip.label],c("dry","fleshy")),piecol=c("darkviolet","darkgreen"),
          cex=0.3)

################################

# models analysing independent & dependent models (first both traits depend on
# each other, than only one is the dependent variable)
# both depend gives us the best fitting model, this is also used
# in the report

library(phytools)
tree <- treedat$phy

fit.new <- fitPagel(tree, x, y)
fit.new

saveRDS(fit.new, file = "fruit_best_fitting_model.rds")
fit.new <- fruit_best_fitting_model
# Pagel's binary character correlation test:

# Assumes "ARD" substitution model for both characters

# Independent model rate matrix:
#                    dehiscent|dry dehiscent|fleshy indehiscent|dry indehiscent|fleshy
# dehiscent|dry          -0.008908         0.003668        0.005241           0.000000
# dehiscent|fleshy        0.002160        -0.007401        0.000000           0.005241
# indehiscent|dry         0.002073         0.000000       -0.005741           0.003668
# indehiscent|fleshy      0.000000         0.002073        0.002160          -0.004234

# Dependent (x & y) model rate matrix:
#                    dehiscent|dry dehiscent|fleshy indehiscent|dry indehiscent|fleshy
# dehiscent|dry          -0.006455         0.004203        0.002251           0.000000
# dehiscent|fleshy        0.062905        -0.147167        0.000000           0.084262
# indehiscent|dry         0.002904         0.000000       -0.006878           0.003973
# indehiscent|fleshy      0.000000         0.000750        0.000784          -0.001533

# Model fit:
#             log-likelihood      AIC
# independent      -2382.695 4773.390
# dependent        -1921.784 3859.568

# Hypothesis test result:
#   likelihood-ratio:  921.823 
#   p-value:  3.1139e-198 

# Model fitting method used was fitMk 

windows()
plot(fit.new, lwd.by.rate=TRUE)


# dependencies, dehiscence depends on type or type depends on dehiscence?
# 
fit.x <- fitPagel(tree, x, y, dep.var="x")
fit.x
saveRDS(fit.x, file = "fruit_depVarDehiscence_model.rds")
#Pagel's binary character correlation test:

#Assumes "ARD" substitution model for both characters

#Independent model rate matrix:
#                   dehiscent|dry dehiscent|fleshy indehiscent|dry indehiscent|fleshy
#dehiscent|dry          -0.008908         0.003668        0.005241           0.000000
#dehiscent|fleshy        0.002160        -0.007401        0.000000           0.005241
#indehiscent|dry         0.002073         0.000000       -0.005741           0.003668
#indehiscent|fleshy      0.000000         0.002073        0.002160          -0.004234

#Dependent (x only) model rate matrix:
#                   dehiscent|dry dehiscent|fleshy indehiscent|dry indehiscent|fleshy
#dehiscent|dry          -0.005866         0.003835        0.002030           0.000000
#dehiscent|fleshy        0.001773        -0.131338        0.000000           0.129565
#indehiscent|dry         0.004435         0.000000       -0.008270           0.003835
#indehiscent|fleshy      0.000000         0.000000        0.001773          -0.001773

#Model fit:
#            log-likelihood      AIC
#independent      -2382.695 4773.390
#dependent        -1972.906 3957.813

#Hypothesis test result:
#  likelihood-ratio:  819.577 
#  p-value:  1.07401e-178 

#Model fitting method used was fitMk

windows()
plot(fit.x, lwd.by.rate=TRUE)

#
fit.y <- fitPagel(tree, x, y, dep.var="y")
fit.y

saveRDS(fit.y, file = "fruit_depVarType_model.rds")
#Pagel's binary character correlation test:

#Assumes "ARD" substitution model for both characters

#Independent model rate matrix:
#                   dehiscent|dry dehiscent|fleshy indehiscent|dry indehiscent|fleshy
#dehiscent|dry          -0.008908         0.003668        0.005241           0.000000
#dehiscent|fleshy        0.002160        -0.007401        0.000000           0.005241
#indehiscent|dry         0.002073         0.000000       -0.005741           0.003668
#indehiscent|fleshy      0.000000         0.002073        0.002160          -0.004234

#Dependent (y only) model rate matrix:
#                   dehiscent|dry dehiscent|fleshy indehiscent|dry indehiscent|fleshy
#dehiscent|dry          -0.005739         0.000590        0.005149           0.000000
#dehiscent|fleshy        0.101064        -0.106213        0.000000           0.005149
#indehiscent|dry         0.003052         0.000000       -0.019129           0.016077
#indehiscent|fleshy      0.000000         0.003052        0.000311          -0.003363

#Model fit:
#            log-likelihood      AIC
#independent      -2382.695 4773.390
#dependent        -2093.463 4198.925

#Hypothesis test result:
#  likelihood-ratio:  578.465 
#  p-value:  2.44279e-126 

#Model fitting method used was fitMk


windows()
plot(fit.y, lwd.by.rate=TRUE)


aic <- setNames(c(fit.new$independent.AIC,
                  fit.x$dependent.AIC,
                  fit.y$dependent.AIC,
                  fit.new$dependent.AIC),
                c("independent","dependent x",
                  "dependent y","dependent x&y"))
aic
# independent   dependent x   dependent y dependent x&y 
# 4773.390      3957.813      4198.925      3859.568 

