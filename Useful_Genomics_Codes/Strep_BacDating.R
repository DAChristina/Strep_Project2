devtools::install_github("xavierdidelot/BactDating")

wd = "C:/Users/dac23/Documents/Downloads" # DIDE
wd = "C:/Users/dac23/Downloads" # library computers
wd = "/home/ron/Downloads" # personal OSs
setwd(wd)

library(tidyverse)
library(BactDating)
library(ape)
# see help(package='BactDating') for more info
# Time-scaled tree with BactDating
# https://xavierdidelot.github.io/BactDating/articles/Staph.html

# 1. Data wrangling ############################################################
library(data.table)
data <- fread("spn_uk_dude.csv")
data <- as.data.frame(data)
data <- data %>% 
  mutate(sample = as.character(sample))

tre <- loadGubbins("trial1")
tre_names <- as.data.frame(tre$tip.label) %>% 
  rename(ID = 'tre$tip.label')
tre_names <- left_join(tre_names, data,
                       by = c("ID" = "sample"))
tre_names <- tre_names %>% 
  mutate(Year = as.Date(Year))

# 2. BacDating #################################################################
# https://xavierdidelot.github.io/BactDating/articles/yourData.html
d = cbind(tre_names$Year, tre_names$Year+1) # 2-d matrix according to the articles above

set.seed(0)
res_pr=bactdate(tre,d,nbIts=1e5,
                showProgress = T)
plot(res_pr,'treeCI',show.tip.label = T)
plot(res_pr,'trace')

# Some info about model selection:
# https://taming-the-beast.org/tutorials/NS-tutorial/

# Further analysis
# https://xavierdidelot.github.io/BactDating/articles/yourData.html
# https://xavierdidelot.github.io/BactDating/articles/Staph.html

# install.packages("coda")
library(coda)
mcmc=as.mcmc.resBactDating(res_pr)
effectiveSize(mcmc)
# Known priors (?)
# mu     sigma     alpha 
# 13.86030 136.62327  27.93412

rooted=initRoot(tre,d[,1]) # Incompatible dimensions because of d as matrix of (74,2)
res_roottotip=roottotip(rooted,d[,1])

modell <- c("mixedgamma", "strictgamma")
res_post <- list()

for (i in 1:length(modell)){
  
  res_post[[i]]=bactdate(rooted,d,nbIts=1e6,
                       initMu = effectiveSize(mcmc)["mu"],
                       initAlpha = effectiveSize(mcmc)["alpha"],
                       initSigma = effectiveSize(mcmc)["sigma"],
                       model = modell[i], # try "strictgamma"
                       showProgress = T)
  
  for (j in 1:length(res_post[[i]])){
    
    # MCMC traces
    png(file = paste(modell[i], "_trace.png", sep = ""),
        width = 1200, height = 650)
    plot(res_post[[i]], 'trace')
    dev.off()
    
    # treeCI
    png(file = paste(modell[i], "_treeCI.png", sep = ""),
        width = 1370, height = 960)
    plot(res_post[[i]],'treeCI',show.tip.label = T)
    dev.off()
    
    # treeRoot
    png(file = paste(modell[i], "_treeRoot.png", sep = ""),
        width = 1200, height = 910)
    plot(res_post[[i]],'treeRoot',show.tip.label=T)
    dev.off()
    
    # axisPhylo
    png(file = paste(modell[i], "_axisPhylo.png", sep = ""),
        width = 1300, height = 700)
    out=extractSample(res_post[[i]],6)
    par(mfrow=c(2,3),mar=c(5,5,0.5,0.5))
    for (k in 1:6) {
      plot(out[[k]])
      axisPhylo(1,backward = F)
  }
    dev.off()
  }
}


par(mfrow=c(1,1))
tree=simcoaltree(2001:2015)
plot(tree)
ape::axisPhylo(backward=F)

obsphy=simobsphy(tree)
r=roottotip(obsphy,2001:2015)

res2=bactdate(obsphy,2001:2015)
plot(res2,'treeCI')
