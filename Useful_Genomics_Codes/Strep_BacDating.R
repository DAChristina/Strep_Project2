
if (!requireNamespace("BactDating", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
  }
  library(devtools)
  devtools::install_github("xavierdidelot/BactDating")
  
  
  library(BactDating)
}
