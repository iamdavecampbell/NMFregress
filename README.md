## Installation; from command line:
From the command line:

git clone https://github.com/iamdavecampbell/NMFregress.git



Then source the files by openning R and running:

anchor_files = list.files("../NMFregress", pattern = ".R$")

lapply(paste0("../NMFregress/",anchor_files), source)


## Install Dependencies:

- Used in inference

install.packages(c("mgcv","dplyr","prodlim","betareg"))

- Only used in stratified bootstrap:

install.packages("doParallel")

## Current version:
Version 0.5 release July 15, 2024

Demo in the Romeo and Juliet Rmd file.