# Mammal-dung beetle association network
# August 2025
# Juan Vargas

library(tidyverse)
library(bipartite)
library(vegan)

### Import data ###
# Import camera-trap data as matrix for vegan
comm <- read.csv("Data/vert_comm_mat.csv", row.names = 1)

# Import pitfall trap data
commb <- read.csv("Data/beet_comm_mat.csv", row.names = 1)  

# Find sites in common, subset data frames
commonsites <- intersect(rownames(comm), rownames(commb))
commb_commonsites <- commb[commonsites,]
comm_commonsites <- comm[commonsites,]

#### JSDM-derived co-occurrence analysis ####

##Here we follow Peres-Neto et al., which quantify and estimate the significance
#of co-occurrence metrics among pairs of species. I will use the posterior for
#the presence variable from the JSDM models to do co-occurrence permutations.
#The model already contains an array with z presence/absence values that can be
#used for the permutation. I need to subset only the sites that were included in
#both models.

# import the best models
vert_jsdm <- readRDS("vert_jsdm.rds")
beet_jsdm <- readRDS("beet_jsdm.rds")

# extract site names and find those in common
vert_jsdm_sites <- dimnames(vert_jsdm$y)[["sites"]]
beet_jsdm_sites <- dimnames(beet_jsdm$y)[["sites"]]

common_jsdm_sites <- beet_jsdm_sites[beet_jsdm_sites %in% vert_jsdm_sites]
# There are 36 sites included in both models
comm_commonsites <- comm[common_jsdm_sites, ]
commb_commonsites <- commb[common_jsdm_sites, ]
# remove species that did not appear at those sites
comm_commonsites <- comm_commonsites[, colSums(comm_commonsites)>0]
commb_commonsites <- commb_commonsites[, colSums(commb_commonsites)>0]
# this leaves 60 beetles species and 24 vert species

dim(vert_jsdm$z.samples)
# dims are 1200 samples x 28 species x 86 sites. We select only the sites in
# common, and only the species in common
vertzsamples <- vert_jsdm$z.samples[, match(names(comm_commonsites),
                                            vert_jsdm$sp.names),
                                    match(common_jsdm_sites,
                                          dimnames(vert_jsdm$y)$sites)]

dim(beet_jsdm$z.samples)
# for beetles the dims are 200 x 67 x 60, we do the same filtering
beetzsamples <- beet_jsdm$z.samples[, match(names(commb_commonsites),
                                                beet_jsdm$sp.names),
                                        match(common_jsdm_sites,
                                              dimnames(beet_jsdm$y)$sites)]

# the resulting arrays have info on 36 sites, for 24 spp of verts and 60 of beetles.

# This gives us occurrence matrices with which we can calculate a p value of
# co-occurrence for all mammal-beetle pairs
nverts <- dim(vertzsamples)[2]
nbeets <- dim(beetzsamples)[2]
coocc_permuts <- array(dim = c(nverts, nbeets, 200))
for (mm in 1:nverts) {
  for (bb in 1:nbeets) {
    # create 200x36 matrix of presences/absences for species mm and bb at all 36
    # sites (200 simulations)
    A <- vertzsamples[1:200, mm, ]
    B <- beetzsamples[1:200, bb, ]
    # calculate joint presences and joint absences, multiply together 
    n_jointpres <- rowSums(A + B == 2, na.rm = T)
    n_jointabs <- rowSums(A + B == 0, na.rm = T)
    Tij <- n_jointabs * n_jointpres
    # Cs <- rowSums(A+B==1, na.rm = T)/nvalidsites
    # Ts <- rowSums(A==B, na.rm = T)/nvalidsites
    # Ss <- rowSums(A+B==2, na.rm = T)/nvalidsites
    # coocc_permuts[mm, bb,1,] <- Cs
    coocc_permuts[mm, bb,] <- Tij
    # coocc_permuts[mm, bb,3,] <- Ss
  }
}
dim(coocc_permuts)


# Base co-occurrence observed togetherness Ts
coocc_indices <- array(dim = c(ncol(comm_commonsites), ncol(commb_commonsites)))
for (i in 1 : ncol(comm_commonsites)) {
  for(j in 1 : ncol(commb_commonsites)) {
    # create the presence matrices
    A <- comm_commonsites[, i] > 0
    B <- commb_commonsites[, j] > 0
    # AB <- cbind(A,B)
    # Cs <- sum(rowSums(AB)==1)/70
    n_jointpres <- sum(A + B == 2)
    n_jointabs <- sum(A + B == 0)
    Tij <- n_jointabs * n_jointpres
    # Ss <- sum(rowSums(AB)==2)/70
    coocc_indices[i, j] <- Tij
  }
}
coocc_indices
# These are the naive co-occurrence indices, but now I want to know which of
# these associations are statistically significant. We calculate a p-value by
# adding the index to the permutation values and calculating the proportion
# greater than observed.
perm_pvals <- matrix(nrow = nverts, ncol = nbeets)
for (mm in 1:nverts) {
  for (bb in 1:nbeets) {
    obsval <- coocc_indices[mm, bb]
    permdist <- c(coocc_permuts[mm, bb, ], obsval)
    perm_pvals[mm, bb] <- sum(permdist > obsval) / length(permdist)
  }
}
perm_pvals
sum(perm_pvals < 0.05)
# According to the models only 41 out of the 1440 potential pairs have
# significant pairwise togetherness statistic


#### Visualization ####
vertnames <- gsub(pattern = "_",
                  replacement = " ",
                  names(comm_commonsites)) %>% str_to_title()
beetnames <- gsub(pattern = "_",
                  replacement = " ",
                  names(commb_commonsites)) 
rownames(perm_pvals) <- vertnames
colnames(perm_pvals) <- beetnames

# bipartite network
plotweb(t(perm_pvals < 0.05),
        sorting = "dec",
        horizontal = TRUE,
        lower_italic = TRUE,
        text_size = 0.5,
        scaling = "relative", 
        lower_color = "dodgerblue",
        higher_color = "goldenrod",
        link_color = "gray"
        )
apply(perm_pvals < 0.05, 2, mean) |> sort()
apply(perm_pvals < 0.05, 1, mean) |> sort()
