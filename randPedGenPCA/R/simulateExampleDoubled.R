# ---- DESCRIPTION --------------------------------------------------------

# Simulation of three populations, of which two (A and B) are separate and
# the third (AB) is crossbred/admixed with continuous migration (from A and B).
# Populations A and B are selected on two different traits, while population AB
# is selected on an index of these two traits (just to build-up differences).

# The simulation can use one common founder population or two founder
# populations that have diverged some time ago.

# The simulation saves:
# * pedigree information (=expected Identity By Descent (IBD), eIBD)
# * IBD haplotypes since the founders (=realised IBD, rIBD)
# * Identity By State (IBS) genotypes
# * Trait genetic and phenotypic values

# We would only really need eIBD, but we save more for comparison and pedagogy.

# ---- INSTALLATION ------------------------------------------------------

pkgs <- c("AlphaSimR", "dplyr", "pedigreeTools", "SIMplyBee")
install.packages(pkg = pkgs)

# ---- SETUP -------------------------------------------------------------

library(package = "AlphaSimR")
library(package = "dplyr")
library(package = "pedigreeTools")

# ---- SIMULATION - ONE OR TWO FOUNDER POPULATION ------------------------

# Simulate founder genomes - one common founder population
founderGenomes <- runMacs(nInd = 100, nChr = 10, segSites = 1100,
                          species = "GENERIC")

# Simulate founder genomes - two founder populations that have a common origin
founderGenomes2 <- runMacs(nInd = 100, nChr = 10, segSites = 1100, split = 100,
                           species = "GENERIC")
# ... uncomment the line below and run the code to use the two-founder population
# founderGenomes <- founderGenomes2

# Initiate simulation parameters
SP <- SimParam$new(founderPop = founderGenomes)
# ... track global pedigree and recombinations
# (recobinations will enable tracking IBD)
SP$setTrackPed(isTrackPed = TRUE)
SP$setTrackRec(isTrackRec = TRUE)
# ... add two complex traits
varG <- matrix(data = c( 1.0, -0.3,
                        -0.3,  1.0), byrow = TRUE, nrow = 2)
SP$addTraitA(nQtlPerChr = 100, mean = c(0, 0), var = diag(varG), corA = varG)
varE <- matrix(data = c(2.0, 0.0,
                        0.0, 2.0), byrow = TRUE, nrow = 2)

# Monitor function
collectData <- function(pop, data = NULL, population, generation) {
  remove <- FALSE
  if (is.null(data)) {
    remove <- TRUE
    data <- vector(mode = "list", length = 3)
    names(data) <- c("pedigree", "haploIBD", "genoIBS")
    data$pedigree <- data.frame(id = NA, population = NA, generation = NA,
                                mid = NA, fid = NA,
                                gv1 = NA, pv1 = NA,
                                gv2 = NA, pv2 = NA)
    data$haploIBD <- matrix(data = NA, ncol = sum(pop@nLoci))
    data$genoIBS <- matrix(data = NA, ncol = sum(pop@nLoci))
  }
  data$pedigree <- rbind(data$pedigree,
                         data.frame(id = pop@id,
                                    population = population,
                                    generation = generation,
                                    mid = pop@mother,
                                    fid = pop@father,
                                    gv1 = pop@gv[, 1],
                                    pv1 = pop@pheno[, 1],
                                    gv2 = pop@gv[, 2],
                                    pv2 = pop@pheno[, 2]))
  data$haploIBD <- rbind(data$haploIBD,
                         pullIbdHaplo(pop = pop))
  data$genoIBS <- rbind(data$genoIBS,
                        pullSegSiteGeno(pop = pop))
  if (remove) {
    data$pedigree <- data$pedigree[-1, ]
    data$haploIBD <- data$haploIBD[-1, ]
    data$genoIBS <- data$genoIBS[-1, ]
  }
  return(data)
}

# Founder population & split
founders <- newPop(rawPop = founderGenomes)
founders <- setPheno(pop = founders, varE = diag(varE))
popA <- founders[1:50]
popB <- founders[51:100]
data <- collectData(pop = popA, data = NULL, population = "A", generation = 0)
data <- collectData(pop = popB, data = data, population = "B", generation = 0)

# Select on each trait and keep the populations separate
for (generation in 1:10) {
  parentsA <- selectInd(pop = popA, nInd = 10, trait = 1)
  parentsB <- selectInd(pop = popB, nInd = 10, trait = 2)
  popA <- randCross(pop = parentsA, nCrosses = 50)
  popB <- randCross(pop = parentsB, nCrosses = 50)
  popA <- setPheno(pop = popA, varE = diag(varE))
  popB <- setPheno(pop = popB, varE = diag(varE))
  data <- collectData(pop = popA, data = data, population = "A", generation = generation)
  data <- collectData(pop = popB, data = data, population = "B", generation = generation)
}

# Continued selection on each trait in each separate population,
# but add also continually admixed population selected on an index
popAB <- randCross(pop = c(parentsA, parentsB), nCrosses = 50)
popAB <- setPheno(pop = popAB, varE = diag(varE))
data <- collectData(pop = popAB, data = data, population = "AB", generation = generation)
economicWeights <- c(1, 1)
selIndexWeights <- smithHazel(econWt = economicWeights, varG = varG, varP = varG + varE)
for (generation in 11:20) {
  parentsA <- selectInd(pop = popA, nInd = 10, trait = 1)
  parentsB <- selectInd(pop = popB, nInd = 10, trait = 2)
  parentsAB <- selectInd(pop = popAB, nInd = 6, trait = selIndex, scale = TRUE, b = selIndexWeights)
  parentsA4AB <- selectInd(pop = popA, nInd = 2, trait = selIndex, scale = TRUE, b = selIndexWeights)
  parentsB4AB <- selectInd(pop = popB, nInd = 2, trait = selIndex, scale = TRUE, b = selIndexWeights)
  parentsAB <- c(parentsAB, parentsA4AB, parentsB4AB)
  popA <- randCross(pop = parentsA, nCrosses = 50)
  popB <- randCross(pop = parentsB, nCrosses = 50)
  popAB <- randCross(pop = parentsAB, nCrosses = 50)
  popA <- setPheno(pop = popA, varE = diag(varE))
  popB <- setPheno(pop = popB, varE = diag(varE))
  popAB <- setPheno(pop = popAB, varE = diag(varE))
  data <- collectData(pop = popA,  data = data, population = "A",  generation = generation)
  data <- collectData(pop = popB,  data = data, population = "B",  generation = generation)
  data <- collectData(pop = popAB, data = data, population = "AB", generation = generation)
}

data$pedigree$population <- factor(data$pedigree$population, levels = c("A", "B", "AB"))
summary(data$pedigree$population)
data$pedigree$gv <- c(selIndex(Y = as.matrix(data$pedigree[, c("gv1", "gv2")]), scale = TRUE, b = selIndexWeights))
data$pedigree$pv <- c(selIndex(Y = as.matrix(data$pedigree[, c("pv1", "pv2")]), scale = TRUE, b = selIndexWeights))

data$pedigree$generationPlotShift <- data$pedigree$generation +
  c(-0.25, +0.25, 0)[data$pedigree$population]

saveRDS(data, "rhpca_large.rds")

# ---- TRAIT TRENDS ------------------------------------------------------

cols <- adjustcolor(col = c("blue", "red", "violet"), alpha.f = 0.5)
names(cols) <- c("A", "B", "AB")

means <- data$pedigree %>%
  group_by(generationPlotShift, population) %>%
  summarise(across(c(gv1, gv2, gv, pv1, pv2, pv), mean))

par(mfrow = c(3, 1))
plot(x = data$pedigree$generationPlotShift, y = data$pedigree$gv1, pch = 21, cex = 0.5,
     col = cols[data$pedigree$population])
points(x = means$generationPlotShift, y = means$gv1, pch = 19, cex = 2,
       col = cols[means$population])
plot(x = data$pedigree$generationPlotShift, y = data$pedigree$gv2, pch = 21, cex = 0.5,
     col = cols[data$pedigree$population])
points(x = means$generationPlotShift, y = means$gv2, pch = 19, cex = 2,
       col = cols[means$population])
plot(x = data$pedigree$generationPlotShift, y = data$pedigree$gv, pch = 21, cex = 0.5,
     col = cols[data$pedigree$population])
points(x = means$generationPlotShift, y = means$gv, pch = 19, cex = 2,
       col = cols[means$population])

plot(x = data$pedigree$generationPlotShift, y = data$pedigree$pv1, pch = 21, cex = 0.5,
     col = cols[data$pedigree$population])
points(x = means$generationPlotShift, y = means$pv1, pch = 19, cex = 2,
       col = cols[means$population])
plot(x = data$pedigree$generationPlotShift, y = data$pedigree$pv2, pch = 21, cex = 0.5,
     col = cols[data$pedigree$population])
points(x = means$generationPlotShift, y = means$pv2, pch = 19, cex = 2,
       col = cols[means$population])
plot(x = data$pedigree$generationPlotShift, y = data$pedigree$pv, pch = 21, cex = 0.5,
     col = cols[data$pedigree$population])
points(x = means$generationPlotShift, y = means$pv, pch = 19, cex = 2,
       col = cols[means$population])

par(mfrow = c(1, 1))
plot(x = data$pedigree$gv1, y = data$pedigree$gv2, pch = 21, cex = 0.5,
     col = cols[data$pedigree$population])
points(x = means$gv1, y = means$gv2, pch = 19, cex = 2,
       col = cols[means$population])
for (pop in levels(means$population)) {
  tmp <- means[means$population == pop, ]
  lines(x = tmp$gv1, y = tmp$gv2, lwd = 2,
        col = cols[pop])
}

plot(x = data$pedigree$pv1, y = data$pedigree$pv2, pch = 21, cex = 0.5,
     col = cols[data$pedigree$population])
points(x = means$pv1, y = means$pv2, pch = 19, cex = 2,
       col = cols[means$population])
for (pop in levels(means$population)) {
  tmp <- means[means$population == pop, ]
  lines(x = tmp$pv1, y = tmp$pv2, lwd = 2,
        col = cols[pop])
}





# add copy of pedigree as unrelated idivudals -----------------------------

meta <- data.frame(data$pedigree)
head(meta)
plot(meta$generation)
summary(meta$generation)
meta10 <- meta |>
  subset(generation>9)
head(meta10)
ped10 <- pedigree(sire = ifelse(meta10$fid %in% meta10$id, meta10$fid, NA),
                  dam = ifelse(meta10$mid %in% meta10$id, meta10$mid, NA),
                  label = meta10$id)

LInv10 <- getLInv(ped10)
ped10

Matrix::writeMM(obj = LInv10,
                file= "~/git_repos/RandPedigreePCA/datasets/LInv10.mtx")
write.table(meta10,
            file="~/git_repos/RandPedigreePCA/datasets/pedMeta10.tsv",
            sep="\t",
            row.names = F,
            quote = F)
id10 <- meta10$id
s10 <- ifelse(meta10$fid %in% meta10$id, meta10$fid, NA)
d10 <- ifelse(meta10$mid %in% meta10$id, meta10$mid, NA)
tail(id10)


ids10num <- sapply(id10, as.numeric, USE.NAMES = F)
s10num <- sapply(s10, as.numeric, USE.NAMES = F)
d10num <- sapply(d10, as.numeric, USE.NAMES = F)



# Run the whole thing again (messy, tidy up leater!) ----------------------

# Simulate founder genomes - one common founder population
founderGenomes <- runMacs(nInd = 100, nChr = 10, segSites = 1100,
                          species = "GENERIC")

# Simulate founder genomes - two founder populations that have a common origin
founderGenomes2 <- runMacs(nInd = 100, nChr = 10, segSites = 1100, split = 100,
                           species = "GENERIC")
# ... uncomment the line below and run the code to use the two-founder population
# founderGenomes <- founderGenomes2

# Initiate simulation parameters
SP <- SimParam$new(founderPop = founderGenomes)
# ... track global pedigree and recombinations
# (recobinations will enable tracking IBD)
SP$setTrackPed(isTrackPed = TRUE)
SP$setTrackRec(isTrackRec = TRUE)
# ... add two complex traits
varG <- matrix(data = c( 1.0, -0.3,
                         -0.3,  1.0), byrow = TRUE, nrow = 2)
SP$addTraitA(nQtlPerChr = 100, mean = c(0, 0), var = diag(varG), corA = varG)
varE <- matrix(data = c(2.0, 0.0,
                        0.0, 2.0), byrow = TRUE, nrow = 2)

# Monitor function
collectData <- function(pop, data = NULL, population, generation) {
  remove <- FALSE
  if (is.null(data)) {
    remove <- TRUE
    data <- vector(mode = "list", length = 3)
    names(data) <- c("pedigree", "haploIBD", "genoIBS")
    data$pedigree <- data.frame(id = NA, population = NA, generation = NA,
                                mid = NA, fid = NA,
                                gv1 = NA, pv1 = NA,
                                gv2 = NA, pv2 = NA)
    data$haploIBD <- matrix(data = NA, ncol = sum(pop@nLoci))
    data$genoIBS <- matrix(data = NA, ncol = sum(pop@nLoci))
  }
  data$pedigree <- rbind(data$pedigree,
                         data.frame(id = pop@id,
                                    population = population,
                                    generation = generation,
                                    mid = pop@mother,
                                    fid = pop@father,
                                    gv1 = pop@gv[, 1],
                                    pv1 = pop@pheno[, 1],
                                    gv2 = pop@gv[, 2],
                                    pv2 = pop@pheno[, 2]))
  data$haploIBD <- rbind(data$haploIBD,
                         pullIbdHaplo(pop = pop))
  data$genoIBS <- rbind(data$genoIBS,
                        pullSegSiteGeno(pop = pop))
  if (remove) {
    data$pedigree <- data$pedigree[-1, ]
    data$haploIBD <- data$haploIBD[-1, ]
    data$genoIBS <- data$genoIBS[-1, ]
  }
  return(data)
}

# Founder population & split
founders <- newPop(rawPop = founderGenomes)
founders <- setPheno(pop = founders, varE = diag(varE))
popA <- founders[1:50]
popB <- founders[51:100]
data <- collectData(pop = popA, data = NULL, population = "A", generation = 0)
data <- collectData(pop = popB, data = data, population = "B", generation = 0)

# Select on each trait and keep the populations separate
for (generation in 1:10) {
  parentsA <- selectInd(pop = popA, nInd = 10, trait = 1)
  parentsB <- selectInd(pop = popB, nInd = 10, trait = 2)
  popA <- randCross(pop = parentsA, nCrosses = 50)
  popB <- randCross(pop = parentsB, nCrosses = 50)
  popA <- setPheno(pop = popA, varE = diag(varE))
  popB <- setPheno(pop = popB, varE = diag(varE))
  data <- collectData(pop = popA, data = data, population = "A", generation = generation)
  data <- collectData(pop = popB, data = data, population = "B", generation = generation)
}

# Continued selection on each trait in each separate population,
# but add also continually admixed population selected on an index
popAB <- randCross(pop = c(parentsA, parentsB), nCrosses = 50)
popAB <- setPheno(pop = popAB, varE = diag(varE))
data <- collectData(pop = popAB, data = data, population = "AB", generation = generation)
economicWeights <- c(1, 1)
selIndexWeights <- smithHazel(econWt = economicWeights, varG = varG, varP = varG + varE)
for (generation in 11:20) {
  parentsA <- selectInd(pop = popA, nInd = 10, trait = 1)
  parentsB <- selectInd(pop = popB, nInd = 10, trait = 2)
  parentsAB <- selectInd(pop = popAB, nInd = 6, trait = selIndex, scale = TRUE, b = selIndexWeights)
  parentsA4AB <- selectInd(pop = popA, nInd = 2, trait = selIndex, scale = TRUE, b = selIndexWeights)
  parentsB4AB <- selectInd(pop = popB, nInd = 2, trait = selIndex, scale = TRUE, b = selIndexWeights)
  parentsAB <- c(parentsAB, parentsA4AB, parentsB4AB)
  popA <- randCross(pop = parentsA, nCrosses = 50)
  popB <- randCross(pop = parentsB, nCrosses = 50)
  popAB <- randCross(pop = parentsAB, nCrosses = 50)
  popA <- setPheno(pop = popA, varE = diag(varE))
  popB <- setPheno(pop = popB, varE = diag(varE))
  popAB <- setPheno(pop = popAB, varE = diag(varE))
  data <- collectData(pop = popA,  data = data, population = "A",  generation = generation)
  data <- collectData(pop = popB,  data = data, population = "B",  generation = generation)
  data <- collectData(pop = popAB, data = data, population = "AB", generation = generation)
}

data$pedigree$population <- factor(data$pedigree$population, levels = c("A", "B", "AB"))
summary(data$pedigree$population)
data$pedigree$gv <- c(selIndex(Y = as.matrix(data$pedigree[, c("gv1", "gv2")]), scale = TRUE, b = selIndexWeights))
data$pedigree$pv <- c(selIndex(Y = as.matrix(data$pedigree[, c("pv1", "pv2")]), scale = TRUE, b = selIndexWeights))

data$pedigree$generationPlotShift <- data$pedigree$generation +
  c(-0.25, +0.25, 0)[data$pedigree$population]

cols <- adjustcolor(col = c("blue", "red", "violet"), alpha.f = 0.5)
names(cols) <- c("A", "B", "AB")

means <- data$pedigree %>%
  group_by(generationPlotShift, population) %>%
  summarise(across(c(gv1, gv2, gv, pv1, pv2, pv), mean))

par(mfrow = c(3, 1))
plot(x = data$pedigree$generationPlotShift, y = data$pedigree$gv1, pch = 21, cex = 0.5,
     col = cols[data$pedigree$population])
points(x = means$generationPlotShift, y = means$gv1, pch = 19, cex = 2,
       col = cols[means$population])
plot(x = data$pedigree$generationPlotShift, y = data$pedigree$gv2, pch = 21, cex = 0.5,
     col = cols[data$pedigree$population])
points(x = means$generationPlotShift, y = means$gv2, pch = 19, cex = 2,
       col = cols[means$population])
plot(x = data$pedigree$generationPlotShift, y = data$pedigree$gv, pch = 21, cex = 0.5,
     col = cols[data$pedigree$population])
points(x = means$generationPlotShift, y = means$gv, pch = 19, cex = 2,
       col = cols[means$population])

plot(x = data$pedigree$generationPlotShift, y = data$pedigree$pv1, pch = 21, cex = 0.5,
     col = cols[data$pedigree$population])
points(x = means$generationPlotShift, y = means$pv1, pch = 19, cex = 2,
       col = cols[means$population])
plot(x = data$pedigree$generationPlotShift, y = data$pedigree$pv2, pch = 21, cex = 0.5,
     col = cols[data$pedigree$population])
points(x = means$generationPlotShift, y = means$pv2, pch = 19, cex = 2,
       col = cols[means$population])
plot(x = data$pedigree$generationPlotShift, y = data$pedigree$pv, pch = 21, cex = 0.5,
     col = cols[data$pedigree$population])
points(x = means$generationPlotShift, y = means$pv, pch = 19, cex = 2,
       col = cols[means$population])

par(mfrow = c(1, 1))
plot(x = data$pedigree$gv1, y = data$pedigree$gv2, pch = 21, cex = 0.5,
     col = cols[data$pedigree$population])
points(x = means$gv1, y = means$gv2, pch = 19, cex = 2,
       col = cols[means$population])
for (pop in levels(means$population)) {
  tmp <- means[means$population == pop, ]
  lines(x = tmp$gv1, y = tmp$gv2, lwd = 2,
        col = cols[pop])
}

plot(x = data$pedigree$pv1, y = data$pedigree$pv2, pch = 21, cex = 0.5,
     col = cols[data$pedigree$population])
points(x = means$pv1, y = means$pv2, pch = 19, cex = 2,
       col = cols[means$population])
for (pop in levels(means$population)) {
  tmp <- means[means$population == pop, ]
  lines(x = tmp$pv1, y = tmp$pv2, lwd = 2,
        col = cols[pop])
}


# Finished 2nd run --

# Get ped from 2nd run ----------------------------------------------------


meta2nd <- data.frame(data$pedigree)
head(meta2nd)
plot(meta2nd$generation)
summary(meta2nd$generation)
meta2nd10 <- meta2nd |>
  subset(generation>9)
head(meta2nd10)
ped2nd10 <- pedigree(sire = ifelse(meta2nd10$fid %in% meta2nd10$id, meta2nd10$fid, NA),
                  dam = ifelse(meta2nd10$mid %in% meta2nd10$id, meta2nd10$mid, NA),
                  label = meta2nd10$id)

id2nd10 <- meta2nd10$id
s2nd10 <- ifelse(meta2nd10$fid %in% meta2nd10$id, meta2nd10$fid, NA)
d2nd10 <- ifelse(meta2nd10$mid %in% meta2nd10$id, meta2nd10$mid, NA)
tail(id2nd10)


ids2nd10num <- sapply(id2nd10, as.numeric, USE.NAMES = F)
s2nd10num <- sapply(s2nd10, as.numeric, USE.NAMES = F)
d2nd10num <- sapply(d2nd10, as.numeric, USE.NAMES = F)






# doubled up version
maxId <- max(ids10num) + 10
ids1010 <- c(ids10num, ids2nd10num + maxId)
s1010 <- c(s10num, s2nd10num + maxId)
d1010 <- c(d10num, d2nd10num + maxId)
ped1010 <- pedigree(
  s1010,
  d1010,
  ids1010
)
LInv1010 <- getLInv(ped1010)
Matrix::writeMM(obj = LInv1010,
                file= "~/git_repos/RandPedigreePCA/datasets/LInv1010.mtx")

