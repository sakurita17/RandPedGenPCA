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
set.seed(987765543)
# Simulate founder genomes - one common founder population
founderGenomes <- runMacs(nInd = 100, nChr = 10, segSites = 1100,
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
popA <- founders[1:100]

data <- collectData(pop = popA, data = NULL, population = "A", generation = 0)


# Select on each trait and keep the populations separate
for (generation in 1:20) {
  print(generation)
  popA <- randCross(pop = popA, nCrosses = 100)
  popA <- setPheno(pop = popA, varE = diag(varE))
  data <- collectData(pop = popA, data = data, population = "A", generation = generation)
}


# ---- eIBD COVARIANCE & PRECISION FACTOR --------------------------------

ped <- pedigree(sire = factor(data$pedigree$fid),
                dam = factor(data$pedigree$mid),
                label = factor(data$pedigree$id))

pedLInv <- randPedPCA::sparse2spam(getLInv(ped = ped))

# ---- SAVE DATA ---------------------------------------------------------
pedMeta <- data$pedigree
pedGeno <- data$genoIBS
#usethis::use_data(pedLInv, overwrite = TRUE)
#usethis::use_data(pedMeta, overwrite = TRUE)
