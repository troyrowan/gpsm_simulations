library(AlphaSimR)
#library(tidyverse)
library(tidyr)
library(dplyr)
library(ggplot2)
library(readr)
library(purrr)
library(cowplot)
library(reshape2)
library(qvalue)
library(stringr)

set.seed(4847558)

args <- commandArgs(TRUE)

#This reads in all of the arguments from the param file
source(args[1])
#Uses already-created founder population (as this takes quite some time)
#Commands below:
# founderPop = runMacs(nInd=2000,
#                      nChr=10,
#                      segSites=5000,
#                      species = "CATTLE")
#founderPop = readRDS("founderpop.RDS")
#This reads in same founder population that we use to do gene drop in Red Angus Pedigree 200K genotypes total for 7601 founder individuals
#Not sure how this'll work as we do selection (particularly WRT burn-in, but we'll see I guess)
founderPop = readRDS(founderpop)
#founderPop = readRDS("100K_cattle_founderpop.RDS")
SP = SimParam$new(founderPop)
SP$addTraitA(nQtlPerChr=nqtl, gamma = gamma, shape = shape)
SP$setGender("yes_sys")
SP$addSnpChip(nSnpPerChr = nsnps)
SP$setVarE(h2 = h2)

pulln = case_when(gens == 5 ~ list(even_5gen, uneven_5gen),
									gens == 10 ~ list(even_10gen, uneven_10gen),
									gens == 20 ~ list(even_20gen, uneven_20gen))

sampler_even = pulln[[1]]
sampler_uneven = pulln[[2]]

# True QTL effects
trait = data.frame(chr = rep(1:chrom, each = nqtl),
                   pos = SP$traits[[1]]@lociLoc,
                   effect = SP$traits[[1]]@addEff,
                   rs = paste(rep(1:chrom, each = nqtl), SP$traits[[1]]@lociLoc, sep = ":"))
write_csv(trait, paste0("generation_genotypes/", testname, "/rep", rep, "/", testname, ".rep", rep,".true_qtl.csv"))

#Map file creation
data.frame(rs = paste0("SNP_", 1:(nsnps*chrom)),
                      pos = SP$snpChips[[1]]@lociLoc,
                      chr = rep(1:chrom, each = nsnps),
                      cm = "0") %>%
                      select(chr, rs, cm, pos) %>%
                      write_tsv(paste0("generation_genotypes/", testname, "/rep", rep, "/", testname,".rep", rep, ".map"),
                                col_names = FALSE)
#Creates file for tracking genetic values and variances in data over generations
pops = data.frame(gen = as.numeric(),
                  pheno_g = as.numeric(),
                  pheno_vg = as.numeric(),
                  gv_g = as.numeric(),
                  gv_vg = as.numeric(),
                  rand_g = as.numeric(),
                  rand_vg = as.numeric())
#Creating empty dataframes for genotypes (even sampling)
gv_geno_even = data.frame(matrix(ncol = chrom*nsnps, nrow = 0))
colnames(gv_geno_even)<-paste("SNP", seq(1,nsnps*chrom), sep = "_")
rand_geno_even = data.frame(matrix(ncol = chrom*nsnps, nrow = 0))
colnames(rand_geno_even)<-paste("SNP", seq(1,nsnps*chrom), sep = "_")
pheno_geno_even = data.frame(matrix(ncol = chrom*nsnps, nrow = 0))
colnames(pheno_geno_even)<-paste("SNP", seq(1,nsnps*chrom), sep = "_")

#Creating empty dataframes for genotypes (uneven sampling)
gv_geno_uneven = data.frame(matrix(ncol = chrom*nsnps, nrow = 0))
colnames(gv_geno_uneven)<-paste("SNP", seq(1,nsnps*chrom), sep = "_")
rand_geno_uneven = data.frame(matrix(ncol = chrom*nsnps, nrow = 0))
colnames(rand_geno_uneven)<-paste("SNP", seq(1,nsnps*chrom), sep = "_")
pheno_geno_uneven = data.frame(matrix(ncol = chrom*nsnps, nrow = 0))
colnames(pheno_geno_uneven)<-paste("SNP", seq(1,nsnps*chrom), sep = "_")

#Duplicates founder population into the three different selection populations that AlphaSimR will use
pheno_pop = newPop(founderPop)
gv_pop = newPop(founderPop)
rand_pop = newPop(founderPop)
#Creates files for tracking true QTL allele frequency changes over time
gv_qtl = bind_cols(trait, data.frame(matrix(nrow=nqtl*chrom, ncol=gens-burnins)))
pheno_qtl =  bind_cols(trait, data.frame(matrix(nrow=nqtl*chrom, ncol=gens-burnins)))
rand_qtl =  bind_cols(trait, data.frame(matrix(nrow=nqtl*chrom, ncol=gens-burnins)))
#Loops through burn-ins and selection generations
for(generation in -burnins:gens){
#Adding in the generation number in Vg G tracking file
  pops[generation + burnins + 1 , "gen"] = generation

  if (generation == -burnins){
    pops[generation + burnins + 1 , "rand_g"] = meanG(rand_pop)
    pops[generation + burnins + 1 , "rand_vg"] = varG(rand_pop)
    pops[generation + burnins + 1 , "pheno_g"] = meanG(pheno_pop)
    pops[generation + burnins + 1 , "pheno_vg"] = varG(pheno_pop)
		pops[generation + burnins + 1 , "gv_g"] = meanG(gv_pop)
    pops[generation + burnins + 1 , "gv_vg"] = varG(gv_pop)


  }

  if (generation < 0){
		rand_pop = selectCross(pop=rand_pop,
                         nFemale=females,
                         nMale=males,
                         use="pheno",
                         nCrosses=burnincrosses)
		pheno_pop = selectCross(pop=pheno_pop,
                         nFemale=females,
                         nMale=males,
                         use="pheno",
                         nCrosses=burnincrosses)
		gv_pop = selectCross(pop=gv_pop,
                         nFemale=females,
                         nMale=males,
                         use="pheno",
                         nCrosses=burnincrosses)
    pops[generation + burnins + 1 , "rand_g"] = meanG(rand_pop) #Updates genetic values in dataframe
    pops[generation + burnins + 1 , "rand_vg"] = varG(rand_pop) #Updates genetic values in dataframe
		pops[generation + burnins + 1 , "pheno_g"] = meanG(pheno_pop) #Updates genetic values in dataframe
    pops[generation + burnins + 1 , "phenho_vg"] = varG(pheno_pop) #Updates genetic values in dataframe
		pops[generation + burnins + 1 , "gv_g"] = meanG(gv_pop) #Updates genetic values in dataframe
    pops[generation + burnins + 1 , "gv_vg"] = varG(gv_pop) #Updates genetic values in dataframe

				print("burnins completed")
  }

  if (generation >= 0){
    rand_pop = selectCross(pop=rand_pop, #Sets genomic EBV via RRBLUP for selection
                         nFemale=females,
                         nMale=males,
                         use="rand",
                         nCrosses=crosses)
    pops[generation + burnins + 1 , "rand_g"] = meanG(rand_pop)
    pops[generation + burnins + 1 , "rand_vg"] = varG(rand_pop)

    rand_geno_even = bind_rows(rand_geno_even,
                             as.data.frame(pullSnpGeno(rand_pop)) %>%
                               sample_n(sampler_even[generation+1,]$keepind))
    #For uneven sampling of genotypes (but from the same selected population)
    rand_geno_uneven = bind_rows(rand_geno_uneven,
                               as.data.frame(pullSnpGeno(rand_pop)) %>%
                                 sample_n(sampler_uneven[generation+1,]$keepind))

      rand_qtl[,(generation+1+ncol(trait))] <- pullQtlGeno(rand_pop) %>%
        as.data.frame() %>%
        colSums()/(2*crosses) %>%
        #t() %>%
        as.vector()
				print(paste0("Random population generation", generation))


    pheno_pop = selectCross(pop=pheno_pop, #Sets genomic EBV via RRBLUP for selection
                         nFemale=females,
                         nMale=males,
                         use="pheno",
                         nCrosses=crosses)
    pops[generation + burnins + 1 , "pheno_g"] = meanG(pheno_pop)
    pops[generation + burnins + 1 , "pheno_vg"] = varG(pheno_pop)

    pheno_geno_even = bind_rows(pheno_geno_even,
                             as.data.frame(pullSnpGeno(pheno_pop)) %>%
                               sample_n(sampler_even[generation+1,]$keepind))
    #For uneven sampling of genotypes (but from the same selected population)
    pheno_geno_uneven = bind_rows(pheno_geno_uneven,
                               as.data.frame(pullSnpGeno(pheno_pop)) %>%
                                 sample_n(sampler_uneven[generation+1,]$keepind))

      rand_qtl[,(generation+1+ncol(trait))] <- pullQtlGeno(pheno_pop) %>%
        as.data.frame() %>%
        colSums()/(2*crosses) %>%
        #t() %>%
        as.vector()
				print(paste0("Phenotypic selection population generation", generation))

    gv_pop = selectCross(pop=gv_pop, #Sets genomic EBV via RRBLUP for selection
                         nFemale=females,
                         nMale=males,
                         use="gv",
                         nCrosses=crosses)
    pops[generation + burnins + 1 , "gv_g"] = meanG(gv_pop)
    pops[generation + burnins + 1 , "gv_vg"] = varG(gv_pop)

    gv_geno_even = bind_rows(gv_geno_even,
                             as.data.frame(pullSnpGeno(gv_pop)) %>%
                               sample_n(sampler_even[generation+1,]$keepind))
    #For uneven sampling of genotypes (but from the same selected population)
    gv_geno_uneven = bind_rows(gv_geno_uneven,
                               as.data.frame(pullSnpGeno(gv_pop)) %>%
                                 sample_n(sampler_uneven[generation+1,]$keepind))

      gv_qtl[,(generation+1+ncol(trait))] <- pullQtlGeno(gv_pop) %>%
        as.data.frame() %>%
        colSums()/(2*crosses) %>%
        #t() %>%
        as.vector()
				print(paste0("TBV selection population generation", generation))
  }
}

#Recoding to allow for PLINK file usage
rand_geno_even =
  rand_geno_even %>%
  mutate_all(funs(case_when(. == 0 ~ "AA",
                            . == 1 ~ "AB",
                            . == 2 ~ "BB")))

rand_geno_uneven =
  rand_geno_uneven %>%
  mutate_all(funs(case_when(. == 0 ~ "AA",
                            . == 1 ~ "AB",
                            . == 2 ~ "BB")))
print("Random selection population generation GT coding changed")

pheno_geno_even =
  pheno_geno_even %>%
  mutate_all(funs(case_when(. == 0 ~ "AA",
                            . == 1 ~ "AB",
                            . == 2 ~ "BB")))
pheno_geno_uneven =
  pheno_geno_uneven %>%
  mutate_all(funs(case_when(. == 0 ~ "AA",
                            . == 1 ~ "AB",
                            . == 2 ~ "BB")))
print("Phenotypic selection population generation GT coding changed")

gv_geno_even =
  gv_geno_even %>%
  mutate_all(funs(case_when(. == 0 ~ "AA",
                            . == 1 ~ "AB",
                            . == 2 ~ "BB")))
gv_geno_uneven =
  gv_geno_uneven %>%
  mutate_all(funs(case_when(. == 0 ~ "AA",
                            . == 1 ~ "AB",
                            . == 2 ~ "BB")))
print("TBV selection population generation GT coding changed")

if ("rand" %in% analyses){
  bind_cols(
  data.frame(famid =  1,
             id = seq(1:nrow(rand_geno_even)),
             sire = 0,
             dam = 0,
             sex = 0,
             generation = rep(sampler_even$gen, sampler_even$keepind)),
    rand_geno_even) %>%
    write_delim(paste0("generation_genotypes/", testname, "/rep", rep, "/", testname, ".rand_pop.even.rep", rep,".ped"),
              delim = " ",
              col_names = FALSE)


  bind_cols(
  data.frame(famid =  1,
             id = seq(1:nrow(rand_geno_uneven)),
             sire = 0,
             dam = 0,
             sex = 0,
             generation = rep(sampler_uneven$gen, sampler_uneven$keepind)),
    rand_geno_uneven) %>%
    write_delim(paste0("generation_genotypes/", testname, "/rep", rep, "/", testname, ".rand_pop.uneven.rep", rep,".ped"),
              delim = " ",
              col_names = FALSE)

  rand_qtl %>%
    write_csv(paste0("gpsm_runs/", testname, "/rep", rep, "/", testname, ".rand_pop.rep", rep, ".qtl_trajectories.csv"))
}

print("Random selection genotypes written")

if ("pheno" %in% analyses){
  bind_cols(
  data.frame(famid =  1,
             id = seq(1:nrow(pheno_geno_even)),
             sire = 0,
             dam = 0,
             sex = 0,
             generation = rep(sampler_even$gen, sampler_even$keepind)),
    pheno_geno_even) %>%
    write_delim(paste0("generation_genotypes/", testname, "/rep", rep, "/", testname, ".pheno_pop.even.rep", rep,".ped"),
                delim = " ",
                col_names = FALSE)
  bind_cols(
  data.frame(famid =  1,
             id = seq(1:nrow(pheno_geno_uneven)),
             sire = 0,
             dam = 0,
             sex = 0,
             generation = rep(sampler_uneven$gen, sampler_uneven$keepind)),
    pheno_geno_uneven) %>%
    write_delim(paste0("generation_genotypes/", testname, "/rep", rep, "/", testname, ".pheno_pop.uneven.rep", rep,".ped"),
                delim = " ",
                col_names = FALSE)

  pheno_qtl %>%
    write_csv(paste0("gpsm_runs/", testname, "/rep", rep, "/", testname, ".pheno_pop.rep", rep, ".qtl_trajectories.csv"))
}

print("Phenotypic selection genotypes written")

if ("gv" %in% analyses){
  bind_cols(
  data.frame(famid =  1,
             id = seq(1:nrow(gv_geno_even)),
             sire = 0,
             dam = 0,
             sex = 0,
             generation = rep(sampler_even$gen, sampler_even$keepind)),
    gv_geno_even) %>%
    write_delim(paste0("generation_genotypes/", testname, "/rep", rep, "/", testname, ".gv_pop.even.rep", rep,".ped"),
                delim = " ",
                col_names = FALSE)
  bind_cols(
  data.frame(famid =  1,
             id = seq(1:nrow(gv_geno_uneven)),
             sire = 0,
             dam = 0,
             sex = 0,
             generation = rep(sampler_uneven$gen, sampler_uneven$keepind)),
    gv_geno_uneven) %>%
    write_delim(paste0("generation_genotypes/", testname, "/rep", rep, "/", testname, ".gv_pop.uneven.rep", rep,".ped"),
                delim = " ",
                col_names = FALSE)

  gv_qtl %>%
    write_csv(paste0("gpsm_runs/", testname, "/rep", rep, "/", testname, ".gv_pop.rep", rep, ".qtl_trajectories.csv"))
}

print("Random selection genotypes written")

plot_grid(
  pops %>%
    select(gen, pheno_g, gv_g, rand_g) %>%
    melt(id = c("gen")) %>%
    ggplot(aes(gen, value, color = variable))+
    geom_point(alpha = 0.3)+
    geom_smooth(size = 1.5)+
    labs(y = "Mean G", x = "Generation", main = "Genetic Gain") +
    scale_color_viridis_d(end = 0.75)+
  cowplot::theme_cowplot(),
  #facet_grid(~burnin, scales = "free_x")

  pops %>%
    select(gen, pheno_vg, gv_vg, rand_vg) %>%
    melt(id = c("gen")) %>%
    ggplot(aes(gen, value, color = variable))+
    geom_point(alpha = 0.3)+
    geom_smooth(size = 1.5)+
    labs(y = "Var G", x = "Generation", main = "Genetic Variance")+
    scale_color_viridis_d(end = 0.75)+
  cowplot::theme_cowplot(),
  ncol = 1)
ggsave(paste0("gpsm_runs/", testname, "/rep", rep, "/figures/", testname,".rep", rep, ".g_vg.trends.png"), width = 8, height = 12)
