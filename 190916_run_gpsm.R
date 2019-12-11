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

#This sources a config file with information about simulation parameters
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
founderPop = readRDS("RAN_founderHaps.50K.7601.rds")
#founderPop = readRDS("100K_cattle_founderpop.RDS")
SP = SimParam$new(founderPop)
SP$addTraitA(nQtlPerChr=nqtl, gamma = gamma, shape = shape)
SP$setGender("yes_sys")
SP$addSnpChip(nSnpPerChr = nsnps)
SP$setVarE(h2 = h2)
#SP$nThreads = 1L

# even_20gen = data.frame(gen = 20:1, keepind = 500)
# even_10gen = data.frame(gen = 10:1, keepind = 1000)
# even_5gen = data.frame(gen = 5:1, keepind = 2000)
# uneven_20gen = data.frame(gen = 20:1, keepind = c(1478, 2408, 1670, 1200, 881, 652, 450, 350, 253, 160, 146, 105, 77, 43, 39, 24, 20, 16, 16, 12))
# uneven_10gen = data.frame(gen = 10:1, keepind = c(3126, 3604, 1696, 819, 423, 190, 76, 27, 25, 14))
# uneven_5gen = data.frame(gen = 5:1, keepind = c(7329, 1919, 560, 144, 48))

#Reads these dataframes from the param file. Use the number of generations to select both the even and uneven samplign schemes
#By doing this we can do two different genotype sampling strategies on a single population
pulln = case_when(gens == 5 ~ list(even_5gen, uneven_5gen),
									gens == 10 ~ list(even_10gen, uneven_10gen),
									gens == 20 ~ list(even_20gen, uneven_20gen))

sampler_even = pulln[[1]]
sampler_uneven = pulln[[2]]

#gens =  deparse(substitute(sampler)) %>% str_split_fixed(., "_", 2)[2] %>% str_replace("gen", "") %>% as.numeric()
trait = data.frame(chr = rep(1:chrom, each = nqtl),
                   pos = SP$traits[[1]]@lociLoc,
                   effect = SP$traits[[1]]@addEff,
                   rs = paste(rep(1:chrom, each = nqtl), SP$traits[[1]]@lociLoc, sep = ":"))
write_csv(trait, paste0("generation_genotypes/", testname, "/", testname, ".true_qtl.csv"))

annotation_file = data.frame(rs = paste0("SNP_", 1:(nsnps*chrom)),
                             pos = SP$snpChips[[1]]@lociLoc,
                             chr = rep(1:chrom, each = nsnps))

write_delim(annotation_file,
            paste0("generation_genotypes/", testname, "/", testname, ".snp.annotation.txt"),
            delim = ",",
            col_names = FALSE)

pheno_pop = gv_pop = rand_pop = newPop(founderPop) #initializes the same base population allows us to select these in parallel

pops = data.frame(gen = as.numeric(),
                  pheno_g = as.numeric(),
                  pheno_vg = as.numeric(),
                  gv_g = as.numeric(),
                  gv_vg = as.numeric(),
                  rand_g = as.numeric(),
                  rand_vg = as.numeric())

gv_geno_even = data.frame(matrix(ncol = chrom*nsnps, nrow = 0))
colnames(gv_geno_even)<-paste("SNP", seq(1,nsnps*chrom), sep = "_")
rand_geno_even = data.frame(matrix(ncol = chrom*nsnps, nrow = 0))
colnames(rand_geno_even)<-paste("SNP", seq(1,nsnps*chrom), sep = "_")
pheno_geno_even = data.frame(matrix(ncol = chrom*nsnps, nrow = 0))
colnames(pheno_geno_even)<-paste("SNP", seq(1,nsnps*chrom), sep = "_")

gv_geno_uneven = data.frame(matrix(ncol = chrom*nsnps, nrow = 0))
colnames(gv_geno_uneven)<-paste("SNP", seq(1,nsnps*chrom), sep = "_")
rand_geno_uneven = data.frame(matrix(ncol = chrom*nsnps, nrow = 0))
colnames(rand_geno_uneven)<-paste("SNP", seq(1,nsnps*chrom), sep = "_")
pheno_geno_uneven = data.frame(matrix(ncol = chrom*nsnps, nrow = 0))
colnames(pheno_geno_uneven)<-paste("SNP", seq(1,nsnps*chrom), sep = "_")

pheno_pop = newPop(founderPop)
gv_pop = newPop(founderPop)
rand_pop = newPop(founderPop)

gv_qtl = bind_cols(trait, data.frame(matrix(nrow=nqtl*chrom, ncol=gens-burnins)))
pheno_qtl =  bind_cols(trait, data.frame(matrix(nrow=nqtl*chrom, ncol=gens-burnins)))
rand_qtl =  bind_cols(trait, data.frame(matrix(nrow=nqtl*chrom, ncol=gens-burnins)))

for(generation in 0:gens){
  pops[generation + 1, "gen"] = generation

  if (generation == 0){
    pops[generation + 1, "pheno_g"] = meanG(newPop(founderPop))
    pops[generation + 1, "pheno_vg"] = varG(newPop(founderPop))
    pops[generation + 1, "gv_g"] = meanG(newPop(founderPop))
    pops[generation + 1, "gv_vg"] = varG(newPop(founderPop))
    pops[generation + 1, "rand_g"] = meanG(newPop(founderPop))
    pops[generation + 1, "rand_vg"] = varG(newPop(founderPop))
  }

  if (generation > 0 & generation <= burnins){
    pheno_pop = selectCross(pop=pheno_pop,
                            nFemale=females,
                            nMale=males,
                            use="pheno",
                            nCrosses=burnincrosses) #specify selection intensity, selection criteria, number of offspring
    pops[generation + 1, "pheno_g"] = meanG(pheno_pop) #Updates genetic values in dataframe
    pops[generation + 1, "pheno_vg"] = varG(pheno_pop) #Updates genetic values in dataframe

    gv_pop = selectCross(pop=gv_pop,
                         nFemale=females,
                         nMale=males,
                         use="pheno",
                         nCrosses=burnincrosses)
    pops[generation + 1, "gv_g"] = meanG(gv_pop) #Updates genetic values in dataframe
    pops[generation + 1, "gv_vg"] = varG(gv_pop) #Updates genetic values in dataframe

    rand_pop = selectCross(pop=rand_pop,
                           nFemale=females,
                           nMale=males,
                           use="pheno",
                           nCrosses=burnincrosses)
    pops[generation + 1, "rand_g"] = meanG(rand_pop) #Updates genetic values in dataframe
    pops[generation + 1, "rand_vg"] = varG(rand_pop) #Updates genetic values in dataframe
  }

  if (generation > burnins){
    pheno_pop = selectCross(pop=pheno_pop,
                            nFemale=females,
                            nMale=males,
                            use="pheno",
                            nCrosses=crosses) #specify selection intensity, selection criteria, number of offspring
    pops[generation + 1, "pheno_g"] = meanG(pheno_pop) #Updates genetic values in dataframe
    pops[generation + 1, "pheno_vg"] = varG(pheno_pop) #Updates genetic values in dataframe

    gv_pop = selectCross(pop=gv_pop, #Sets genomic EBV via RRBLUP for selection
                         nFemale=females,
                         nMale=males,
                         use="gv",
                         nCrosses=crosses)
    pops[generation + 1, "gv_g"] = meanG(gv_pop)
    pops[generation + 1, "gv_vg"] = varG(gv_pop)

    rand_pop = selectCross(pop=rand_pop,
                           nFemale=females,
                           nMale=males,
                           use="rand",
                           nCrosses=crosses)
    pops[generation + 1, "rand_g"] = meanG(rand_pop)
    pops[generation + 1, "rand_vg"] = varG(rand_pop)

    gv_geno_even = bind_rows(gv_geno_even,
                        as.data.frame(pullSnpGeno(gv_pop)) %>%
                          sample_n(sampler_even[generation-burnins,]$keepind))

    pheno_geno_even = bind_rows(pheno_geno_even,
                           as.data.frame(pullSnpGeno(pheno_pop)) %>%
                             sample_n(sampler_even[generation-burnins,]$keepind))

    rand_geno_even = bind_rows(rand_geno_even,
                          as.data.frame(pullSnpGeno(rand_pop)) %>%
                            sample_n(sampler_even[generation-burnins,]$keepind))

#For uneven sampling of genotypes (but from the same selected population)
    gv_geno_uneven = bind_rows(gv_geno_uneven,
                    as.data.frame(pullSnpGeno(gv_pop)) %>%
                      sample_n(sampler_uneven[generation-burnins,]$keepind))

    pheno_geno_uneven = bind_rows(pheno_geno_uneven,
                       as.data.frame(pullSnpGeno(pheno_pop)) %>%
                         sample_n(sampler_uneven[generation-burnins,]$keepind))

    rand_geno_uneven = bind_rows(rand_geno_uneven,
                      as.data.frame(pullSnpGeno(rand_pop)) %>%
                        sample_n(sampler_uneven[generation-burnins,]$keepind))

gv_qtl[,(generation-burnins+ncol(trait))] <- pullQtlGeno(gv_pop) %>% as.data.frame() %>% colSums()/(2*crosses) %>% t() %>% as.vector()
pheno_qtl[,(generation-burnins+ncol(trait))] <- pullQtlGeno(pheno_pop) %>% as.data.frame() %>% colSums()/(2*crosses) %>% t() %>% as.vector()
rand_qtl[,(generation-burnins+ncol(trait))] <- pullQtlGeno(rand_pop) %>% as.data.frame() %>% colSums()/(2*crosses) %>% t() %>% as.vector()

    # list(gv_pop, rand_pop) %>%
    #   set_names(c("gv_pop", "rand_pop")) %>%
    # iwalk(
    # ~pullSnpGeno(.x) %>%
    #   as.data.frame() %>%
    #   filter(row_number() %% pulleach == 0) %>%
    #   t() %>%
    #   as.data.frame() %>%
    #   mutate(rs = paste(rep(1:10, each = nqtl), seq(1:nqtl), sep = ":"),
    #          a1 = "a1",
    #          a2 = "a2") %>%
    #   select(rs, a1, a2, everything()) %>%
    #   write_csv(paste0("generation_genotypes/", testname, ".", .y, ".replicate", run+1, ".genotypes.mgf"), append = TRUE))

    # list(gv_pop, pheno_pop, rand_pop) %>%
    #   set_names(c("gv_pop", "pheno_pop", "rand_pop")) %>%
    #   iwalk(
    #     ~pheno(.x) %>%
    #       as.data.frame() %>%
    #       filter(row_number() %% pulleach == 0) %>%
    #       write_tsv(paste0("generation_genotypes/", testname, "/", testname, ".", .y, ".trait_phenotypes.txt"), append = TRUE))
  }
}

#These will first check to see if we're doing GV/Pheno/Random selection, then write genotypes for both the even and uneven sampling strategies

if ("gv" %in% analyses){
  gv_geno_even %>%
    t() %>%
    as.data.frame() %>%
    mutate(rs = rownames(.),
           a1 = "a1",
           a2 = "a2") %>%
    select(rs, a1, a2, everything()) %>%
    write_csv(paste0("generation_genotypes/", testname, "/", testname, ".gv_pop.even.genotypes.mgf"),
              col_names = FALSE)
  gv_geno_uneven %>%
    t() %>%
    as.data.frame() %>%
    mutate(rs = rownames(.),
           a1 = "a1",
           a2 = "a2") %>%
    select(rs, a1, a2, everything()) %>%
    write_csv(paste0("generation_genotypes/", testname, "/", testname, ".gv_pop.uneven.genotypes.mgf"),
            col_names = FALSE)
  gv_qtl %>%
	write_csv(paste0("gpsm_runs/", testname, "/", testname, ".gv_pop.qtl_trajectories.csv"))
}

if ("pheno" %in% analyses){
  pheno_geno_even %>%
    t() %>%
    as.data.frame() %>%
    mutate(rs = rownames(.),
           a1 = "a1",
           a2 = "a2") %>%
    select(rs, a1, a2, everything()) %>%
    write_csv(paste0("generation_genotypes/", testname, "/", testname, ".pheno_pop.even.genotypes.mgf"),
              col_names = FALSE)

  pheno_geno_uneven %>%
    t() %>%
    as.data.frame() %>%
    mutate(rs = rownames(.),
           a1 = "a1",
           a2 = "a2") %>%
    select(rs, a1, a2, everything()) %>%
    write_csv(paste0("generation_genotypes/", testname, "/", testname, ".pheno_pop.uneven.genotypes.mgf"),
            col_names = FALSE)

  pheno_qtl %>%
    write_csv(paste0("gpsm_runs/", testname, "/", testname, ".pheno_pop.qtl_trajectories.csv"))
}
if ("rand" %in% analyses){
  rand_geno_even %>%
    t() %>%
    as.data.frame() %>%
    mutate(rs = rownames(.),
           a1 = "a1",
           a2 = "a2") %>%
    select(rs, a1, a2, everything()) %>%
    #filter(!is.na(V1)) %>%
    write_csv(paste0("generation_genotypes/", testname, "/", testname, ".rand_pop.even.genotypes.mgf"),
              col_names = FALSE)

  rand_geno_uneven %>%
    t() %>%
    as.data.frame() %>%
    mutate(rs = rownames(.),
           a1 = "a1",
           a2 = "a2") %>%
    select(rs, a1, a2, everything()) %>%
    #filter(!is.na(V1)) %>%
    write_csv(paste0("generation_genotypes/", testname, "/", testname, ".rand_pop.uneven.genotypes.mgf"),
              col_names = FALSE)

  rand_qtl %>%
    write_csv(paste0("gpsm_runs/", testname, "/", testname, ".rand_pop.qtl_trajectories.csv"))}

rep(sampler_even$gen, sampler_even$keepind) %>%
  as.data.frame() %>%
  write_tsv(paste0("generation_genotypes/", testname, "/", testname, ".even.generation_phenotypes.txt"),
            col_names = FALSE)

rep(sampler_uneven$gen, sampler_uneven$keepind) %>%
  as.data.frame() %>%
  write_tsv(paste0("generation_genotypes/", testname, "/", testname, ".uneven.generation_phenotypes.txt"),
            col_names = FALSE)



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
ggsave(paste0("gpsm_runs/", testname, "/figures/", testname, ".g_vg.trends.png"), width = 8, height = 12)
