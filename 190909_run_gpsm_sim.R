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
args <- commandArgs(TRUE)
#This reads in all of the arguments from the param file
source(args[1])
#Uses already-created founder population (as this takes quite some time)
#Commands below:
# founderPop = runMacs(nInd=2000,
#                      nChr=10,
#                      segSites=5000,
#                      species = "CATTLE")
founderPop = readRDS("founderpop.RDS")

SP = SimParam$new(founderPop)
SP$addTraitA(nQtlPerChr=nqtl, gamma = gamma, shape = shape)
SP$setGender("yes_sys")
SP$addSnpChip(nSnpPerChr = nsnps)
SP$setVarE(h2 = h2)

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
                  run = as.numeric(),
                  pheno_g = as.numeric(),
                  pheno_vg = as.numeric(),
                  gv_g = as.numeric(),
                  gv_vg = as.numeric(),
                  rand_g = as.numeric(),
                  rand_vg = as.numeric())

gv_geno = data.frame(matrix(ncol = chrom*nsnps, nrow = 0))
colnames(gv_geno)<-paste("SNP", seq(1,nsnps*chrom), sep = "_")
rand_geno = data.frame(matrix(ncol = chrom*nsnps, nrow = 0))
colnames(gv_geno)<-paste("SNP", seq(1,nsnps*chrom), sep = "_")
pheno_geno = data.frame(matrix(ncol = chrom*nsnps, nrow = 0))
colnames(pheno_geno)<-paste("SNP", seq(1,nsnps*chrom), sep = "_")

for(run in 0:(reps-1)){
  pheno_pop = newPop(founderPop)
  gv_pop = newPop(founderPop)
  rand_pop = newPop(founderPop)

  gv_qtl = data.frame(matrix(nrow=nqtl*chrom, ncol=gens-burnins))
  pheno_qtl =  data.frame(matrix(nrow=nqtl*chrom, ncol=gens-burnins))
  rand_qtl =  data.frame(matrix(nrow=nqtl*chrom, ncol=gens-burnins))

  for(generation in 0:gens){
    pops[gens*run + generation + 1, "gen"] = generation
    pops[gens*run + generation + 1, "run"] = run + 1

    if (generation == 0){
      pops[gens*run + generation + 1, "pheno_g"] = meanG(newPop(founderPop))
      pops[gens*run + generation + 1, "pheno_vg"] = varG(newPop(founderPop))
      pops[gens*run + generation + 1, "gv_g"] = meanG(newPop(founderPop))
      pops[gens*run + generation + 1, "gv_vg"] = varG(newPop(founderPop))
      pops[gens*run + generation + 1, "rand_g"] = meanG(newPop(founderPop))
      pops[gens*run + generation + 1, "rand_vg"] = varG(newPop(founderPop))
    }

    if (generation > 0 & generation <= burnins){
      pheno_pop = selectCross(pop=pheno_pop,
                              nFemale=females,
                              nMale=males,
                              use="pheno",
                              nCrosses=crosses) #specify selection intensity, selection criteria, number of offspring
      pops[gens*run + generation + 1, "pheno_g"] = meanG(pheno_pop) #Updates genetic values in dataframe
      pops[gens*run + generation + 1, "pheno_vg"] = varG(pheno_pop) #Updates genetic values in dataframe

      gv_pop = selectCross(pop=gv_pop,
                           nFemale=females,
                           nMale=males,
                           use="pheno",
                           nCrosses=crosses)
      pops[gens*run + generation + 1, "gv_g"] = meanG(gv_pop) #Updates genetic values in dataframe
      pops[gens*run + generation + 1, "gv_vg"] = varG(gv_pop) #Updates genetic values in dataframe

      rand_pop = selectCross(pop=rand_pop,
                             nFemale=females,
                             nMale=males,
                             use="pheno",
                             nCrosses=crosses)
      pops[gens*run + generation + 1, "rand_g"] = meanG(rand_pop) #Updates genetic values in dataframe
      pops[gens*run + generation + 1, "rand_vg"] = varG(rand_pop) #Updates genetic values in dataframe
    }

    if (generation > burnins){
      pheno_pop = selectCross(pop=pheno_pop,
                              nFemale=females,
                              nMale=males,
                              use="pheno",
                              nCrosses=crosses) #specify selection intensity, selection criteria, number of offspring
      pops[gens*run + generation + 1, "pheno_g"] = meanG(pheno_pop) #Updates genetic values in dataframe
      pops[gens*run + generation + 1, "pheno_vg"] = varG(pheno_pop) #Updates genetic values in dataframe


      gv_pop = selectCross(pop=gv_pop, #Sets genomic EBV via RRBLUP for selection
                           nFemale=females,
                           nMale=males,
                           use="gv",
                           nCrosses=crosses)
      pops[gens*run + generation + 1, "gv_g"] = meanG(gv_pop)
      pops[gens*run + generation + 1, "gv_vg"] = varG(gv_pop)

      rand_pop = selectCross(pop=rand_pop,
                             nFemale=females,
                             nMale=males,
                             use="rand",
                             nCrosses=crosses)
      pops[gens*run + generation + 1, "rand_g"] = meanG(rand_pop)
      pops[gens*run + generation + 1, "rand_vg"] = varG(rand_pop)

      gv_geno = bind_rows(gv_geno,
                          as.data.frame(pullSnpGeno(gv_pop)) %>%
                            filter(row_number() %% pulleach == 0))
      gv_qtl[,(generation-burnins)] <- pullQtlGeno(gv_pop) %>% as.data.frame() %>% colSums()/(2*crosses) %>% t() %>% as.vector()

      pheno_geno = bind_rows(pheno_geno,
                             as.data.frame(pullSnpGeno(pheno_pop)) %>%
                               filter(row_number() %% pulleach == 0))
      pheno_qtl[,(generation-burnins)] <- pullQtlGeno(pheno_pop) %>% as.data.frame() %>% colSums()/(2*crosses) %>% t() %>% as.vector()

      rand_geno = bind_rows(rand_geno,
                            as.data.frame(pullSnpGeno(rand_pop)) %>%
                              filter(row_number() %% pulleach == 0))
      rand_qtl[,(generation-burnins)] <- pullQtlGeno(rand_pop) %>% as.data.frame() %>% colSums()/(2*crosses) %>% t() %>% as.vector()

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

      list(gv_pop, pheno_pop, rand_pop) %>%
        set_names(c("gv_pop", "pheno_pop", "rand_pop")) %>%
        iwalk(
          ~pheno(.x) %>%
            as.data.frame() %>%
            filter(row_number() %% pulleach == 0) %>%
            write_tsv(paste0("generation_genotypes/", testname, "/", testname, ".", .y, ".replicate", run+1, ".trait_phenotypes.txt"), append = TRUE))
    }
  }

  if ("gv" %in% analyses){
    gv_geno %>%
      t() %>%
      as.data.frame() %>%
      mutate(rs = rownames(.),
             a1 = "a1",
             a2 = "a2") %>%
      select(rs, a1, a2, everything()) %>%
      write_csv(paste0("generation_genotypes/", testname, "/", testname, ".gv_pop.replicate", run+1, ".genotypes.mgf"),
                col_names = FALSE)
  }

  if ("pheno" %in% analyses){
    pheno_geno %>%
      t() %>%
      as.data.frame() %>%
      mutate(rs = rownames(.),
             a1 = "a1",
             a2 = "a2") %>%
      select(rs, a1, a2, everything()) %>%
      write_csv(paste0("generation_genotypes/", testname, "/", testname, ".pheno_pop.replicate", run+1, ".genotypes.mgf"),
                col_names = FALSE)
  }
  if ("rand" %in% analyses){
    rand_geno %>%
      t() %>%
      as.data.frame() %>%
      mutate(rs = rownames(.),
             a1 = "a1",
             a2 = "a2") %>%
      select(rs, a1, a2, everything()) %>%
      filter(!is.na(V1)) %>%
      write_csv(paste0("generation_genotypes/", testname, "/", testname, ".rand_pop.replicate", run+1, ".genotypes.mgf"),
                col_names = FALSE)
  }

}

rep(1:(generation-burnins), each=crosses/pulleach) %>%
  as.data.frame() %>%
  write_tsv(paste0("generation_genotypes/", testname, "/", testname, ".generation_phenotypes.txt"),
            col_names = FALSE)

bind_cols(trait, gv_qtl) %>%
  write_csv("gpsm_runs/", testname, "/", testname, "gv_pop.qtl_trajectories.csv")

bind_cols(trait, pheno_qtl) %>%
  write_csv("gpsm_runs/", testname, "/", testname, "pheno_pop.qtl_trajectories.csv")

bind_cols(trait, rand_qtl) %>%
  write_csv("gpsm_runs/", testname, "/", testname, "rand_pop.qtl_trajectories.csv")
plot_grid(
  pops %>%
    select(gen, run, pheno_g, gv_g, rand_g) %>%
    melt(id = c("gen", "run")) %>%
    ggplot(aes(gen, value, color = variable))+
    geom_point(alpha = 0.3)+
    geom_smooth(size = 1.5)+
    labs(y = "Mean G", x = "Generation", main = "Genetic Gain") +
    scale_color_viridis_d(end = 0.75),
  #facet_grid(~burnin, scales = "free_x")

  pops %>%
    select(gen, run, pheno_vg, gv_vg, rand_vg) %>%
    melt(id = c("gen", "run")) %>%
    ggplot(aes(gen, value, color = variable))+
    geom_point(alpha = 0.3)+
    geom_smooth(size = 1.5)+
    labs(y = "Var G", x = "Generation", main = "Genetic Variance")+
    scale_color_viridis_d(end = 0.75),
  ncol = 1)
ggsave(paste0("gpsm_runs/", testname, "/figures/", testname, ".g_vg.trends.png"), width = 8, height = 12)
