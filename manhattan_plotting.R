source("plotting_functions.R")
library(ggplot2)
library(stringr)
library(readr)
library(dplyr)
library(tidyr)
library(cowplot)
library(qvalue)
library(viridis)
#Arg 1 = Scenario name
#Arg 2 = Selected population


arg <- commandArgs(TRUE)

choose_name <- function(xx){
	if (xx == "gv_pop"){
		selection = "True Breeding Value Selection"
	}
	if (xx == "pheno_pop"){
		selection = "Phenotypic Selection"
	}
	if (xx == "rand_pop"){
		selection = "Random Selection"
	}
return(selection)
}

scenario_number = str_replace(arg[1], "scenario", "")
selection_type = choose_name(arg[2])

true_qtl = read_csv(paste0("generation_genotypes/", arg[1], "/", arg[1], ".true_qtl.csv"), col_types = cols(rs = col_character()))
gpsm = readgwas(paste0("gpsm_runs/", arg[1], "/", arg[1], ".", arg[2], ".gpsm.assoc.txt")) %>% mutate(rs = paste(chr, pos, sep = ":"))

print(filter(gpsm, q < 0.1)$rs)
print(true_qtl$rs)

ggqq(gpsm$p_score)+
	ggtitle(paste("Scenario", scenario_number , selection_type))
ggsave(paste0("gpsm_runs/", arg[1], "/figures/", arg[1], ".", arg[2], ".gpsm.qq.png"),
width = 8,
height = 8)



ggmanhattan(gpsm,
						value = q,
						sigsnps = true_qtl$rs)+
							ggtitle(paste("Scenario", scenario_number , selection_type))
ggsave(paste0("gpsm_runs/", arg[1], "/figures/", arg[1], ".", arg[2], ".gpsm.manhattan.png"),
width = 11,
height = 8)
