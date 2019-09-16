source("plotting_functions.R")
library(ggplot2)
library(stringr)
#Arg 1 = Scenario row_number
#Arg 2 = Selected population
#Arg 3 = Replicate


args <- commandArgs(TRUE)

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
selection_type = choose_name(args[2])
replicate = paste0("replicate", arg[3])

true_qtl = read_csv(paste0("generation_genotypes/", arg[1], "/", arg[1], ".true_qtl.csv"))
gpsm = readgwas(paste0("gpsm_runs/", arg[1], "/", arg[1], ".", arg[2], ".", replicate, ".gpsm.assoc.txt"))

ggqq(gpsm$p_score)+
	ggtitle(paste("Scenario", scenario_number , selection_type))
ggsave(paste0("gpsm_runs/", arg[1], "/figures/", arg[1], ".", replicate, ".gpsm.qq.png"),
width = 8,
height = 8)



ggmanhattan(gpsm,
						value = q,
						sigsnps = true_qtl$rs)+
							ggtitle(paste("Scenario", scenario_number , selection_type))
ggsave(paste0("gpsm_runs/", arg[1], "/figures/", arg[1], ".", replicate, ".gpsm.manhattan.png"),
width = 11,
height = 8)
