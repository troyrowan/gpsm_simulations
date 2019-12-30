source("plotting_functions.R")
library(ggplot2)
library(stringr)
library(readr)
library(dplyr)
library(tidyr)
library(qvalue)
library(viridis)
library(reshape2)
library(cowplot)
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

true_qtl =
  read_csv(
    paste0("generation_genotypes/", arg[1], "/", arg[1], ".true_qtl.csv"),
           col_types = cols(rs = col_character()))
#Example of file names for reading:
#gpsm_runs/scenario335/scenario335.rand_pop.uneven.gpsm.assoc.txt
even_gpsm =
  readgwas(
    paste0("gpsm_runs/", arg[1], "/", arg[1], ".", arg[2], ".even.gpsm.assoc.txt")) %>%
  mutate(rs = paste(chr, pos, sep = ":"))

uneven_gpsm =
  readgwas(
    paste0("gpsm_runs/", arg[1], "/", arg[1], ".", arg[2], ".uneven.gpsm.assoc.txt")) %>%
  mutate(rs = paste(chr, pos, sep = ":"))


trajectories =
  read_csv(
		paste0("gpsm_runs/", arg[1] , "/", arg[1], ".", arg[2], ".qtl_trajectories.csv"),
           col_types = cols(rs = col_character())) %>%
  left_join(uneven_gpsm)

print(head(uneven_gpsm))

plot_grid(
ggqq(even_gpsm$p_score)+
	ggtitle(paste("Scenario", scenario_number , selection_type, "Even Sampling"))+
cowplot::theme_cowplot(),
ggqq(uneven_gpsm$p_score)+
	ggtitle(paste("Scenario", scenario_number , selection_type, "Uneven Sampling"))+
cowplot::theme_cowplot(),
nrow = 2)
ggsave(paste0("gpsm_runs/", arg[1], "/figures/", arg[1], ".", arg[2], ".gpsm.qq.png"),
width = 8,
height = 8)

plot_grid(
	ggmanhattan(even_gpsm,
						value = p,
						true_qtl = true_qtl)+
							ggtitle(paste("Scenario", scenario_number , selection_type, "Even Sampling")),
	ggmanhattan(uneven_gpsm,
							value = p,
							true_qtl = true_qtl)+
								ggtitle(paste("Scenario", scenario_number , selection_type, "Uneven Sampling")),
							nrow = 2)
ggsave(paste0("gpsm_runs/", arg[1], "/figures/", arg[1], ".", arg[2], ".gpsm.manhattan.png"),
width = 16,
height = 8)


plot_grid(
	trajectories %>%
	  top_n(., 20, abs(effect)) %>%
	  melt(id = c("chr", "pos", "effect", "rs", "af", "p_score", "q"), value.name = "AF") %>%
	  mutate(generation = as.numeric(str_remove(variable, "^X"))) %>%
	  ggplot(aes(x = generation, y = AF, group = rs, color = effect))+
	  geom_smooth(se = FALSE)+
	  scale_color_viridis_c()+
		ylim(c(0,1))+
	  labs(x = "Generation", title = paste("Scenario", scenario_number , selection_type, "Top QTL"))+
  cowplot::theme_cowplot(),
	trajectories %>%
		filter(-log10(q) > 1) %>%# only looking at significant or nearly significant loci
		top_n(., 20, -log10(q)) %>%
	  melt(id = c("chr", "pos", "effect", "rs", "af", "p_score", "q"), value.name = "AF") %>%
	  mutate(generation = as.numeric(str_remove(variable, "^X"))) %>%
	  ggplot(aes(x = generation, y = AF, group = rs, color = -log10(q)))+
	  geom_smooth(se = FALSE)+
	  scale_color_viridis_c()+
		ylim(c(0,1))+
	  labs(x = "Generation", title = paste("Scenario", scenario_number , selection_type, "Top GPSM Hits"))+
  cowplot::theme_cowplot(),
	nrow = 2)
ggsave(
	paste0("gpsm_runs/", arg[1], "/figures/", arg[1], ".", arg[2], ".qtl_trajectories.png"),
	width = 8,
	height = 11)
