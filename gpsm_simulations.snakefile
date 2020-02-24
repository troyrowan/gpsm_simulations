pops = ["gv_pop",  "rand_pop"]
rule gpsm_testing:
	input:
		# endpoint = expand("gpsm_runs/{run}/{run}.{population}.{sampling}.qtl_trajectories.csv",
		# run = ["scenario" + str(xx) for xx in list(range(1,274))],
		# population = pops)
		# manhattan = expand("gpsm_runs/{run}/figures/{run}.{population}.{sampling}.qtl_trajectories.png",
		# run = ["scenario" + str(xx) for xx in list(range(100, 187))],
		# population = pops)
		manhattan = expand("gpsm_runs/{run}/rep{rep}/{run}.{population}.{sampling}.rep{rep}.gpsm.assoc.txt",
		run = ["scenario" + str(xx) for xx in list(range(370,388))],
		population = pops,
		rep = list(range(1,11)),
		sampling = ["even", "uneven"])
		#sims = expand("generation_genotypes/{run}/rep{rep}/{run}.{population}.{sampling}.rep{rep}.genotypes.mgf.gz",
		#run = ["scenario" + str(xx) for xx in list(range(997, 1000))],
		#population = pops,
		#rep = list(range(1,6)),
		#sampling = ["even", "uneven"])
		# manhattan = expand("generation_genotypes/{run}/{run}.{population}.{sampling}.genotypes.mgf.gz",
		# run = ["scenario" + str(xx) for xx in list(range(334, 370))],
		# population = pops,
		# sampling = ["even", "uneven"])
		# manhattan = expand("gpsm_runs/{run}/{run}.{population}.{sampling}.gpsm.assoc.txt",
		# run = ["scenario" + str(xx) for xx in list(range(334,370))],
		# population = pops,
		# sampling = ["even", "uneven"])


#Probably not the best way of doing most this: I skipped a step or two by manually making my param files and founder pop outside of Snakemake, but they could easily be their own rule.

rule run_sims:
	input:
		#founderpop = "100K_cattle_founderpop.RDS",
		param = "gpsm_runs/{run}/rep{rep}/{run}.config.R", #Config files created from my running sheet of simulation scenarios
		founderpop = "founderpops/founderpop{rep}.RDS" #Manually created cattle founder population with 10 chromosomes, 1K SNPs per chromosome
	priority: 100 #Priority ensures that all iterations of this rule will be run before any subsequent rules
	params:#Example of param where this is an even more intermediate file
		genos = expand("generation_genotypes/{{run}}/rep{{rep}}/{{run}}.{population}.{sampling}.rep{{rep}}",
		population = pops,
		sampling = ["even", "uneven"]) #Example of an expand command
	output:
		pedfile = temp(expand("generation_genotypes/{{run}}/rep{{rep}}/{{run}}.{population}.{sampling}.rep{{rep}}.ped",
						population = pops,
						sampling = ["even", "uneven"])),
		mapfile = temp("generation_genotypes/{run}/rep{rep}/{run}.rep{rep}.map"),
		#true_qtl = "generation_genotypes/{run}/rep{rep}/{run}.true_qtl.csv",
		genetic_gain = "gpsm_runs/{run}/rep{rep}/figures/{run}.rep{rep}.g_vg.trends.png",
		qtl_trajectories = expand("gpsm_runs/{{run}}/rep{{rep}}/{{run}}.{population}.rep{{rep}}.qtl_trajectories.csv",
						population = pops)
	shell:
		"Rscript Runsims.R {input.param}"#Notice how to wrap Stdout and point it to log file specified above
# 		#This script may be a bit funky in that it runs three simulations simultaneously (BV, phenotypic, and random selection, but on the same starting population), so this isn't a simple "one in one out"
rule convert_plink:
	input:
		pedfile = "generation_genotypes/{run}/rep{rep}/{run}.{population}.{sampling}.rep{rep}.ped",
		mapfile = "generation_genotypes/{run}/rep{rep}/{run}.rep{rep}.map"
	params:
		oprefix = "generation_genotypes/{run}/rep{rep}/{run}.{population}.{sampling}.rep{rep}"
	output:
		genotypes = expand("generation_genotypes/{{run}}/rep{{rep}}/{{run}}.{{population}}.{{sampling}}.rep{{rep}}.{suffix}",
						suffix = ["bed", "bim", "fam"]),
		extras = temp(expand("generation_genotypes/{{run}}/rep{{rep}}/{{run}}.{{population}}.{{sampling}}.rep{{rep}}.{suffix}",
						suffix = ["nosex", "log"]))
	shell:
		"plink --ped {input.pedfile} --map {input.mapfile} --make-bed --out {params.oprefix}"

rule make_grm:
	input:
		genotypes = "generation_genotypes/{run}/rep{rep}/{run}.{population}.{sampling}.rep{rep}.bed"
	threads: 8
	params:
		inprefix = "generation_genotypes/{run}/rep{rep}/{run}.{population}.{sampling}.rep{rep}",
		grm = "gpsm_runs/{run}/rep{rep}/{run}.{population}.{sampling}.rep{rep}.grm.sXX.txt", #Again this could be an output of the rule, but I choose to gzip in the same shell command, so it's no longer the final output of the rule
		outdir = "gpsm_runs/{run}/rep{rep}",
		oprefix = "{run}.{population}.{sampling}.rep{rep}.grm"
	output:
		grm = temp("gpsm_runs/{run}/rep{rep}/{run}.{population}.{sampling}.rep{rep}.grm.sXX.txt.gz"), #These files are giant and only needed for intermediate step, so wrapping in temp causes them to be deleted when they're no longer needed
		#log = "gpsm_runs/{run}/{run}.{population}.{sampling}.grm.log.txt"
	shell:
		"gemma98 -bfile {params.inprefix} -gk 2 -outdir {params.outdir} -o {params.oprefix}; pigz {params.grm}"

rule run_gpsm:
	input:
		grm = "gpsm_runs/{run}/rep{rep}/{run}.{population}.{sampling}.rep{rep}.grm.sXX.txt.gz",
		genotypes = "generation_genotypes/{run}/rep{rep}/{run}.{population}.{sampling}.rep{rep}.bed"
	threads: 8
	params:
		inprefix = "generation_genotypes/{run}/rep{rep}/{run}.{population}.{sampling}.rep{rep}",
		outdir = "gpsm_runs/{run}/rep{rep}",
		oprefix = "{run}.{population}.{sampling}.rep{rep}.gpsm"
	output:
		"gpsm_runs/{run}/rep{rep}/{run}.{population}.{sampling}.rep{rep}.gpsm.assoc.txt",
		"gpsm_runs/{run}/rep{rep}/{run}.{population}.{sampling}.rep{rep}.gpsm.log.txt"
	shell:
		"gemma98 -bfile {params.inprefix} -k {input.grm} -lmm 4 -maf 0.01 -outdir {params.outdir} -o {params.oprefix}"

rule generate_reports:
	input:
		assoc = expand("gpsm_runs/{{run}}/rep{rep}/{{run}}.{population}.{sampling}.rep{rep}.gpsm.assoc.txt",
		rep = list(range(1,11)),
		population = ["gv_pop", "rand_pop"],
		sampling = ["even", "uneven"], variable = ["age", "generation"]),
		#true_qtl = "generation_genotypes/{run}/{run}.true_qtl.csv"
	params:

	output:

	shell: #No actual input or output files in this shell command, but it still requires input files and will create stated outputfiles.
		#"(Rscript manhattan_plotting.R {params.run} {params.pop}) > {log}"
		"Rscript sampling_manhattan_plotting.R {params.run} {params.pop}"
