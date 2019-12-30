pops = ["gv_pop",  "rand_pop"]
rule gpsm_testing:
	input:
		# endpoint = expand("gpsm_runs/{run}/{run}.{population}.{sampling}.qtl_trajectories.csv",
		# run = ["scenario" + str(xx) for xx in list(range(1,274))],
		# population = pops)
		# manhattan = expand("gpsm_runs/{run}/figures/{run}.{population}.{sampling}.qtl_trajectories.png",
		# run = ["scenario" + str(xx) for xx in list(range(100, 187))],
		# population = pops)
		manhattan = expand("gpsm_runs/{run}/figures/{run}.{population}.gpsm.manhattan.png",
		run = ["scenario" + str(xx) for xx in list(range(334, 370))],
		population = pops,
		)
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
		param = "gpsm_runs/{run}/{run}.config.R", #Config files created from my running sheet of simulation scenarios
		founderpop = "founderpop.RDS" #Manually created cattle founder population with 10 chromosomes, 1K SNPs per chromosome
	threads: 10 #Example of how specifying threads works with Snakemake
	priority: 100 #Priority ensures that all iterations of this rule will be run before any subsequent rules
	params:#Example of param where this is an even more intermediate file
		genos = expand("generation_genotypes/{{run}}/{{run}}.{population}.{sampling}.genotypes.mgf",
		population = pops,
		sampling = ["even", "uneven"]) #Example of an expand command
	# benchmark: #Specifys where to write benchmarks for this rule
	# 	"benchmarks/run_sims/{run}.benchmark.txt"
	# log: #log file written here
	# 	"logs/run_sims/{run}.log"
	output:
		#gen_phenotypes = "generation_genotypes/{run}/{run}.{sampling}.generation_phenotypes.txt",
		##real_phenotypes = expand("generation_genotypes/{{run}}/{{run}}.{population}.{sampling}.trait_phenotypes.txt", population = pops), #expands over three populations, but uses whatever "run" environment tells it to
		genotypes = expand("generation_genotypes/{{run}}/{{run}}.{population}.{sampling}.genotypes.mgf.gz",
		population = pops,
		sampling = ["even", "uneven"]), #Final output, RScript will only output .mgf file (see params), final result after gzipping will be this file though
		# true_qtl = "generation_genotypes/{run}/{run}.true_qtl.csv",
		# genetic_gain = "gpsm_runs/{run}/figures/{run}.g_vg.trends.png",
		# snps = "generation_genotypes/{run}/{run}.snp.annotation.txt",
		# qtl_trajectories = expand("gpsm_runs/{{run}}/{{run}}.{population}.qtl_trajectories.csv",
		# population = pops)
	shell:#two shell commands happen here: first creates all output files, then gzips the outputted genotypes (referrred to in params)
		"Rscript 190916_run_gpsm.R {input.param}; pigz {params.genos}"#Notice how to wrap Stdout and point it to log file specified above
# 		#This script may be a bit funky in that it runs three simulations simultaneously (BV, phenotypic, and random selection, but on the same starting population), so this isn't a simple "one in one out"
rule make_grm:
	input:
		genotypes = "generation_genotypes/{run}/{run}.{population}.{sampling}.genotypes.mgf.gz",
		gen_phenotypes = "generation_genotypes/{run}/{run}.{sampling}.generation_phenotypes.txt"
	threads: 8
	params:
		grm = "gpsm_runs/{run}/{run}.{population}.{sampling}.grm.sXX.txt", #Again this could be an output of the rule, but I choose to gzip in the same shell command, so it's no longer the final output of the rule
		outdir = "gpsm_runs/{run}",
		oprefix = "{run}.{population}.{sampling}.grm"
	benchmark:
		"benchmarks/make_grm/{run}.{population}.{sampling}.benchmark.txt"
	log:
		"logs/make_grm/{run}.{population}.{sampling}.log"
	output:
		grm = temp("gpsm_runs/{run}/{run}.{population}.{sampling}.grm.sXX.txt.gz"), #These files are giant and only needed for intermediate step, so wrapping in temp causes them to be deleted when they're no longer needed
		log = "gpsm_runs/{run}/{run}.{population}.{sampling}.grm.log.txt"
	shell:
		"(gemma98 -g {input.genotypes} -p {input.gen_phenotypes} -gk 2 -n 1 -outdir {params.outdir} -o {params.oprefix}; pigz {params.grm}) > {log}"

rule run_gpsm:
	input:
		genotypes = "generation_genotypes/{run}/{run}.{population}.{sampling}.genotypes.mgf.gz",
		gen_phenotypes = "generation_genotypes/{run}/{run}.{sampling}.generation_phenotypes.txt",
		grm = "gpsm_runs/{run}/{run}.{population}.{sampling}.grm.sXX.txt.gz",
		snps = "generation_genotypes/{run}/{run}.snp.annotation.txt"
	threads: 8
	params:
		outdir = "gpsm_runs/{run}",
		oprefix = "{run}.{population}.{sampling}.gpsm"
	benchmark:
		"benchmarks/run_gpsm/{run}.{population}.{sampling}.benchmark.txt"
	log:
		"logs/run_gpsm/{run}.{population}.{sampling}.log"
	output:
		"gpsm_runs/{run}/{run}.{population}.{sampling}.gpsm.assoc.txt",
		"gpsm_runs/{run}/{run}.{population}.{sampling}.gpsm.log.txt"
	shell:
		"(gemma98 -g {input.genotypes} -p {input.gen_phenotypes} -k {input.grm} -a {input.snps} -n 1 -lmm 4 -maf 0.01 -outdir {params.outdir} -o {params.oprefix}) > {log}"

rule plot_manhattans:
	input:
		assoc = expand("gpsm_runs/{{run}}/{{run}}.{{population}}.{sampling}.gpsm.assoc.txt", sampling = ["even", "uneven"]),
		#true_qtl = "generation_genotypes/{run}/{run}.true_qtl.csv"
	params:
		run = "{run}",  #These strings are the only input that this script needs to find relevant inputs and plot them
		pop = "{population}"
	benchmark:
		"benchmarks/plot_manhattans/{run}.{population}.benchmark.txt"
	log:
		"logs/plot_manhattans/{run}.{population}.log"
	output:
		qq = "gpsm_runs/{run}/figures/{run}.{population}.gpsm.qq.png",
		manhattan = "gpsm_runs/{run}/figures/{run}.{population}.gpsm.manhattan.png",
		#qtls = "gpsm_runs/{run}/figures/{run}.{population}.qtl_trajectories.png"
	shell: #No actual input or output files in this shell command, but it still requires input files and will create stated outputfiles.
		#"(Rscript manhattan_plotting.R {params.run} {params.pop}) > {log}"
		"(Rscript sampling_manhattan_plotting.R {params.run} {params.pop}) > {log}"
