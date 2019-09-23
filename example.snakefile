pops = ["gv_pop", "pheno_pop", "rand_pop"]

rule gpsm_testing:
	input:
		manhattan = expand("gpsm_runs/{run}/figures/{run}.{population}.gpsm.manhattan.png",
		run = ["scenario" + str(xx) for xx in list(range(73,185))],
		population = pops)

rule run_sims:
	input:
		param = "gpsm_runs/{run}/{run}.config.R",
		founderpop = "founderpop.RDS"
	params:
		genos = expand("generation_genotypes/{{run}}/{{run}}.{population}.genotypes.mgf", population = pops)
	output:
		gen_phenotypes = "generation_genotypes/{run}/{run}.generation_phenotypes.txt",
		real_phenotypes = expand("generation_genotypes/{{run}}/{{run}}.{population}.trait_phenotypes.txt", population = pops),
		genotypes = expand("generation_genotypes/{{run}}/{{run}}.{population}.genotypes.mgf.gz", population = pops),
		true_qtl = "generation_genotypes/{run}/{run}.true_qtl.csv",
		genetic_gain = "gpsm_runs/{run}/figures/{run}.g_vg.trends.png",
		snps = "generation_genotypes/{run}/{run}.snp.annotation.txt",
		qtl_trajectories = expand("gpsm_runs/{{run}}/{{run}}.{population}.qtl_trajectories.csv", population = pops)
	shell:
		"Rscript 190916_run_gpsm.R {input.param}; pigz {params.genos}"

rule make_grm:
	input:
		genotypes = "generation_genotypes/{run}/{run}.{population}.genotypes.mgf.gz",
		gen_phenotypes = "generation_genotypes/{run}/{run}.generation_phenotypes.txt"
	params:
		grm = "gpsm_runs/{run}/{run}.{population}.grm.sXX.txt",
		outdir = "gpsm_runs/{run}",
		oprefix = "{run}.{population}.grm"
	output:
		grm = temp("gpsm_runs/{run}/{run}.{population}.grm.sXX.txt.gz"),
		log = "gpsm_runs/{run}/{run}.{population}.grm.log.txt"
	shell:
		"gemma98 -g {input.genotypes} -p {input.gen_phenotypes} -gk 2 -n 1 -outdir {params.outdir} -o {params.oprefix}; pigz {params.grm}"

rule run_gpsm:
	input:
		genotypes = "generation_genotypes/{run}/{run}.{population}.genotypes.mgf.gz",
		gen_phenotypes = "generation_genotypes/{run}/{run}.generation_phenotypes.txt",
		grm = "gpsm_runs/{run}/{run}.{population}.grm.sXX.txt.gz",
		snps = "generation_genotypes/{run}/{run}.snp.annotation.txt"
	params:
		outdir = "gpsm_runs/{run}",
		oprefix = "{run}.{population}.gpsm"
	output:
		"gpsm_runs/{run}/{run}.{population}.gpsm.assoc.txt",
		"gpsm_runs/{run}/{run}.{population}.gpsm.log.txt"
	shell:
		"gemma98 -g {input.genotypes} -p {input.gen_phenotypes} -k {input.grm} -a {input.snps} -n 1 -lmm 4 -maf 0 -outdir {params.outdir} -o {params.oprefix}"

rule plot_manhattans:
	input:
		assoc = "gpsm_runs/{run}/{run}.{population}.gpsm.assoc.txt",
		true_qtl = "generation_genotypes/{run}/{run}.true_qtl.csv"
	params:
		run = "{run}",
		pop = "{population}"
	output:
		qq = "gpsm_runs/{run}/figures/{run}.{population}.gpsm.qq.png",
		manhattan = "gpsm_runs/{run}/figures/{run}.{population}.gpsm.manhattan.png"
	shell:
		"Rscript manhattan_plotting.R {params.run} {params.pop}"
