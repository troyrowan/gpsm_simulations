pops = ["gv_pop", "pheno_pop", "rand_pop"]
rule gpsm_testing:
	input:
		endpoint = expand("gpsm_runs/{run}/figures/{run}.{population}.replicate{replicate}.manhattan.png",
		run = "scenario14",
		population = pops,
		replicate = "1")



rule run_sims:
	input:
		param = "gpsm_runs/{run}/{run}.config.txt"
	params:
		genos = expand("generation_genotypes/{{run}}/{{run}}.{population}.replicate{{replicate}}.genotypes.mgf", population = pops)
	benchmark:
		"benchmarks/run_sims/{run}.replicate{replicate}.benchmark.txt"
	log:
		"logs/run_sims/{run}.replicate{replicate}.log"
	output:
		gen_phenotypes = "generation_genotypes/{run}/{run}.generation_phenotypes.txt",
		real_phenotypes = expand("generation_genotypes/{{run}}/{{run}}.{population}.replicate{{replicate}}.trait_phenotypes.txt", population = pops),
		genotypes = expand("generation_genotypes/{{run}}/{{run}}.{population}.replicate{{replicate}}.genotypes.mgf.gz", population = pops),
		true_qtl = "generation_genotypes/{run}/{run}.true_qtl.csv",
		genetic_gain = "gpsm_runs/{run}/figures/{run}.g_vg.trends.png"
	shell:
		"(Rscript 190909_run_gpsm_sim.R {input.param}; pigz {params.genos}) > {log}"

rule make_grm:
	input:
		genotypes = "generation_genotypes/{run}/{run}.{population}.replicate{replicate}.genotypes.mgf.gz",
		gen_phenotypes = "generation_genotypes/{run}/{run}.generation_phenotypes.txt"
	params:
		outdir = "gpsm_runs/{run}",
		oprefix = "{run}.{population}.replicate{replicate}.grm"
	benchmark:
		"benchmarks/make_grm/{run}.benchmark.txt"
	log:
		"logs/make_grm/{run}.log"
	output:
		"gpsm_runs/{run}/{run}.{population}.replicate{replicate}.grm.sXX.txt",
		"gpsm_runs/{run}/{run}.{population}.replicate{replicate}.grm.log.txt"
	shell:
		"(gemma98 -g {input.genotypes} -p {input.gen_phenotypes} -gk 2 -n 1 -outdir {params.outdir} -o {params.oprefix}) > {log}"

rule run_gpsm:
	input:
		genotypes = "generation_genotypes/{run}/{{run}.{population}.replicate{replicate}.genotypes.mgf.gz",
		gen_phenotypes = "generation_genotypes/{run}/{run}.generation_phenotypes.txt",
		grm = "gpsm_runs/{run}/{run}.{population}.replicate{replicate}.grm.sXX.txt"
	params:
		outdir = "gpsm_runs/{run}",
		oprefix = "{run}.{population}.replicate{replicate}.gpsm"
	benchmark:
		"benchmarks/run_gpsm/{run}.benchmark.txt"
	log:
		"logs/run_gpsm/{run}.log"
	output:
		"gpsm_runs/{run}/{run}.{population}.replicate{replicate}.gpsm.assoc.txt",
		"gpsm_runs/{run}/{run}.{population}.replicate{replicate}.gpsm.log.txt"
	shell:
		"(gemma98 -g {input.genotypes} -p {input.gen_phenotypes} -k {input.grm} -n 1 -maf 0 -outdir {params.outdir} -o {params.oprefix}) > {log}"

rule plot_manhattans:
	input:
		assoc = "gpsm_runs/{run}/{run}.{population}.replicate{replicate}.gpsm.assoc.txt",
		true_qtl = "generation_genotypes/{run}/{run}.true_qtl.csv"
	params:
		run = "{run}",
		pop = "{population}",
		rep = "{replicate}"
	benchmark:
		"benchmarks/plot_manhattans/{run}.{population}.replicate{replicate}.benchmark.txt"
	log:
		"logs/plot_manhattans/{run}.{population}.replicate{replicate}.log"
	output:
		qq = "gpsm_runs/{run}/figures/{run}.{population}.replicate{replicate}.qq.png",
		manhattan = "gpsm_runs/{run}/figures/{run}.{population}.replicate{replicate}.manhattan.png"
	shell:
		"(Rscript manhattan_plotting.R {params.run} {params.pop} {params.replicate}) > {log}"
