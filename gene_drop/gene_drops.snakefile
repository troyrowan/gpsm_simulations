rule target:
	input:
		#expand("run{run}/RAN.genedrop.run{run}.gpsm.assoc.txt", run = list(range(1,11)))
		expand("run{run}/RAN.genedrop.run{run}.gpsm.assoc.txt", run = list(range(41,51)))
		#expand("founderhaps/RAN_founderHaps.fullpop.chr{chr}.rds", chr = list(range(1,11)))

#rule create_founders:
#	input:
#		script = "create_founders.R"
#	params:
#		chrom = "{chr}"
#	output:
#		haps = "founderhaps/RAN_founderHaps.fullpop.chr{chr}.rds"
#	shell:
#		"Rscript create_founders.R {params.chrom}"

rule run_genedrop:
	input:
		#founderpop = "founderhaps/RAN_founderHaps.truesnps.7601.rds"
		founderpop = "founderhaps/RAN_founderHaps.200K.7601.rds"
		#founderpop = "founderhaps/RAN_founderHaps.50K.7601.rds"
		#founderpop = "founderhaps/RAN_founderHaps.50K.halfpop.rds"
		#founderpop = "founderhaps/RAN_founderHaps.fullpop.rds"
	params:
		run = "{run}"
	output:
		mgf = "run{run}/RAN.genedrop.run{run}.mgf",
		phen = "run{run}/RAN.genedrop.run{run}.phenotypes.txt",
		snps = "run{run}/RAN.genedrop.run{run}.snps.txt"
	shell:
		"Rscript run_genedrop.R {input.founderpop} {params.run}"
rule build_grm:
	input:
		mgf = "run{run}/RAN.genedrop.run{run}.mgf",
		phen = "run{run}/RAN.genedrop.run{run}.phenotypes.txt",
		snps = "run{run}/RAN.genedrop.run{run}.snps.txt"
	params:
		oprefix = "RAN.genedrop.run{run}.grm",
		outdir = "run{run}"
	output:
		grm = temp("run{run}/RAN.genedrop.run{run}.grm.sXX.txt")
	shell:
		"~/gemma0.98 -g {input.mgf} -p {input.phen} -a {input.snps} -gk 2 -outdir {params.outdir} -o {params.oprefix}"
rule run_gpsm:
	input:
		mgf = "run{run}/RAN.genedrop.run{run}.mgf",
		phen = "run{run}/RAN.genedrop.run{run}.phenotypes.txt",
		snps = "run{run}/RAN.genedrop.run{run}.snps.txt",
		grm = "run{run}/RAN.genedrop.run{run}.grm.sXX.txt"
	params:
		oprefix = "RAN.genedrop.run{run}.gpsm",
		outdir = "run{run}"
	output:
		assoc = "run{run}/RAN.genedrop.run{run}.gpsm.assoc.txt"
	shell:
		"~/gemma0.98 -g {input.mgf} -p {input.phen} -a {input.snps} -k {input.grm} -lmm 4 -outdir {params.outdir} -o {params.oprefix}"

rule run_lm:
        input:
                mgf = "run{run}/RAN.genedrop.run{run}.mgf",
                phen = "run{run}/RAN.genedrop.run{run}.phenotypes.txt",
                snps = "run{run}/RAN.genedrop.run{run}.snps.txt",
                grm = "run{run}/RAN.genedrop.run{run}.grm.sXX.txt"
        params:
                oprefix = "RAN.genedrop.run{run}.lm",
                outdir = "run{run}"
        output:
                assoc = "run{run}/RAN.genedrop.run{run}.lm.assoc.txt"
        shell:
                "~/gemma0.98 -g {input.mgf} -p {input.phen} -a {input.snps} -lm 4 -outdir {params.outdir} -o {params.oprefix}"
