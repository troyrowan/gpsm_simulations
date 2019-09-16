with open("/data/tnr343/gpsm_sims/GPSM_scenarios.csv", "r") as sc:
	header = sc.readline().strip().split(",")
	header.append("analyses")
	for xx in sc:
		run = xx.split(",")[0]
		with open("/data/tnr343/gpsm_sims/gpsm_runs/"+ run + "/" + run + ".config.R", "w") as config:
			scenario = xx.strip().split(",")
			scenario[0] = '\"'+ scenario[0]+'\"'
			scenario.append("""c("gv", "pheno", "rand")""")
			outstuff = [(" = ").join([header[yy], scenario[yy]]) for yy in list(range(0,len(header)))]
			[config.write(xx + "\n") for xx in outstuff]
