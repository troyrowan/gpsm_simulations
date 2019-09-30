from sys import argv
import os


script, file = argv
#GPSM_scenarios_1_72
command =  "mkdir gpsm_runs/scenario{" + file.split("_")[2] + ".." + file.split("_")[3].strip(".csv") + "}; mkdir gpsm_runs/scenario{" + file.split("_")[2] + ".." + file.split("_")[3].strip(".csv") + "}/figures/; mkdir generation_genotypes/scenario{" + file.split("_")[2] + ".." + file.split("_")[3].strip(".csv") + "}"
print(command)
os.system(command)

with open(file, "r") as sc:
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
