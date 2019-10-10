from sys import argv
import os


script, file = argv
#GPSM_scenarios_1_72
#command =  "mkdir gpsm_runs/scenario{" + file.split("_")[2] + ".." + file.split("_")[3].strip(".csv") + "}; mkdir gpsm_runs/scenario{" + file.split("_")[2] + ".." + file.split("_")[3].strip(".csv") + "}/figures/; mkdir generation_genotypes/scenario{" + file.split("_")[2] + ".." + file.split("_")[3].strip(".csv") + "}"
# print(command)
# os.system(command)

with open(file, "r") as sc:
	header = sc.readline().strip().split(",")
	header.append("analyses")
	header.append("samplingname")
	for xx in sc:
		run = xx.split(",")[0]
		with open("/data/tnr343/gpsm_sims/gpsm_runs/"+ run + "/" + run + ".config.R", "w") as config:
			scenario = xx.strip().split(",")
			scenario[0] = '\"'+ scenario[0]+'\"'
			scenario.append("""c("gv", "pheno", "rand")""")
			scenario.append('\"'+ scenario[11]+'\"')
			outstuff = [(" = ").join([header[yy], scenario[yy]]) for yy in list(range(0,len(header)))]
			config.writelines(["even_20gen = data.frame(gen = 1:20, keepind = 500)\n",
			"even_10gen = data.frame(gen = 1:10, keepind = 1000)\n",
			"even_5gen = data.frame(gen = 1:5, keepind = 2000)\n",
			"uneven_20gen = data.frame(gen = 1:20, keepind = c(12, 16, 16, 20, 24, 39, 43, 77, 105, 146, 160, 253, 350, 450, 652, 881, 1200, 1670, 2408, 1478))\n",
			"uneven_10gen = data.frame(gen = 1:10, keepind = c(14, 25, 27, 76, 190, 423, 819, 1696, 3604, 3126))\n",
			"uneven_5gen = data.frame(gen = 1:5, keepind = c(48, 144, 560, 1919, 7329))\n"])
			[config.write(xx + "\n") for xx in outstuff]
