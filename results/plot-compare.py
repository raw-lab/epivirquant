#!/usr/bin/env python

import re
from pathlib import Path
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go

def main():

	re_path = re.compile(r"panel(\w)/([\w\d]+)(?:-fsize([\d]+))?(?:-sig([\d.]+))?-GSL-(\d+)")
	re_mu = re.compile(r"Mean size over \d+ images:\s+([\d.]+) nm", re.MULTILINE)
	re_count = re.compile(r"Objects identified per image:\s+([\d.]+)", re.MULTILINE)

	# Load the data
	data = dict()
	for tool in ["cyclops", "2025-08-30"]:
		for filename in sorted(Path(f"{tool}").glob("out-panel*/*/Step-3_Corr/corrLog.txt")):
			match = re_path.search(str(filename))
			if match:
				panel, method, fsize, sig, i = match.groups()
				if fsize == "9":  # Skip fsize 9 for cyclops
					continue
				mu = re_mu.search(filename.read_text())
				if mu:
					fsize = fsize if fsize else "Default"
					sig = sig if sig else "Default"
					if tool not in data:
						data[tool] = dict()
					if panel not in data[tool]:
						data[tool][panel] = dict()
					if method not in data[tool][panel]:
						data[tool][panel][method] = dict()
					if fsize not in data[tool][panel][method]:
						data[tool][panel][method][fsize] = dict()
					if sig not in data[tool][panel][method][fsize]:
						data[tool][panel][method][fsize][sig] = dict(
							iteration=list(),
							mu=list(),
							count=list()
						)
					data[tool][panel][method][fsize][sig]["iteration"] += [int(i)]
					data[tool][panel][method][fsize][sig]["mu"] += [float(mu.group(1))]
					count = re_count.search(filename.read_text())
					if count:
						data[tool][panel][method][fsize][sig]["count"] += [float(count.group(1))]
	# Create a boxplot
	sns.set_theme(style="whitegrid")
	plt.figure(figsize=(10, 6))
	figPanel = dict()
	figCount = dict()
	# Iterate through the data and plot each series
	for tool in data:
		with open(f"data-{tool}.tsv", "w") as f:
			print("panel", "method", "fsize", "sigma", "iteration", "mu_size", "mu_count", sep="\t", file=f)
			for panel in data[tool]:
				for method in data[tool][panel]:
					for fsize in data[tool][panel][method]:
						for sig in data[tool][panel][method][fsize]:
							df = data[tool][panel][method][fsize][sig]
							# Sort by iteration
							sorted_list = sorted(zip(df["iteration"], df["mu"], df["count"]))
							df = dict(
								iteration=[x[0] for x in sorted_list],
								mu=[x[1] for x in sorted_list],
								count=[x[2] for x in sorted_list])							
							# save to tsv
							for i in range(1, len(df["iteration"])):
								print(panel, method, fsize, sig, df["iteration"][i], df["mu"][i], df["count"][i], sep="\t", file=f)
							# Size plot
							label = f"{tool}-{panel}-{method}-f{fsize}-s{sig}"
							if panel not in figPanel:
								figPanel[panel] = go.Figure()
							figPanel[panel].add_trace(go.Scatter(x=df["iteration"], y=df["mu"], name=label))
							#sns.lineplot(x='iteration', y='mu', data=df, label=label)

							# Count plot
							label = f"{tool}-{panel}-{method}-f{fsize}-s{sig}"
							if panel not in figCount:
								figCount[panel] = go.Figure()
							figCount[panel].add_trace(go.Scatter(x=df["iteration"], y=df["count"], name=label))
							# seaborn plot	

	for p,f in figPanel.items():
		f.update_layout(title=f'Panel {p}: Comparison of Sizes by Category', xaxis_title='Iteration', yaxis_title='Size (nm) average', legend_title='Category')
		f.show()
	for p,f in figCount.items():
		f.update_layout(title=f'Panel {p}: Comparison of Counts by Category', xaxis_title='Iteration', yaxis_title='Object Count', legend_title='Category')
		f.show()
	# Set title and labels
	plt.title('Comparison of Values by Category')
	plt.xlabel('Iteration')
	plt.ylabel('Size (nm) average')

	# Show the plot
	plt.tight_layout()
	plt.savefig("comparison_plot.png")
	return 0

if __name__ == "__main__":
	exit(main())
