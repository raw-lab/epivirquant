#!/usr/bin/env python

import re
from pathlib import Path
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

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
	fig5 = dict()
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

							# Define line style based on parameters
							#line_ = None
							line_width = 1
							if panel  == "A":
								if sig == "0.50":
									line_ = dict(color='blue', width=line_width, dash='solid')
								if sig == "0.75":
									line_ = dict(color='green', width=line_width, dash='dot')
								if sig == "1.00":
									line_ = dict(color='red', width=line_width, dash='dash')
							if panel  == "B":
								if fsize == "3":
									line_ = dict(color='purple', width=line_width, dash='solid')
								if fsize == "5":
									line_ = dict(color='orange', width=line_width, dash='dash')
								if fsize == "7":
									line_ = dict(color='black', width=line_width, dash='dot')
							if panel  == "C":
								if method in ["tPSF", "gam"]:
									line_ = dict(color='blue', width=line_width, dash='solid')
								if method in ["cs2", "gau"]:
									line_ = dict(color='red', width=line_width, dash='dash')
								if method in ["cs3", "hyb"]:
									line_ = dict(color='black', width=line_width, dash='dot')


							# Size plot
							label = f"{tool}-{panel}-{method}-f{fsize}-s{sig}"
							if panel not in figPanel:
								figPanel[panel] = go.Figure()
							figPanel[panel].add_trace(go.Scatter(x=df["iteration"], y=df["mu"], name=label, line=line_))
							#sns.lineplot(x='iteration', y='mu', data=df, label=label)

							# Count plot
							label = f"{tool}-{panel}-{method}-f{fsize}-s{sig}"
							if panel not in figCount:
								figCount[panel] = go.Figure()
							figCount[panel].add_trace(go.Scatter(x=df["iteration"], y=df["count"], name=label, line=line_))

							# Figure 5 plot
							if tool not in fig5:
								fig5[tool] = make_subplots(rows=3, cols=2, subplot_titles=("a", "b", "c", "d", "e", "f"))
							row_ = ["A", "B", "C"].index(panel) + 1
							fig5[tool].add_trace(go.Scatter(x=df["iteration"], y=df["mu"], line=line_), row=row_, col=1)
							fig5[tool].add_trace(go.Scatter(x=df["iteration"], y=df["count"], line=line_), row=row_, col=2)

							# seaborn plot	

	for p,f in figPanel.items():
		f.update_layout(title=f'Panel {p}: Comparison of Sizes by Category', xaxis_title='Iteration', yaxis_title='Size (nm) average', legend_title='Category')
		f.write_html(f'comparison_size_panel{p}.html')
		f.write_image(f'comparison_size_panel{p}.jpg', scale=2)
	for p,f in figCount.items():
		f.update_layout(title=f'Panel {p}: Comparison of Counts by Category', xaxis_title='Iteration', yaxis_title='Object Count', legend_title='Category')
		f.write_html(f'comparison_count_panel{p}.html')
		f.write_image(f'comparison_count_panel{p}.jpg', scale=2)
	for t,f in fig5.items():
		f.update_layout(
			template = "plotly_white",
			height=900, width=800,
			showlegend=False,
			font=dict(size=20, color="Black"),
			xaxis_range=[0,300],
			)
		f.update_xaxes(
			gridcolor="rgba(00, 00, 00, 0.00)",
			#tickangle=90
			)
		f.update_yaxes(
			showline=True,
			linewidth=2,
			linecolor='black',
			gridcolor="rgba(38,38,38,0.15)",
			title=dict(text='Percentage of proteins')
			)
		f.update_xaxes(
			zeroline=True,
			zerolinewidth=2,
			zerolinecolor='black'
			)
		f.update_yaxes(
			zeroline=True,
			zerolinewidth=2,
			zerolinecolor='black'
			)

		for i, annotation in enumerate(f['layout']['annotations']):
			if i % 2 == 0:
				annotation['x'] = 0.025
			else:
				annotation['x'] = 0.575
		#f.layout.annotations[0].update(x=0.025)
		#f.update_annotations(xanchor='left')
		#f.update_layout_annotations(xanchor='left')
		f.write_html(f'{t}/figure_5.html')
		f.write_image(f'{t}/figure_5.jpg', scale=2, height=1280, width=960)
	# Set title and labels
	plt.title('Comparison of Values by Category')
	plt.xlabel('Iteration')
	plt.ylabel('Size (nm) average')

	# Show the plot
	plt.tight_layout()
	#plt.savefig("comparison_plot.png")
	return 0

if __name__ == "__main__":
	exit(main())
