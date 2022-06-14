import numpy as np
import csv
import matplotlib.pyplot as plt


def fill_unique_chrms(chrms, n):
	
	def check(lst, n):
		if len(lst[0]) == n:
			return list(map(int,lst[1]))
		new_vals = [0]*n
		for i in range(n):
			if str(i+1) in lst[0]:
				new_vals[i] = lst[1][np.where(lst[0] == str(i+1))[0]][0]
			else:
				continue
		return list(map(int,new_vals))
	
	pos = check(chrms[0],n)
	neg = check(chrms[1],n)
	
	return [sum(x) for x in zip(pos, neg)]
		

def gral_table(lengths, contents, chrms, output):
	header = ["Replicon", "Size (bp)", "GC content (%)", "CDS", "tRNA", "rRNA"]
	ttrnas = fill_unique_chrms(chrms[0], len(lengths))
	rrnas = fill_unique_chrms(chrms[1], len(lengths))
	cdss = fill_unique_chrms(chrms[2], len(lengths))
	
	csv_file = open(output, "w")
	writer = csv.writer(csv_file)
	writer.writerow(header)
	
	for i in range(len(lengths)):
		avg = contents[i][0]/contents[i][1]
		line = map(str, [i+1, lengths[i], avg, cdss[i], ttrnas[i], rrnas[i]])
		writer.writerow(line)
		
	avg_total = np.sum([x[0] for x in contents])/np.sum([x[1] for x in contents])
	footer = map(str, ["Total", np.sum(lengths), avg_total, np.sum(cdss), np.sum(ttrnas), np.sum(rrnas)])
	
	writer.writerow(footer)
	csv_file.close()
	
	return
	
def cogs_classif(hist, output):
	
	header1 = ["","Cellular Processes and Signaling"]+[""]*9+["Information Storage and Processing"]+[""]*5+["Metabolism"]+[""]*7+["Poorly Characterized","","",""]
	header2 = ["Replicon","D","M","N","O","T","U","V","W","Y","Z","A","B","J","K","L","X","C","E","F","G","H","I","P","Q","R","S","Unclassified"]
	
	csv_file = open(output, 'w')
	writer = csv.writer(csv_file)
	writer.writerow(header1)
	writer.writerow(header2)
	
	for i in range(len(hist.columns) - 2):
		line = map(str,[i+1] + [hist[hist["COG Category"] == c]["chr"+str(i+1)].item() for c in header2[1:]])
		writer.writerow(line)
		
	footer = map(str,["Total"] + [hist[hist["COG Category"] == c]["Frequency"].item() for c in header2[1:]])
	writer.writerow(footer)
	csv_file.close()
	
	return
	
def draw_histogram(hist, output):
	colors = ["#697BB7", '#304E9D', '#6A9AB5', '#1D82B0', '#4CA9B5', '#1C6875', '#6DBFA4', '#1A936F', '#4EB160', '#1C7633', '#E34949', '#CD1B1B', '#AD5B9F', '#A3378C', '#8F76B4',
	          '#533D91', '#B5D27B', '#83AD29', '#BCB868', '#988F1C', '#EDC28A', '#B7822A', '#DD8950', '#C65F17', '#696969', '#999999', '#EAEAEA']
	ax = hist.plot.bar(x="COG Category", y="Frequency", rot=90, legend=False, color=colors)
	for p in ax.patches:
		ax.annotate(str(p.get_height()), (p.get_x()+p.get_width()/2, p.get_height()+5), ha='center', va='center')
	plt.tight_layout()
	ax.figure.savefig(output+"_COG_Histogram.png", dpi=150)
	cogs_classif(hist, output+"_COG_Classification.csv")
	
	return
