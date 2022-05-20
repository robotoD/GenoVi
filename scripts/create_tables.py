import numpy as np
import csv


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
