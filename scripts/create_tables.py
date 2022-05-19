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
	writer = csv.writer(csv_file, delimiter='\t')
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
