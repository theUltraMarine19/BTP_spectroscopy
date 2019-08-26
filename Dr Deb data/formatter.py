import os
import sys
import csv

file = open("1A-tissue-100X-600-D00 20x20x1-15x2sec-2-r.txt", 'r')
all_lines = file.readlines()
ctr = 0
data = [[] for i in range(len(all_lines)-1)]
final_rows = [[] for i in range(len(all_lines))]
for i in range(len(all_lines)):
	if i==0:
		wavelens = [float(x) for x in all_lines[i].strip().split('\t')]
		header = [1, 2]
		header.extend(wavelens)
		final_rows[0] = header
	else:
		tmp = all_lines[i].strip().split('\t')
		data[ctr] = [float(x) for x in tmp[2:]]
		posx = float(tmp[0])
		posy = float(tmp[1])
		final_rows[ctr+1] = [posx, posy]
		final_rows[ctr+1].extend(data[ctr])
		ctr += 1

with open("paraffin_skin2.csv", 'w') as writeFile:
	writer = csv.writer(writeFile)
	writer.writerows(final_rows)

file.close()
writeFile.close()
