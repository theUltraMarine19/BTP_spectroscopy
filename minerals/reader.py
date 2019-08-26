import os
import sys

start_min = 0
start_max = 20000
selected = [1,3,14,23,26]
all_files = os.listdir('.')
# print(all_files)
ctr = 1
for f in all_files:
	if f in ['desktop.ini', 'reader.py', 'data.txt', 'data1.txt']:
		continue
	rf = open(f, 'r')
	all_lines = rf.readlines()
	flag1 = flag2 = 0
	for i in range(2,len(all_lines)):
		if (float(all_lines[i].split(',')[0])) == 490.59:
			flag1 = i
		if (float(all_lines[i].split(',')[0])) == 950.34:
			flag2 = i
			break

	# if flag == 2:
	# 	print("found", f)
	# else :
	# 	print("not f", f)
	# if flag1 + flag2 != 0:
		# print((950.34 - 490.59)/(flag2 - flag1), ctr, flag2 - flag1)
	# print(f)
	min_w = float(all_lines[2].split(',')[0])
	second = float(all_lines[3].split(',')[0])
	max_w = float(all_lines[len(all_lines) - 1].split(',')[0])
	
	# print(min_w, max_w, second - min_w, (max_w - min_w)/(len(all_lines)-2))
	if min_w > start_min:
		start_min = min_w
	if max_w < start_max:
		start_max = max_w

	if ctr in selected:
		intensity_values = []
		for j in range(flag1, flag2+1):
			intensity_values.append(float(all_lines[j].split(',')[1]))

	
		# print(ctr, "========", intensity_values)

		filtered_values = []
		for iter in range(0, int(len(intensity_values)/4)):
			filtered_values.append(sum(intensity_values[4*iter:4*iter+3]))

		print(ctr, "========", filtered_values)
	ctr += 1

# print(start_min, start_max)


# intensity_values = [[] for x in range(len(all_files)-2)]
# ctr = 0
# include = []
# for f in all_files:
# 	if f in ['desktop.ini', 'reader.py', 'data.txt']:
# 		continue
# 	rf = open(f, 'r')
# 	# print(f)
# 	all_lines = rf.readlines()
# 	for i in range(2,len(all_lines)):
# 		if (float(all_lines[i].split(',')[0])) == 490.59:
# 				include.append(ctr)
# 		if float(all_lines[i].split(',')[0]) >= start_min and float(all_lines[i].split(',')[0]) <= start_max:
# 			intensity_values[ctr].append(float(all_lines[i].split(',')[1]))
			

# 	# print(intensity_values[ctr][1:10])
# 	ctr += 1

# print(include)
# for k in range(len(intensity_values)):
# 	if k in include:
# 		print(intensity_values[k][0:100])

