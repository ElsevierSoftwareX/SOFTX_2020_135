#!/usr/bin/python

def contains(theString, theQueryValue):
	return theString.find(theQueryValue) > -1

def modifyFile(file):
	if os.path.exists(os.getcwd() + "/" + file):
		openFile = open(file)
		writeFile = open(file+".new","w")
		allLines = openFile.readlines()
		i=0
		found = 0
		lineno = 0
		for line in allLines:
			lineno = lineno+1
			if contains(line.lower(), "function"):
				# cutting out comments with c and !
				if line.strip().startswith("!") or line.lower().startswith("c"):
					writeFile.write(line)
					continue
				# cutting out strings !
				if line.strip().startswith("'") or line.strip().startswith('"'):
					writeFile.write(line)
					continue
				# writing down functions that already are defined correct
				if line.strip().lower().startswith("function"):
					writeFile.write(line)
					continue
				# writing down the end of functions
				if contains(line.lower(), "end function"):
					writeFile.write(line)
					continue
				# we have to do something
				else:
					word_list = re.split(" ", line); print "wl=",word_list 
					if len(word_list[0])> 0:
						continue
					#counting white spaces
					istart = 0
					iend = 0
					newLine = ""
					for i in range(0, len(word_list)):
						if len(word_list[i]) == 0:
							istart = istart+1
							newLine = newLine + " " + word_list[i]
						else:
							break
					iend = istart
					length = 0;
					for i in range(istart, len(word_list)):
						if word_list[i].lower() == "function":
							iend = i
							break
					iendstart = iend
					#special treatment of recursive functions
					if word_list[iend-1].lower() == "recursive":
						iendstart = iend -1
					type = word_list[istart]
					for i in range(istart+1, iendstart):
						type = type + " " + word_list[i]
					for i in range(iendstart, len(word_list)):
						newLine = newLine + " " + word_list[i]
					func = word_list[iend+1]
					func_list = re.split("\(", func)
					name = func_list[0]
					found = 1
					print "Doing tricks for file " + file + ", function " + name + ". (line " + str(lineno) + ")\n[",line, " is now ",newLine, " with type ",type
					line = newLine
			elif found == 1:
				if contains(line.lower(), "implicit"):
					found = 2
			elif found == 2:
				if not contains(line.lower(), "include"):
					line = "  " + type + " " + name + "\n" + line
					#attach whie spaces
					for i in range(0, istart):
						line = " " + line
					found = 0
			writeFile.write(line)
	else:
		print "Could not open file ",file

if __name__ == "__main__":
	import sys
	import os
	import re
	i = 0
	for arg in sys.argv:
		if i == 0:
			i = 1
		else:
			modifyFile(arg)

