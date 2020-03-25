#!/usr/bin/python
# we need to check if the name of the routine is used as parameter, that could be a bug in tapenade and preprocessing

def contains(theString, theQueryValue):
	return theString.find(theQueryValue) > -1

def modifyFile(file):
	if os.path.exists(os.getcwd() + "/" + file):
		openFile = open(file)
		allLines = openFile.readlines()
		i=0
		found = 0
		lineno = 0
                subroutine_name = ""
                newname=""
                out=""
		for line in allLines:
			lineno = lineno+1
			if contains(line.lower(), "subroutine"):
				# cutting out comments with c and !
				if line.strip().startswith("!") or line.lower().startswith("c"):
					out=out+line
					continue
				# cutting out strings !
				if line.strip().startswith("'") or line.strip().startswith('"'):
                                    out=out+line
				    continue
				# find the name of th
				if line.strip().lower().startswith("subroutine"):
                                    spl=line.strip().lower().split("(")
                                    print spl
                                    spl[1].replace(")","")
                                    subroutine_name = spl[0].split(" ")[1]
                                    line="subroutine "+subroutine_name+"("
                                    print "Found subroutine ",subroutine_name,"("
                                    for s in spl[1].split(","):
                                        print s,subroutine_name
                                        if s.replace(")","")==subroutine_name:
                                            found=1
                                            newname=s.replace("_ad","")
                                            print newname
                                        line=line+s
                                    print line
                                    out=out+line+")\n"
				    continue
				# writing down the end of functions
				if contains(line.lower(), "end subroutine"):
					out=out+line
					continue
			elif found == 1:
                            if contains(line.lower(),subroutine_name):
                                line.replace(subroutine_name,newname)
			out=out+line
                if found == 1 :
		    writeFile = open(file+".new","w")
                    writeFile.write(out)

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

