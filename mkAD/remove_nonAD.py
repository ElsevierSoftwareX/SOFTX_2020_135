#!/usr/bin/python
# -*- coding: iso-8859-1 -*-

import sys
import os
import re
def contains(theString, theQueryValue):
  return theString.find(theQueryValue) > -1
targetfile = sys.argv[1]

if os.path.exists(os.getcwd() + "/" + targetfile):
    openFile = open(targetfile)
    writeFile = open(targetfile+".new","w")
    allLines = openFile.readlines()
    lines = set(allLines)
    found=False
    in_nonAD=False
    out=""
    counter=0
    for line in allLines:
        if line.startswith("#ifndef AD") or in_nonAD:
            found=True
            in_nonAD=True
            if (line.startswith("#if")):
                counter+=1
            elif(line.startswith("#endif")):
                counter-=1
            if (counter==0):
                counter=0
                in_nonAD=False
        else:
            out+=line
    writeFile.write(out);
    if found:
      print "#ifndef AD block deleted in "+ targetfile+"."  
else:
    print "#ifndef AD replacement error: " + targetfile +" does not exist!"
