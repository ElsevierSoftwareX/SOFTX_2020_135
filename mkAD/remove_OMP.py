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
    i=0
    for line in allLines:
        if line.replace(' ','').startswith("C$OMPparallel") or line.replace(' ','').startswith("!$OMPparallel"):
          writeFile.write('      Tlevel_1=Tlevel_0\n');
        if line.startswith("C$OMP") or line.startswith("!$OMP") or line.startswith("C$ ") or line.startswith("!$ "):
          line = "      call dummy(\'"+line.strip()+"\')\n"
          i+=1
#          print  line.strip()
        writeFile.write(line);
    if i>0:
      print i ,"OMP replacements in "+ targetfile+"."  
else:
    print "OMP replacement error: " + targetfile +" does not exist!"
