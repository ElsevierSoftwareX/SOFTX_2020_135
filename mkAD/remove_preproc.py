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
        if line.replace(' ','').lower().startswith('#include'):
          line = '      INCLUDE '+line.lstrip().lstrip('#').lstrip().split(None,1)[1].replace('\"',"\'")
        elif line.startswith("#"):
          line = "      call dummy(\'"+line.strip()+"\')\n"
#          print  line.strip()
          i+=1
        writeFile.write(line);
    if i>0:
      print i ,"preprocessor replacements in "+ targetfile+"."  
else:
    print "Preprocessor replacement error: " + targetfile +" does not exist!"
