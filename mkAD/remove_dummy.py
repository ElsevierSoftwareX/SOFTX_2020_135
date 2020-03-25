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
        if line.replace(' ','').lower().startswith("calldummy(\'") or line.replace(' ','').lower().startswith("calldummy_ori(\'"):
          word_list = re.split(r'\'',line)
          line = word_list[1]+"\n"
          i+=1
        if line.replace(' ','').lower() != 'tlevel_1=tlevel_0\n':
          writeFile.write(line);
    if i>0:
      print i ,"OMP/preprocessor restores in "+ targetfile+"."
else:
    print "OMP/preprocessor restore error: " + targetfile +" does not exist!"

