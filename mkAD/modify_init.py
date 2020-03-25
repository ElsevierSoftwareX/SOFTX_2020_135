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
    print "Begin parsing head file to include init call"
    allLines = openFile.readlines()
    lines = set(allLines)
    search = 0
    for line in allLines:
      if(search == 2):
            # schon gefunden
            writeFile.write(line)
      else:
        if search == 0:
          if contains(line,"IMPLICIT NONE"):
            # we begin search for final INCLUDE
            writeFile.write(line)
            search = 1
          else:
            writeFile.write(line)
        else:
          if search == 1:
            if line.startswith("!"):
              # just a comment
              writeFile.write(line)
            else:
              if contains(line,"INCLUDE") or contains(line,"INTEGER") or contains(line,"DOUBLE") or contains(line,"EXTERNAL"):
                # just an include
                writeFile.write(line)
              else:
                # HERE we introduce the call
                writeFile.write("  CALL INIT()\n");
                writeFile.write(line)
                search = 2
    openFile.close()
    writeFile.close()
else:
    print "Data " + targetfile +" does not exist!"
