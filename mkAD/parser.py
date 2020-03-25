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
    writeFile = open("diffsizes.f90","w")
    print "Begin writing diffsizes.mod"
    writeFile.write("      module diffsizes\n\n")
    writeFile.write("        implicit none\n\n")
    allLines = openFile.readlines()
    lines = set(allLines)
    for line in lines:
        if contains(line,"ISIZE"):
            #print line
            word_list = re.split(r'\s+',line)
            #print word_list[2]
            isize = word_list[2]
            i1 = 5
            i2 = isize.find("OF")
            col = isize[i1:i2]
            array = isize[i2+2:]
            print col+ "+"+array
            #ISIZEcolOFarray --> SIZE(array, col)
            print "integer "+word_list[2] +" = SIZE("+array+", "+col+")"
            writeFile.write("        integer "+word_list[2] +"\n")
    writeFile.write("      contains\n")
    writeFile.write("        subroutine init()\n\n")
    writeFile.write("          implicit none\n\n")
    for line in lines:
        if contains(line,"ISIZE"):
            #print line
            word_list = re.split(r'\s+',line)
            #print word_list[2]
            isize = word_list[2]
            i1 = 5
            i2 = isize.find("OF")
            col = isize[i1:i2]
            array = isize[i2+2:]
            print col+ "+"+array
            #ISIZEcolOFarray --> SIZE(array, col)
            print "init: "+word_list[2] +" = SIZE("+array+", "+col+")"
            writeFile.write("          "+word_list[2] +" = SIZE("+array+", "+col+")\n")
    writeFile.write("        end subroutine init\n")
    writeFile.write("      end module diffsizes\n")
    openFile.close()
    writeFile.close()
else:
    print "Data " + targetfile +" does not exist!"
