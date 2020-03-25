#!/usr/bin/python

def contains(theString, theQueryValue):
   return theString.find(theQueryValue) > -1

def modifyFile(file):
   if os.path.exists(os.getcwd() + "/" + file):
      openFile = open(file)
      writeFile = open(file+".new","w")
      allLines = openFile.readlines()
      lineno = 0
      found=False
      for line in allLines:
         lineno = lineno+1
         if contains(line.lower().strip(), "use g_") or contains(line.lower().strip(), "include \'g_"):
            line=line.split(",")[0]
            newLine = line.replace("g_","")
            line = newLine+"\n" + line+"\n"
            print "Adding include in " + file + ", string " + newLine
            found=True
         elif contains(line.lower().strip(), "use gps_") or contains(line.lower().strip(), "include \'gps_"):
            line=line.split(",")[0]
            newLine = line.replace("gps_","")
            line = newLine+"\n" + line+"\n"
            print "Adding include in " + file + ", string " + newLine
            found=True
         elif (contains(line.lower().strip(), "use ") or contains(line.lower().strip(), "include ")) and contains(line.lower().strip(),"_ad"):
            line=line.split(",")[0]
            newLine = line.lower().replace("_ad","")
            line = newLine+"\n" + line+"\n"
            print "Adding include in " + file + ", string " + newLine
            found=True
         elif found==True and contains(line.lower().strip(),"&") :
            line=""
         else:
            found=False
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

