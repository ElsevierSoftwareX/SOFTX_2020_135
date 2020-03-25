#!/usr/bin/python
import re

def contains(theString, theQueryValue):
   return theString.find(theQueryValue) > -1

def modifyFile(file):
   if os.path.exists(os.getcwd() + "/" + file):
      openFile = open(file)
      writeFile = open(file+".new","w")
      allLines = openFile.readlines()
      lineno = 0
      found=False
      functions=[]
      lines=[]
      notdefined=[]
      for line in allLines:
         lineno = lineno+1
         if re.search(r'^double precision :: temp$',line.lower().strip()) or re.search(r'double precision :: temp_ad$',line.lower().strip()):
            found=True
            print file," found temp"
            break
      for line in allLines:
         if found:
           if (contains(line.lower().strip(),"temp") and not re.search(r'temp_ad\d\(',line.lower().strip())):
               line=line.replace("temp_ad","temporary_ad")
         adtemp_list=re.findall(r'temp_ad\d\(',line.lower().strip())
         for l in adtemp_list:
             line=line.replace(l,"temp_ad(")
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

