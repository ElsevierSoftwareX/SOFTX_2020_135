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
      sname=re.compile("^subroutine (\w+)\(")
      s=["asdffgh"]
      for line in allLines:
         lineno = lineno+1
         if line.lower().strip().startswith("subroutine"):
             s=sname.findall(line.strip().lower())
             l=line.split('(')
             print s,l
             if(contains(l[1],s[0])):
                 found=True
                 line=l[0]+"("+l[1].replace(s[0],s[0]+"v")
         elif not found and contains(line,s[0]):
             if contains(line,s[0]):
                 found=True
                 line=line.replace(s[0],s[0]+"v")
         elif found and contains(line.lower(),"implicit none"):
             line=line+"  double precision :: "+s[0]+"v\n"
         elif found:
             line=line.replace(s[0],s[0]+"v")
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

