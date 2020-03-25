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
      functions=[]
      lines=[]
      notdefined=[]
      for line in allLines:
         lineno = lineno+1
         if contains(line.lower().strip(), "external "):
            try:
                funname=line.lower().split("::")[1]
            except:
                funname=line.lower().split(" ")[-1]
            func=funname.strip().split(",")
            functions.append(func)
            lines.append(lineno)
#            print "Adding ",func," from ",lineno, "to list"
     
#      print functions
#      print lines

      g_fun=re.compile('g_\w+\(')
      if lines:
        for line in allLines[lines[-1]:]:
          if (contains(line.lower().strip(),"g_")):
              fun=g_fun.findall(line.lower())
              for f in fun:
                  s=f.split("(")
#                  print fun,s
                  if [s[0]] not in functions and [s[0].strip("g_")] in functions and s[0] not in notdefined:
                      notdefined.append(s[0])
      if notdefined:
          print notdefined, "is not defined"
      for n in notdefined:
          allLines.insert(lines[-1],"\tdouble precision, external :: "+n+"\n")
      for l in allLines:
          writeFile.write(l)
                        
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

