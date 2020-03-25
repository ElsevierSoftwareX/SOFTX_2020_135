# MIT License
#
# Copyright (c) 2020 SHEMAT-Suite
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

#====================================================
#Generate version.inc
#====================================================
string(TIMESTAMP _configuration_time "%Y-%m-%d %H:%M:%S [UTC]" UTC)
configure_file(version.inc.in generated/version.inc @ONLY)
#====================================================


#====================================================
#This section is for updating forward/input/read_check.f90 with newly implemented reading paramters
#====================================================
#Read all read_*.f90 files and put the parameter names in @read_cur_use
FOREACH(fname ${read_files})
   file(STRINGS ${fname} cur_file REGEX "^[^c^C^!].*found.*(\#|key_char//') [a-zA-Z0-9 _]+.*")
   #message("${cur_file}")
   foreach(str ${cur_file})
      #      message("str=${str}")
      string(REGEX REPLACE  ".*'(\#| )([a-zA-Z0-9_]*)(:|'| ).*" "\\2" cur ${str})
      list(APPEND read_cur_use ${cur})
   endforeach()
ENDFOREACH()
list(REMOVE_DUPLICATES read_cur_use)
list(SORT read_cur_use)

#Read all parameters we are currently checking in read_check.f90, store them in @read_check_use
file(STRINGS forward/input/read_check.f90 check_file_content REGEX "key_char//")
foreach(str ${check_file_content})
   string(REGEX REPLACE  ".*'(\#| )([a-zA-Z0-9_]*)(:|'| ).*" "\\2" cur ${str})
   list(APPEND read_check_use ${cur})
endforeach()
list(REMOVE_DUPLICATES read_check_use)
list(SORT read_check_use)

#Now remove all the parameters we are currently checking from the @read_cur_use list
list(REMOVE_ITEM read_cur_use ${read_check_use})

# Do we need to update read_check.f90 with missing parameters in read_cur_use? 
if (read_cur_use)
   message("We have to update read_check.f90 with ${read_cur_use}")
   file(STRINGS forward/input/read_check.f90 check_file_content)
   set(automatic_section ${check_file_content})
   list(FILTER automatic_section INCLUDE REGEX ".*automatic generated.*")
   list(FIND check_file_content "${automatic_section}" automatic_element)
   MATH(EXPR automatic_element "${automatic_element}+1")
   foreach(item ${read_cur_use})
      list(INSERT check_file_content ${automatic_element} "       f_entry = f_entry + locstr(line,key_char\/\/' ${item}')")
      MATH(EXPR automatic_element "${automatic_element}+1")
   endforeach()
   file(WRITE forward/input/read_check.f90 "")
   foreach(str ${check_file_content})
      file(APPEND forward/input/read_check.f90 "${str}\n")
   endforeach()
endif()
#====================================================
