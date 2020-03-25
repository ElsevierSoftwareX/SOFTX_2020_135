#!/bin/csh
set file0 = "treelist0.txt"
set file1 = "filelist.txt"
set file2 = "treelist.txt"
set file3 = "treelist.tmp"
set file4 = "treelist.tmp2"
set file5 = "treelevels.txt"
set file6 = "filelist_notneeded.txt"
set base_file = "forward_compute.f90"
#set solve_file = "solve.f90"

grep -i -H "^[[:space:]]*subroutine" *.f* |  sed "s/^\(..*\)\:.*subroutine\s\s*\(\w\w*\).*/\1 \2/ig" | sort | uniq > $file1
grep -i -H "^[[:space:]]\+[^[:punct:]]*function" *.f* |  sed "s/^\(..*\)\:.*function\s\s*\(\w\w*\).*/\1 \2/ig" | sort | uniq >> $file1
grep -i -H "^[[:space:]]*module" *.f* |  sed "s/^\(..*\)\:.*module\s\s*\(\w\w*\).*/\1 \2/ig" | sort | uniq >> $file1

echo $base_file > $file3
cp -f $file3 $file2
touch $file5
rm -f $file5
set iter = 1

echo "Constructing tree of required subroutines, functions, and modules."
echo -n "   Iteration step " 
while ( `more $file3 | wc -l ` > 0 && $iter < 150 )
  cp -f $file2 $file4
  less $file3 | sed "s/^/$iter /g" >> $file5
  echo -n "$iter "
  set n3 = `more $file3 | wc -l `
  set i3 = 1

  while ( $i3 <= $n3 ) 
    set file = `more $file3 | head -$i3 | tail -1`
#    if ( "$file" == "$solve_file" ) then
#	set file = "$base_file"
#    endif

    set calllist = `grep -i -h "^[[:space:]]\+\bcall\b" $file | sed "s/^\s\+call \+\(\w\+\)\b.*/\1/ig" | sort | uniq `
    foreach call ( $calllist)
#      echo "call" $file $call
      grep -i "[[:space:]]\+\b$call\b" $file1 | cut -d " " -f 1 >> $file4
    end

    set extrnlist = ( `grep -i -h "^[[:space:]]\+.*\bexternal\b" $file | sed "s/^\s\+external \+\([[:alnum:][:space:],_]\+\).*/\1/ig" | sed "s/ *, *\([[:alnum:]_]\)/\n\1/g" | sed "s/,//g" | sort | uniq ` `grep -i -h -A1 "^[[:space:]]\+\bexternal\b.*\&" $file | grep -v -i "^[[:space:]]\+\bexternal\b" | grep -v -i "^--" | sed "s/ *//g" | sed "s/,/\n/g" | sort | uniq  ` )
    foreach extrn ( $extrnlist)
#      echo "external" $file $extrn
      grep -i "[[:space:]]\+\b$extrn\b" $file1 | cut -d " " -f 1 >> $file4
    end

    set modllist = `grep -i -h "^[[:space:]]\+\buse\b" $file | sed "s/^\s\+use \+\(\w\+\)\b.*/\1/ig" | sort | uniq `
    foreach modl ( $modllist)
      grep -i "[[:space:]]\+\b$modl\b" $file1 | cut -d " " -f 1 >> $file4
    end

    grep -i -h "^[[:space:]\#]\+\binclude\b" $file | sed "s/^[[:space:]\#]\+include \+.\([[:alnum:]_\.]\+\)\b.*/\1/ig" | sort | uniq  >> $file4

    @ i3 ++
  end

  more $file4 | sort | uniq > $file3
  diff $file2 $file3 | grep -i "^> " | cut -d " " -f 2 > $file4
  cp -f $file3 $file2
  cp -f $file4 $file3
  @ iter ++
end

ls -1 *.f* *.inc | sort > $file0
diff $file2 $file0 | grep -i "^> " | cut -d " " -f 2 | sed "s/^/rm -f /g" > $file6

rm -f $file3 $file4
echo " --- Done."
