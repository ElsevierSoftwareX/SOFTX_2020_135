#!/bin/csh

#comment out all calls in solve.f90
sed "s/^\([[:space:]]\+call[[:space:]]\+\)/\!\1/ig" solve.f90 | grep -v "newt\.inc" | grep -i -v "implicit[[:space:]]\+none" > tmp; mv tmp solve.f90

#commenting out calls to forward_newton, forward_nitsol, g_alloc_arrays, g_dealloc_arrays in forward_wrapper.f90
set file = "forward_wrapper.f90"; sed "s/^\([[:space:]]\+call[[:space:]]\+\)\(forward_newton\|forward_nitsol\)\b/\!\1\2/ig" $file | sed "s/^\([[:space:]]\+.*\+call[[:space:]]\+\)\(g_alloc_arrays\|g_dealloc_arrays\)\b/\!\1\2/ig" | sed -re "s/^( +stop.*)/\!\1/i" > tmp; mv tmp $file

#commenting out write-only calls
#commenting out continuations of commented out lines
foreach file ( `ls -1 *.f*` ) 
    #comment out write_* calls
    if ( `grep -i -c "^[[:space:]]\+call[[:space:]]\+write_" $file` > 0  ) then 
	sed -i "s/^\([[:space:]]\+call[[:space:]]\+\)\(write_\)/\!\1\2/ig" $file
    endif
    #coment out read_* calls
    if ( `grep -i -c "^[[:space:]]\+call[[:space:]]\+read_" $file` > 0  ) then 
	sed -i "s/^\([[:space:]]\+call[[:space:]]\+\)\(read_\)/\!\1\2/ig" $file
    endif
    
    #if there are multiline commands in the files
    if ( `grep -c !".*"'&'"[[:space:]]*"'$' $file` > 0 ) then
    #replace all ! comments with QQQQQQQQQ to find multiline comments
	sed "s/\!/QQQQQQQQQ/g" $file > tmp
	set problemline = `grep -h -n -A1 '^QQQQQQQQQ.*\&[[:space:]]*$' tmp | grep -v "^[[:digit:]]\+.QQQQQQQQQ" | grep -v "^--" | tail -1 | sed "s/^\([[:digit:]]\+\).*/\1/g"`

	while ( "$problemline" != "" )
	    echo "Multiline comment in $file (line ${problemline})."
	    sed "${problemline}s/^/QQQQQQQQQ/g" tmp | sed "s/QQQQQQQQQ/\!/g" > $file
	    sed "s/\!/QQQQQQQQQ/g" $file > tmp
	set problemline = `grep -h -n -A1 '^QQQQQQQQQ.*\&[[:space:]]*$' tmp | grep -v "^[[:digit:]]\+.QQQQQQQQQ" | grep -v "^--" | tail -1 | sed "s/^\([[:digit:]]\+\).*/\1/g"`
	end
   endif
end
rm -f tmp
