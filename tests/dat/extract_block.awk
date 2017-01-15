#! /usr/bin/awk -f

BEGIN{
    nblocki = ARGV[2]
    nblocke = ARGV[3]
    if (nblocki == "" || nblocke == ""){
	print "missing initial or end block"
	exit
    }
    ARGC=2
}
/^# *DCP/{
    name = $4
}
$0 !~ /^ *\|/{
    inblock = ""
    next
}
/^ *\|/{
    if (inblock == ""){
	nblock++
	nl = 0
    }
    inblock = 1
    nl++ 
    if (nblock >= nblocki && nblock <= nblocke){
	if (nl > 1){
	    gsub(/f/,"",$6)
	    print
	} else {
	    printf("\n# DCP: %s\n",name)
	}
    }
}
