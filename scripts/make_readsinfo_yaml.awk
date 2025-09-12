#!/bin/awk -f 

BEGIN {
  FS = ","
}

NR == 1 {
	path = $0
	print "reads:"
}

NR > 1 {
	delete reads
	print "  " $1 ":"
	for (i = 2; i <= NF; i++) {
		if ($i ~ "=") {
			break
		}
		reads[i] = path "/" $i
	}
	j = i
	while (i++ < NF) {
		params[tolower($i)] = 1
	}
	if ("split=true" in params || "nosplit=false" in params) {
		print "    interleaved:"
	} else if ("format=paired" in params) {
		print "    paired:"
	} else {
		print "    single:"
	}
	for (i = 2; i < j; i++) {
		print "      - " reads[i]
	}
}

