BEGIN {RS = ">" ; FS = "\n" ; ORS = ""} $2 {print ">"$0}
