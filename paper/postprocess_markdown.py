import sys
import string

msg = "This file was automatically generated from latex using pandoc\n\n"

sys.stdout.write(msg)

f = open(sys.argv[1])
for line in f.readlines():
    # replace with backticks
    line = string.replace(line, '$$', '`')
    
    sys.stdout.write(line)