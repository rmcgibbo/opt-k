import sys
import string
import re
import urllib

msg = "This file was automatically generated from latex using pandoc\n\n"

sys.stdout.write(msg)

# match inline math or displayed math
math_re = re.compile('[$]+([^$]*)[$]+')

base_url = 'http://latex.codecogs.com/gif.latex'

f = open(sys.argv[1])
for line in f.readlines():
    m1 = math_re.search(line)
    if m1:
        # replace latex with a url that renders the math
        query = urllib.quote(m1.group(1))
        repl = '![equation](%s)' % (base_url + '?' + query)
        line =  math_re.sub(repl, line)
    
    sys.stdout.write(line)