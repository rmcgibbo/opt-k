import sys
import string
import re
import urllib

msg = "This file was automatically generated from latex using pandoc\n\n"

sys.stdout.write(msg)
inline_math = re.compile('[$]+([^$]*)[$]+')

base_url = 'http://latex.codecogs.com/gif.latex?'

f = open(sys.argv[1])
for line in f.readlines():
    # replace with backticks
    m1 = inline_math.search(line)
    if m1:
        #print >> sys.stderr, m1.group(1)
        query = urllib.quote(m1.group(1))
        line = '![equation](' + inline_math.sub(base_url + query, line) + ')'
        
    #line = string.replace(line, '$$', '`')
    #line = string.replace(line, '$', '`')
    
    #![equation](http://latex.codecogs.com/gif.latex?1%2Bsin%28mc%5E2%29%0D%0A)

    
    
    sys.stdout.write(line)