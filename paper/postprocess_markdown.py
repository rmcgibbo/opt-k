import sys
import string
import re
import urllib

msg = "This file was automatically generated from latex using pandoc\n\n"

sys.stdout.write(msg)

# match inline math or displayed math
math_re = re.compile('[$]+([^$]*)[$]+')
fig_re = re.compile('\((figs/\S*\.png)')
base_url = 'http://latex.codecogs.com/gif.latex'

github_url = "https://raw.github.com/rmcgibbo/opt-k/master/paper/"

f = open(sys.argv[1])
for line in f.readlines():
    m1 = fig_re.search(line)    
    if m1:
        print >> sys.stderr, 'inside m1'
        line = fig_re.sub(github_url + m1.group(1), line)

    while True:
        m1 = math_re.search(line)
        if m1:
            # replace latex with a url that renders the math
            query = urllib.quote(m1.group(1))
            repl = '![equation](%s)' % (base_url + '?' + query)
            line, _ =  math_re.subn(repl, line, 1)
        else:
            break
    
    sys.stdout.write(line)
