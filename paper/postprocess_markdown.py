import sys
import string
import re
import urllib
import itertools

msg = "This file was automatically generated from latex using pandoc\n\n"

sys.stdout.write(msg)

# match inline math or displayed math
math_re = re.compile('[$]+([^$]*)[$]+')
fig_re = re.compile('\((figs/\S*\.png)')
label_re = re.compile(r'\\label{\S*}')
fig_location_re = re.compile('\[[ch!]+\]')
wierd_re = re.compile('\[fig:\S+\]')

base_url = 'http://latex.codecogs.com/gif.latex'

github_url = "https://raw.github.com/rmcgibbo/opt-k/master/paper/"

f = open(sys.argv[1])
text = f.read()
counter = itertools.count()
while True:
    print >> sys.stderr, counter.next()
    matched = False

    m0 = fig_location_re.search(text)
    if m0:
        matched = True
        text, _ = fig_location_re.subn('', text)

    m0 = wierd_re.search(text)
    if m0:
        matched = True
        text, _ = wierd_re.subn('', text)

    m1 = fig_re.search(text)
    if m1:
        matched = True
        text, _ = fig_re.subn('(' + github_url + m1.group(1), text, 1)

    m2 = label_re.search(text)
    if m2:
        matched = True
        text, _ = label_re.subn('', text)
    
    m3 = math_re.search(text)
    if m3:
        matched = True
        # replace latex with a url that renders the math
        query = urllib.quote(m3.group(1))
        repl = '![equation](%s)' % (base_url + '?' + query)
        text, _ =  math_re.subn(repl, text, 1)

    if not matched:
        break

sys.stdout.write(text)
