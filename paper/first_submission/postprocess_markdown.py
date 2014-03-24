import sys
import string
import re
import urllib
import itertools

msg = "This file was automatically generated from latex using pandoc\n\n"

sys.stdout.write(msg)

base_url = 'http://latex.codecogs.com/gif.latex'
github_url = "https://raw.github.com/rmcgibbo/opt-k/master/paper/"

# match inline math or displayed math
math_re = re.compile('[$]+([^$]*)[$]+')
fig_re = re.compile('\((figs/\S*\.png)')
label_re = re.compile(r'\\label{\S*}')
fig_location_re = re.compile('\[[ch!]+\]')
wierd_re = re.compile('\[fig:\S+\]')
caption_re = re.compile('!\[(.+)\]\(' + github_url)

f = open(sys.argv[1])
text = f.read()
counter = itertools.count()
while True:
    #print >> sys.stderr, counter.next()
    matched = False

    m0 = caption_re.search(text)
    if m0:
        matched = True
        # caption re matches the first part of the ![caption](img_link)
        # bit, and moves the caption out front, displaying it in italics
        # (the underlines)
        text, _ = caption_re.subn('_' + m0.group(1) + '_\n' + '![](' + github_url, text)

    m0 = fig_location_re.search(text)
    if m0:
        # just matches [h!] and removes it
        matched = True
        text, _ = fig_location_re.subn('', text)

    m0 = wierd_re.search(text)
    if m0:
        # removes "[fig:<random_junk>]" that gets inseted by pandoc?
        matched = True
        text, _ = wierd_re.subn('', text)

    m1 = fig_re.search(text)
    if m1:
        # remap the url in figures so that they point to the github
        # raw url for that img. Right now we're only detecting the img
        # links because they're in the fig/ directory and end in .png
        matched = True
        text, _ = fig_re.subn('(' + github_url + m1.group(1), text, 1)

    m2 = label_re.search(text)
    if m2:
        # remove the "\label{eq:<whatever>}" statement from inside
        # math equations
        matched = True
        text, _ = label_re.subn('', text)
    
    m3 = math_re.search(text)
    if m3:
        matched = True
        # replace latex with a url that renders the math
        # looks for $$<stuff>$$ or $<stuff>$
        query = urllib.quote(m3.group(1))
        repl = '![](%s)' % (base_url + '?' + query)
        text, _ =  math_re.subn(repl, text, 1)

    if not matched:
        break

sys.stdout.write(text)
