#!/usr/bin/env python

import sys
import optparse
import commands
import os
import glob

#def insertPlot(base, findex):
    
    
    
    
#######################
# Get options
#######################

parser = optparse.OptionParser("usage: %prog [options]\
<input directory> \n")
parser.add_option ('--o', dest='outdir', type='string',
                   default = './',
                   help="directory for output png")
parser.add_option ('--t', dest='title', type='string',
                   default = 'Tracking Plots',
                   help="Title of html page")
options, args = parser.parse_args()

if len(args) != 1:
    print "Please specify input dir.  Exiting."
    sys.exit()

indir  = args[0]+"/"
outdir = options.outdir+"/"
title = options.title

##############################################
# Check dir
##############################################
if not os.path.isdir(indir) :
    print "Cannot find %s.  Exiting." % infile
    sys.exit()

#if os.path.isfile(indir+"index.html") :
#    print "%s already exists.  Exiting." % indir+"index.html"

os.system("rm -f "+indir+"*small*.png")
files = glob.glob(indir+'*png')

for file in files:
    base = os.path.splitext(file)[0]
    print file, base
    command = "convert %s.png -resize 300x400 %s_small.png" % (base, base)
    os.system(command)


os.system("rm -f "+indir+"index.html")
os.system("touch "+indir+"index.html")
findex = open(indir+"index.html","w")

print >> findex, "<h1> "+title+"</h1>"
print >> findex, "<hr>\n"
print >> findex, "<table>"

nfile = 0
for file in files:
    base = os.path.basename(file).split(".")[0]
    print base

    if nfile == 0: print >> findex, "<tr>"
    elif nfile % 3 == 0: print >> findex, "</tr>\n<tr>"
    
    print >> findex, '<td><a href="%s.png"><img src="%s_small.png"></a></td>' % (base, base)  
    nfile += 1

print >> findex, "</tr>\n</table>\n"



findex.close() 
