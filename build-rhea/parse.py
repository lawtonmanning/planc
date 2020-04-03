from os import listdir
from os.path import isfile, join
import sys
from re import search, match
from csv import writer

name = sys.argv[1]

lines = []
for line in open(name):
    info = match(r".*m::(\d+)::n::(\d+)::t::(\d+)::pr::(\d+)::pc::(\d+)",line)
    time = match(r".*NMF took (.*) secs",line)
    if info:
        m=info.group(1)
        n=info.group(2)
        t=info.group(3)
        pr=info.group(4)
        pc=info.group(5)
        p=int(pr)*int(pc)
    if time:
        lines.append([m,n,t,p,pr,pc,time.group(1)])


with open(name + '-output','w') as f:
    write = writer(f,delimiter=' ')
    write.writerow(['m','n','it','p','pr','pc','t'])
    for line in lines:
        write.writerow(line)