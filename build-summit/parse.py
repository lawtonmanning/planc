from os import listdir
from os.path import isfile, join
import sys
from re import search, match
from csv import writer
from math import log2, floor

name = sys.argv[1]
totals = {}

for line in open(name):
  total = match(r"total_(.+?)::.*min::(.+?)::avg::(.+?)::max::(.+?)$",line)
  if total:
    cat = total.group(1)
    times = [total.group(2), total.group(3), total.group(4)]
    totals[cat] = times
  
  size = match(r".*globalm::(\d+?)::globaln::(\d+?)::.*",line)

  if size:
    m = size.group(1)
    n = size.group(2)
  
  procs = match(r".*pr::(\d+?)::pc::(\d+?)::.*",line)

  if procs:
    pr = int(procs.group(1))
    pc = int(procs.group(2))
    p = pr*pc


with open(name + '-output', 'w') as f:
  write = writer(f,delimiter=',')
  header = ['m','n','p','pr','pc']
  data = [m,n,p,pr,pc]
  for total in totals.keys():
    header.append(total + '-min')
    header.append(total + '-avg')
    header.append(total + '-max')
    data = data + totals[total]
  write.writerow(header)
  write.writerow(data)