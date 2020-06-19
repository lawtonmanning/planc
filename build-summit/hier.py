from os import listdir
from os.path import isfile, join
import sys
from re import search, match
from csv import writer
from math import log2, floor

name = sys.argv[1]
totals = {};

count = 0

ms = []
ns = []
begin = True
end = False
nodes = []
nodeline = []
for line in open(name):
  total = match(r"total_(.+?)::.*min::(.+?)::avg::(.+?)::max::(.+?)$",line)
  if total:
    cat = total.group(1)
    times = [total.group(2), total.group(3), total.group(4)]
    if cat not in totals:
      totals[cat] = [times]
    else:
      totals[cat].append(times)
    begin = True
  
  size = match(r".*globalm::(\d+)::globaln::(\d+)",line)
  if size and begin:
    ms.append(size.group(1))
    ns.append(size.group(2))
    count += 1
    begin = False

  if end:
    nodes.append(line.split())

  idx = match(r"^idx*",line)
  if idx:
    nodeline = line.split()
    end = True



with open(name + '-nodes', 'w') as n, open(name + '-levels', 'w') as l:
  nwrite = writer(n,delimiter=',')
  lwrite = writer(l,delimiter=',')
  first = nodeline + ['m','n']
  for total in totals.keys():
    first.append(total + '-min')
    first.append(total + '-avg')
    first.append(total + '-max')
  nwrite.writerow(first)
  first[0] = "level"
  lwrite.writerow(first)
  
  levels = [floor(log2(int(node[0])+1)) for node in nodes]
  nodes_per_level = [levels.count(level) for level in list(set(levels))]
  full_levels = [count == 2**idx for idx,count in enumerate(nodes_per_level)]
  lrows = [[0.0]*len(first)]*(max(levels)+1)
  for idx in range(count):
    line = nodes[idx] + [ms[idx],ns[idx]]
    for total,times in totals.items():
      line = line + times[idx]
    nwrite.writerow(line)
    lrows[levels[idx]] = [x+float(y) for x,y in zip(lrows[levels[idx]],line)]
  
  for idx,row in enumerate(lrows):
    if full_levels[idx]:
      row[0] = idx
      lwrite.writerow(row)
       
  
