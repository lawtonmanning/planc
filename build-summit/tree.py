import sys

if __name__ == '__main__':
    tree = sys.argv[1]
    dictionary = sys.argv[2]
    perm = sys.argv[3]
    
    nodes = {}
    with open(tree) as f:
        for line in f.readlines():
            items = [int(item) for item in line.split()]
            idx = items[0]
            nodes[idx] = items[1:]
    
    idxs = []
    with open(perm) as f:
        oidxs = [int(line) for line in f.readlines()]
        idxs = [0]*len(oidxs)
        for idx,val in enumerate(oidxs):
          idxs[val] = idx

    with open(dictionary) as f:
        words = f.read().splitlines()
        for idx in nodes:
            nodes[idx] = [words[idxs[i]] for i in nodes[idx]]

    for idx in nodes:
        print(idx, nodes[idx])


    
            
