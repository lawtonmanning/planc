import sys

if __name__ == '__main__':
    tree = sys.argv[1]
    dictionary = sys.argv[2]
    
    nodes = {}
    with open(tree) as f:
        for line in f.readlines():
            items = [int(item) for item in line.split()]
            idx = items[0]
            nodes[idx] = items[1:]
    
    with open(dictionary) as f:
        words = f.read().splitlines()
        for idx in nodes:
            nodes[idx] = [words[i] for i in nodes[idx]]

    for node in nodes:
        print(node, nodes[node])

    
            