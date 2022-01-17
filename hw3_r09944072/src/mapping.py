from contextlib import redirect_stdout
import os
import sys

def read(a): 
    d = {}
    with open(a, encoding = "big5hkscs") as f:
        content = f.readlines() 
    for i in range (len(content)):
        text = content[i].strip('\t')
        l = text.split() # l[0]: 國字、l[1]: 注音
        Zhuyin_list = l[1].split()
        
        for key in Zhuyin_list:
            if key[0] not in d:
                d[key[0]] = []
            d[key[0]].append(l[0])
    return d

def write (d):
    with open(sys.argv[2], 'w', encoding = "big5hkscs") as f:
        with redirect_stdout(f):
            for K, V in d.items():
                print (K + '\t' + ' '.join(V))
                for v in V:
                    print(v + '\t' + v)

d = read(sys.argv[1])
write(d)

