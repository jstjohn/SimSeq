
import sys
stdin = sys.stdin


for line in stdin:
    line = line.strip()
    if not line.startswith("#"):
        data = map(int,line.split()) #split the line and convert it to ints
        a=data[2:7]
        c=data[7:12]
        g=data[12:17]
        t=data[17:22]
        total=sum(a)
        oa=[x/float(total) for x in a]
        total=sum(c)
        oc=[x/float(total) for x in c]
        total=sum(g)
        og=[x/float(total) for x in g]
        total=sum(t)
        ot=[x/float(total) for x in t]
        res=map(str,data[:2])
        res+=map(str,oa)
        res+=map(str,oc)
        res+=map(str,og)
        res+=map(str,ot)
        print("\t".join(res))
    else:
        print(line)
