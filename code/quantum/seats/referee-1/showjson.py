import json,sys
d=json.load(open(sys.argv[1]))
def show(o,p=''):
  if isinstance(o,dict):
    for k,v in o.items(): show(v,p+'/'+k)
  elif isinstance(o,list):
    print(p,'list',len(o), str(o[:2])[:400])
  else: print(p,str(o)[:300])
show(d)
