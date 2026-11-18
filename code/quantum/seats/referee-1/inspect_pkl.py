import pickle, sys
for f in ['<LAB>/sideways/build-5/phase2/outward_C.pkl','<LAB>/engine/cert/lin_tube2_k1_3.pkl']:
    d = pickle.load(open(f,'rb'))
    print(f, type(d), list(d.keys()) if isinstance(d,dict) else '')
    for k in ['sm','cm']:
        if k in d:
            s = d[k]; print(' ',k, list(s.keys()))
            for kk in s:
                v = s[kk]
                if isinstance(v,list): print('   ',kk,'len',len(v), str(v[0])[:150] if v else '')
                else: print('   ',kk, str(v)[:150])
    for k in ['xc','V0','W0','params']:
        if k in d: print(' ',k, str(d[k])[:300])
