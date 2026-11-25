import pickle
D=pickle.load(open('big.pkl','rb')); F=pickle.load(open('family.pkl','rb'))
def sg(x): return (x>0)-(x<0)
laws={'S0':(D['S0'],'-++-',4),'S1':(D['S1'],'+--+',4),'mu2':(D['mu2'],'--++',4),
      'U1':(F['U1'],'++--',4),'U0':(F['U0'],'-++-+--+',8)}
for nm,(c,pat,M) in laws.items():
    bad=[n for n in range(len(c)) if sg(c[n])!={'+':1,'-':-1}[pat[n%M]]]
    print(nm,len(c)-1,'exceptions:',bad, [c[n] for n in bad])
