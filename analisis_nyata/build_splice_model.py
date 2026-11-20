"""Ekstrak model kompetensi splicing dari Source Data Shen dkk. 2025 (Nat Commun 16:7376)."""
import openpyxl, json, statistics as st
wb=openpyxl.load_workbook('/tmp/xl/650787bd-96e3-4fd7-8199-ec8584de30fd.xlsx',read_only=True,data_only=True)
def rows(name): 
    return [list(r) for r in wb[name].iter_rows(values_only=True)]
def mean(v): 
    v=[x for x in v if isinstance(x,(int,float))]; return st.mean(v) if v else None
out={}
def table(title,r,labcols):
    i=[k for k,x in enumerate(r) if x and str(x[0]).strip()==title][0]
    res=[]
    for x in r[i+2:]:
        if x[0] is None or isinstance(x[0],str) and x[0].startswith('Figure'): break
        res.append(x)
    return res
r=rows('Figure 2')
def sect(label):
    i=[k for k,x in enumerate(r) if x[0]==label][0]; res=[]
    for x in r[i+2:]:
        if x[0] is None or (isinstance(x[0],str) and x[0].startswith('Figure')): break
        res.append(x)
    return res
# 2b: sym exons, cols: a,b, 37C x3 (cols 2-4), 55C x3 (5-7)
out['fig2b']={int(x[0]):mean(x[5:8]) for x in sect('Figure 2b')}
out['fig2c_5exon']={int(x[0]):mean(x[5:8]) for x in sect('Figure 2c')}
out['fig2d_3exon']={int(x[1]):mean(x[5:8]) for x in sect('Figure 2d')}
out['fig2b_37']={int(x[0]):mean(x[2:5]) for x in sect('Figure 2b')}
# 2e arms: label, 37 x3, 55 x3 -> cols1-3,4-6 ; some rows have 5 values; use last 3 non-None as 55C if 6 else
def arm(x):
    v=[c for c in x[1:] if c is not None]
    return v
out['fig2e_arm']={}
for x in sect('Figure 2e'):
    v=arm(x); out['fig2e_arm'][int(x[0])]=mean(v[-2:])
out['fig2f_arm_MFI']={}
for x in sect('Figure 2f'):
    if x[0]=='Ctrl': continue
    out['fig2f_arm_MFI'][int(x[0])]=mean([c for c in x[1:] if c is not None])
sp={}
for x in sect('Figure 2g'):
    v=[c for c in x[1:] if c is not None]; sp[x[0]]=mean(v[-2:])
out['fig2g_spacer']=sp
r6=rows('Figure 6')
def sect6(label):
    i=[k for k,x in enumerate(r6) if x[0]==label][0]; res=[]
    for x in r6[i+2:]:
        if x[0] is None or (isinstance(x[0],str) and x[0].startswith('Figure')): break
        res.append(x)
    return res
out['fig6b_scarless']={f"{int(x[0])}/{int(x[1])}":x[2] for x in sect6('Figure 6b')}
out['fig6d_5exon']={x[2]:mean(x[3:6]) for x in sect6('Figure 6d')}
out['fig6f_split']={}
for x in sect6('Figure 6f'):
    out['fig6f_split'][int(x[0])]={'exon5':str(x[1]).split()[0],'eff':mean(x[2:5])}
json.dump(out,open('data/real/splice_model_data.json','w'),indent=1)
for k,v in out.items(): print(k,{a:(round(b,1) if isinstance(b,float) else b) for a,b in (v.items())})
