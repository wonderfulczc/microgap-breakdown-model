from pathlib import Path
import csv,sys,re
ROOT=Path(__file__).resolve().parents[1]; errors=[]; warnings=[]
def load(p):
 try:
  raw=p.read_bytes(); raw.decode('utf-8'); return list(csv.DictReader(raw.decode('utf-8').splitlines()))
 except Exception as e: errors.append((str(p),0,'',f'UTF-8/CSV read failure: {e}','save as valid UTF-8 CSV')); return []
S=load(ROOT/'docs/source_registry.csv'); P=load(ROOT/'docs/parameter_registry.csv'); sids={x['source_id'] for x in S}
def err(i,r,msg,fix): errors.append(('docs/parameter_registry.csv',i+2,r.get('parameter_id',''),msg,fix))
def warn(i,r,msg): warnings.append(('docs/parameter_registry.csv',i+2,r.get('parameter_id',''),msg))
ids=[r['parameter_id'] for r in P]
for x in set(ids):
 if ids.count(x)>1: errors.append(('docs/parameter_registry.csv',0,x,'duplicate parameter_id','renumber uniquely'))
pe={'direct','derived','inferred_from_figure','assumption','unresolved'}; ce={'high','medium','low','unknown'}; ie={'frozen','provisional','blocked','not_required_yet'}
for i,r in enumerate(P):
 if r['source_id'] and r['source_id'] not in sids: err(i,r,'source_id not in source registry','add source or correct ID')
 if r['provenance_type'] in {'direct','derived','inferred_from_figure'} and not r['source_location'].strip(): err(i,r,'missing source_location','add PDF page/equation/table/figure')
 if r['provenance_type']=='unresolved' and (r['value_original'].strip() or r['value_SI'].strip()): err(i,r,'unresolved row has definite value','clear values or change provenance')
 if bool(r['value_SI'].strip()) != bool(r['unit_SI'].strip()): err(i,r,'value_SI and unit_SI presence mismatch','fill or clear both')
 if r['confidence'] not in ce: err(i,r,'invalid confidence','use enumerated value')
 if r['provenance_type'] not in pe: err(i,r,'invalid provenance_type','use enumerated value')
 if r['implementation_status'] not in ie: err(i,r,'invalid implementation_status','use enumerated value')
 txt=' '.join(r.values()).lower()
 if r['provenance_type']=='direct' and any(x in txt for x in ['estimated','assumed','approximately inferred']): err(i,r,'direct row contains inference/assumption wording','change wording or provenance')
 if r['provenance_type']=='inferred_from_figure' and r['confidence']=='high' and 'reason' not in r['notes'].lower(): err(i,r,'high-confidence figure inference lacks reason','lower confidence or document reason')
 if r['canonical_name'].startswith('sp3_a') and r['unit_original']!='cm^-1 Torr^-1': err(i,r,'SP3 A_j unit family mismatch','use Table 3 cm^-1 Torr^-1')
 if r['canonical_name'].startswith('helmholtz_a') and r['unit_original']!='cm^-2 Torr^-2': err(i,r,'Helmholtz A_j unit family mismatch','use Table 2 cm^-2 Torr^-2')
 # stable scalar parsing where the field claims a simple numeric token
 for fld in ['value_original','value_SI']:
  v=r[fld].strip()
  if v and re.fullmatch(r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?',v):
   try: float(v)
   except: err(i,r,f'{fld} numeric parse failed','use stable decimal/scientific notation')
# conflicting simple numeric values per symbol/case
G={}
for i,r in enumerate(P):
 k=(r['symbol'],r['case_scope']); v=r['value_SI'].strip()
 if v and re.fullmatch(r'[+-]?[\d.]+(?:e[+-]?\d+)?',v,re.I): G.setdefault(k,set()).add(v)
for k,v in G.items():
 if len(v)>1: warnings.append(('docs/parameter_registry.csv',0,str(k),f'conflicting values: {sorted(v)}'))
for s in S:
 if not s['source_id'].startswith('SRC-'): errors.append(('docs/source_registry.csv',0,s['source_id'],'invalid source_id format','use SRC-0001'))
 if s['local_file'] and (not s['file_sha256'] or len(s['file_sha256'])!=64): errors.append(('docs/source_registry.csv',0,s['source_id'],'missing/invalid SHA256','recompute SHA256'))
for w in warnings: print('WARNING | file=%s | line=%s | id=%s | %s'%w)
for e in errors: print('ERROR | file=%s | line=%s | id=%s | reason=%s | fix=%s'%e)
print(f'Validated {len(S)} sources and {len(P)} parameters: {len(errors)} error(s), {len(warnings)} warning(s).')
sys.exit(1 if errors else 0)
