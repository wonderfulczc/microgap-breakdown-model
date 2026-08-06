from pathlib import Path
import csv,yaml
REQUIRED=('alpha','beta','I_CM0','T0','source_id','source_location')
def load_config(path):
    p=Path(path); cfg=yaml.safe_load(p.read_text(encoding='utf-8'))
    for name in ('primary','figure3_caption_variant','slow_decay','later_transition'):
        if name not in cfg or any(k not in cfg[name] for k in REQUIRED): raise ValueError(f'incomplete config case {name}')
    return cfg
def params(case): return {'alpha':float(case['alpha']),'beta':float(case['beta']),'I0':float(case['I_CM0']),'T0':float(case['T0'])}
def cross_check_registry(cfg,registry):
    rows=list(csv.DictReader(Path(registry).open(encoding='utf-8'))); p=cfg['primary']; checks={'growth_rate':p['alpha'],'decay_rate':p['beta'],'current_moment_scale':p['I_CM0'],'peak_time':p['T0']}
    for suffix,val in checks.items():
        row=next((r for r in rows if r['canonical_name']==suffix),None)
        if row is None or abs(float(row['value_SI'])-float(val))>1e-12*max(1,abs(float(val))): raise ValueError(f'registry mismatch for {suffix}')

