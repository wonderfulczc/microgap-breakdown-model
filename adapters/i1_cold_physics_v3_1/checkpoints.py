"""Accepted-only atomic envelope; verifies config, bytes and output prefix offsets."""
import hashlib,json,os,shutil
from pathlib import Path
VERSION='I1-ACCEPTED-CHECKPOINT-3.1'
class CheckpointRefusal(ValueError):pass
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def canonical(obj):return json.dumps(obj,sort_keys=True,separators=(',',':'),allow_nan=False).encode()
def identity(config):return hashlib.sha256(canonical(config)).hexdigest()
def syncfile(path,data):
    with path.open('xb') as f:f.write(data);f.flush();os.fsync(f.fileno())
def seal(payload,meta,output_dir,destination,config):
    destination=Path(destination);stage=destination.with_name(destination.name+'.incomplete')
    if destination.exists() or stage.exists():raise CheckpointRefusal('REFUSE_CHECKPOINT_OVERWRITE')
    stage.mkdir(parents=True)
    syncfile(stage/'payload.bin',Path(payload).read_bytes())
    outputs={}
    for name,offset in meta['output_offsets'].items():
        if Path(name).name!=name:raise CheckpointRefusal('UNSAFE_OUTPUT_NAME')
        data=(Path(output_dir)/name).read_bytes()
        if type(offset) is not int or not 0<=offset<=len(data):raise CheckpointRefusal('INVALID_OUTPUT_OFFSET')
        syncfile(stage/name,data[:offset]);outputs[name]={'offset':offset,'sha256':sha(stage/name)}
    manifest=dict(version=VERSION,kind='ACCEPTED_CHECKPOINT',identity=identity(config),config=config,
      accepted=meta,payload_sha256=sha(stage/'payload.bin'),output_prefixes=outputs,
      rebuild={'Poisson':'same immutable matrix/config; fresh RHS/zero guess','SP3':'same gas/metal matrix/config; zero guess; warm components preserved',
      'physical_history':'serialized; no density/time-only restart','scope':'same-build serial fixed prescribed voltage'},
      provider='ColdFullProvider',capacitance_mode=meta['capacitance_mode'])
    syncfile(stage/'manifest.json',canonical(manifest))
    syncfile(stage/'COMPLETE',sha(stage/'manifest.json').encode())
    fd=os.open(stage,os.O_RDONLY);os.fsync(fd);os.close(fd)
    stage.rename(destination)
    fd=os.open(destination.parent,os.O_RDONLY);os.fsync(fd);os.close(fd)
    return verify(destination,config)
def verify(path,config):
    p=Path(path)
    try:
        if p.is_symlink() or p.name.endswith('.incomplete'):raise CheckpointRefusal('INCOMPLETE_CHECKPOINT')
        if (p/'COMPLETE').read_text()!=sha(p/'manifest.json'):raise CheckpointRefusal('INCOMPLETE_OR_CORRUPT_MANIFEST')
        m=json.loads((p/'manifest.json').read_text())
        if m['version']!=VERSION or m['kind']!='ACCEPTED_CHECKPOINT':raise CheckpointRefusal('CHECKPOINT_VERSION_OR_KIND_MISMATCH')
        if m['identity']!=identity(config) or m['config']!=config:raise CheckpointRefusal('CHECKPOINT_CONFIG_IDENTITY_MISMATCH')
        if sha(p/'payload.bin')!=m['payload_sha256']:raise CheckpointRefusal('CHECKPOINT_PAYLOAD_SHA_MISMATCH')
        if m['capacitance_mode'] not in ('TOTAL_CHARGE','SPLIT_GEOMETRIC_CAPACITANCE'):raise CheckpointRefusal('CHECKPOINT_CAPACITANCE_MODE')
        if m['accepted']['capacitance_mode']!=m['capacitance_mode'] or config.get('capacitance_mode',m['capacitance_mode'])!=m['capacitance_mode']:
            raise CheckpointRefusal('CHECKPOINT_CAPACITANCE_IDENTITY_MISMATCH')
        if (p/'payload.bin').is_symlink():raise CheckpointRefusal('CHECKPOINT_PAYLOAD_SYMLINK')
        for name,item in m['output_prefixes'].items():
            if Path(name).name!=name or (p/name).is_symlink() or (p/name).stat().st_size!=item['offset'] or sha(p/name)!=item['sha256']:raise CheckpointRefusal('CHECKPOINT_OUTPUT_PREFIX_MISMATCH')
        return m
    except (OSError,KeyError,json.JSONDecodeError) as e:raise CheckpointRefusal('INCOMPLETE_CHECKPOINT') from e
def restore_outputs(path,config,out):
    m=verify(path,config);out=Path(out);out.mkdir(exist_ok=False,parents=True)
    for name in m['output_prefixes']:shutil.copyfile(Path(path)/name,out/name)
    return m
