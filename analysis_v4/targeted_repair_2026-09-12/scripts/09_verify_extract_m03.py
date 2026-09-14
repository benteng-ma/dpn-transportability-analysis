from pathlib import Path
import hashlib, json, shutil, tarfile

ROOT = Path(__file__).resolve().parents[1]
OLD=ROOT.parent/"all_extensions_2026-09-11"
arc=ROOT/"07_resource_recovery/downloads/GSE295206_RAW.tar"
out=ROOT/"07_resource_recovery/M03_selected_inputs"; out.mkdir(parents=True,exist_ok=True)
fail=OLD/"logs/GSE295206_SELECTED_DOWNLOAD_FAILURES.tsv"
wanted=[]
for i,line in enumerate(fail.read_text(encoding="utf-8").splitlines()):
    if i==0 or not line.strip(): continue
    wanted.append(line.split("\t")[0])
assert len(wanted)==48 and len(set(wanted))==48
h=hashlib.sha256()
with arc.open("rb") as f:
    for b in iter(lambda:f.read(16*1024*1024),b""): h.update(b)
members=[]; extracted=[]
with tarfile.open(arc,"r:") as tf:
    bybase={Path(m.name).name:m for m in tf.getmembers() if m.isfile()}
    members=[{"name":m.name,"size":m.size} for m in bybase.values()]
    for name in wanted:
        if name not in bybase: continue
        m=bybase[name]
        if Path(m.name).name!=name or ".." in Path(m.name).parts: raise RuntimeError("unsafe tar member")
        src=tf.extractfile(m)
        dest=out/name
        with src, dest.open("wb") as g: shutil.copyfileobj(src,g,16*1024*1024)
        extracted.append({"name":name,"expected_size":m.size,"actual_size":dest.stat().st_size,"sha256":hashlib.sha256(dest.read_bytes()).hexdigest()})
audit={"archive":arc.name,"archive_bytes":arc.stat().st_size,"archive_sha256":h.hexdigest(),"tar_member_files":len(members),"requested":len(wanted),"extracted":len(extracted),"all_sizes_match":all(x["expected_size"]==x["actual_size"] for x in extracted),"missing":sorted(set(wanted)-{x['name'] for x in extracted})}
(ROOT/"logs/09_verify_extract_m03.json").write_text(json.dumps(audit,indent=2),encoding="utf-8")
(ROOT/"07_resource_recovery/M03_TAR_MEMBERS.json").write_text(json.dumps(members,indent=2),encoding="utf-8")
with (ROOT/"07_resource_recovery/M03_SELECTED_FILE_INTEGRITY.tsv").open("w",encoding="utf-8",newline="") as f:
    f.write("name\texpected_size\tactual_size\tsha256\n")
    for x in extracted: f.write(f"{x['name']}\t{x['expected_size']}\t{x['actual_size']}\t{x['sha256']}\n")
if audit["archive_bytes"]!=9698846720 or len(extracted)!=48 or not audit["all_sizes_match"]: raise SystemExit(1)
print(json.dumps(audit,indent=2))
