from __future__ import annotations

import os
import hashlib,json,zipfile
from datetime import datetime,timezone
from pathlib import Path

P = Path(os.environ.get("DPN_PROJECT_ROOT", Path(__file__).resolve().parents[3]))
C=P/"analysis_v4/closeout_20260912T062015Z"
ZIP=P/"analysis_v4/DPN_v4_CLOSEOUT_REVIEW_BUNDLE_2026-09-12_v1.zip"
SHA=Path(str(ZIP)+".sha256")
VERIFY=P/"analysis_v4/DPN_v4_CLOSEOUT_REVIEW_BUNDLE_2026-09-12_v1.verification.json"
MAN=C/"MANIFEST_SHA256.tsv"
def digest(p):
 h=hashlib.sha256()
 with p.open("rb") as f:
  for b in iter(lambda:f.read(1048576),b""):h.update(b)
 return h.hexdigest()
if ZIP.exists() or SHA.exists() or VERIFY.exists(): raise SystemExit("Refusing to overwrite an existing closeout archive")
pre={"sealed_utc":datetime.now(timezone.utc).isoformat(),"source":str(C),"zip":str(ZIP),"manifest_self_reference":False,"writers_confirmed_closed":True,"public_release":False}
(C/"09_verification/PACKAGE_PREFLIGHT.json").write_text(json.dumps(pre,indent=2),encoding="utf-8")
files=sorted(p for p in C.rglob("*") if p.is_file() and p!=MAN)
rows=[(p.relative_to(C).as_posix(),p.stat().st_size,digest(p)) for p in files]
with MAN.open("w",encoding="utf-8",newline="") as f:
 f.write("relative_path\tbytes\tsha256\n")
 for r,n,h in rows:f.write(f"{r}\t{n}\t{h}\n")
with zipfile.ZipFile(ZIP,"w",compression=zipfile.ZIP_DEFLATED,compresslevel=6,allowZip64=True) as z:
 for p in sorted([*files,MAN]):z.write(p,p.relative_to(C).as_posix())
zip_hash=digest(ZIP); SHA.write_text(f"{zip_hash}  {ZIP.name}\n",encoding="ascii")
with zipfile.ZipFile(ZIP,"r") as z:
 bad_crc=z.testzip(); names=set(z.namelist()); expected={r for r,_,_ in rows}|{"MANIFEST_SHA256.tsv"}
 hash_bad=[]
 for r,_,h in rows:
  got=hashlib.sha256(z.read(r)).hexdigest()
  if got!=h:hash_bad.append(r)
report={"status":"PASS" if bad_crc is None and names==expected and not hash_bad else "FAIL","zip":str(ZIP),"zip_bytes":ZIP.stat().st_size,"zip_sha256":zip_hash,"zip_members":len(names),"manifest_payloads":len(rows),"crc_bad_member":bad_crc,"member_set_exact":names==expected,"hash_mismatches":hash_bad,"manifest_lists_itself":any(r=="MANIFEST_SHA256.tsv" for r,_,_ in rows)}
VERIFY.write_text(json.dumps(report,indent=2),encoding="utf-8")
print(json.dumps(report,indent=2))
