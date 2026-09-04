"""Copy members out of a concatenated split zip whose local headers sit past
the first slice. unzip assumes single-disk offsets and fails on them."""
import os, re, struct, sys, zipfile, zlib

# run against the archive folder; the joined archive is the three delivered
# slices concatenated in order, cat *.z01 *.z02 *.zip > merged_008271.zip
os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "02.inputs", "SAR"))
Z = "merged_008271.zip"
parts = [os.path.getsize(f"BDP_LakeChilwa_esar-ds_ASA_APP_1P_008271.{e}") for e in ("z01", "z02")]
# the reader adds the size of the earlier slices to every offset as if all
# entries sat on the last slice, so both signs are tried and the local header
# signature plus the filename decides
disk_starts = sorted({0, parts[0], parts[0] + parts[1], -parts[0], -(parts[0] + parts[1])})
man = {m.group(1): int(m.group(2)) for m in
       (re.match(r"(ASA_APP\S+\.N1)\s+(\d+) Bytes", l.strip()) for l in open("BDP_LakeChilwa_esar-ds_manifest_008271.txt"))
       if m}
zf = zipfile.ZipFile(Z)
fh = open(Z, "rb")
for info in zf.infolist():
    n = info.filename
    if n not in man or (os.path.exists(n) and os.path.getsize(n) == man[n]):
        continue
    found = None
    for ds in disk_starts:
        fh.seek(info.header_offset + ds)
        h = fh.read(30)
        if h[:4] != b"PK\x03\x04":
            continue
        nlen, xlen = struct.unpack("<HH", h[26:30])
        if fh.read(nlen).decode() == n:
            found = info.header_offset + ds + 30 + nlen + xlen
            break
    if found is None:
        print("NOT FOUND", n, flush=True); continue
    fh.seek(found)
    remaining = info.compress_size
    dec = zlib.decompressobj(-15) if info.compress_type == zipfile.ZIP_DEFLATED else None
    with open(n, "wb") as out:
        while remaining > 0:
            chunk = fh.read(min(64 << 20, remaining)); remaining -= len(chunk)
            out.write(dec.decompress(chunk) if dec else chunk)
        if dec: out.write(dec.flush())
    ok = os.path.getsize(n) == man[n]
    print(("ok  " if ok else "BAD ") + n, os.path.getsize(n), flush=True)
