#!/usr/bin/env python3
"""Check the sizes of the actual source archive and its installed assets."""
import sys
import tarfile
from pathlib import Path

archive = Path(sys.argv[1])
with tarfile.open(archive) as source:
    files = [item for item in source.getmembers() if item.isfile()]
    docs = sum(item.size for item in files if "/inst/doc/" in item.name
               or "/man/" in item.name
               or item.name.split("/", 1)[-1] in {"README.md", "NEWS.md"})
    data = sum(item.size for item in files if "/inst/extdata/" in item.name
               or "/data/" in item.name)
sizes = {"source archive": (archive.stat().st_size, 10_000_000),
         "documentation": (docs, 5_000_000), "data": (data, 5_000_000)}
for name, (size, limit) in sizes.items():
    print(f"{name}: {size:,} bytes (limit {limit:,})")
if any(size > limit for size, limit in sizes.values()):
    raise SystemExit("Package size guard failed")
