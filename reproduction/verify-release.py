"""Verify the files shipped in the publication snapshot (no downloads)."""
from pathlib import Path
import hashlib

ROOT = Path(__file__).resolve().parent.parent
manifest = ROOT / 'reproduction' / 'release-files.sha256'
failures = []
count = 0
for line in manifest.read_text().splitlines():
    expected, relative = line.split('  ', 1)
    path = ROOT / relative
    count += 1
    if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != expected:
        failures.append(relative)
if failures:
    raise SystemExit('Changed or missing release files:\n' + '\n'.join(failures))
print(f'Verified {count} publication files.')
