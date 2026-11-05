"""Fetch only the three paper datasets and verify original file SHA-256 hashes."""
from pathlib import Path
import argparse
import datetime
import hashlib
import io
import json
import urllib.request
import zipfile

BASE = Path(__file__).resolve().parent
MANIFEST = json.loads((BASE / 'datasets.json').read_text())
ROOT = BASE / 'engine' / 'applications' / 'data' / 'raw'


def check(spec):
    for name, expected in spec['files'].items():
        path = ROOT / spec['name'] / name
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            raise RuntimeError(f'Missing or changed dataset file: {path}. Run reproduction/fetch-data.py.')


def fetch(spec):
    target = ROOT / spec['name']
    try:
        check(spec)
        print(spec['name'], 'already verified')
        return
    except RuntimeError:
        pass
    target.mkdir(parents=True, exist_ok=True)
    print('Downloading', spec['name'], flush=True)
    request = urllib.request.Request(spec['url'], headers={'User-Agent': 'distkre-reproducibility/0.1'})
    with urllib.request.urlopen(request, timeout=120) as response:
        content = response.read()
    found = {}

    def members(blob):
        with zipfile.ZipFile(io.BytesIO(blob)) as archive:
            for member in archive.namelist():
                name = Path(member).name
                if name in spec['files']:
                    found[name] = archive.read(member)
                elif name.lower().endswith('.zip'):
                    members(archive.read(member))

    members(content)
    if set(found) != set(spec['files']):
        raise RuntimeError('Archive layout differs from the pinned dataset manifest.')
    for name, data in found.items():
        if hashlib.sha256(data).hexdigest() != spec['files'][name]:
            raise RuntimeError(f'Upstream data changed: {name}; no files replaced.')
    for name, data in found.items():
        path = target / name
        temporary = path.with_suffix(path.suffix + '.partial')
        temporary.write_bytes(data)
        temporary.replace(path)
    (target / 'download-provenance.json').write_text(json.dumps({
        'url': spec['url'], 'retrieved_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'archive_sha256': hashlib.sha256(content).hexdigest(), 'files': spec['files'],
    }, indent=2) + '\n')
    check(spec)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true', help='Verify existing files without downloading.')
    parser.add_argument('datasets', nargs='*', help='heart, beijing, spambase (default: all)')
    args = parser.parse_args()
    selected = args.datasets or [x['name'] for x in MANIFEST]
    if set(selected) - {x['name'] for x in MANIFEST}:
        parser.error('Dataset names must be heart, beijing, or spambase.')
    for spec in MANIFEST:
        if spec['name'] in selected:
            check(spec) if args.check else fetch(spec)
