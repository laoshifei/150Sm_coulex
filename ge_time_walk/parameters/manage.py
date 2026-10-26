#!/usr/bin/env python3
"""Verify or install the selected frozen timing package. Never changes active.txt.

No ROOT, fitting, raw event reading, or third-party Python packages are used.
Run from any directory: python3 /path/to/ge_time_walk/parameters/manage.py verify
"""
import argparse
import hashlib
import json
import re
import shutil
import tempfile
from pathlib import Path

BASE = Path(__file__).resolve().parent
PROJECT = BASE.parent.parent
REQUIRED = {
    'calibration.json', 'metadata.json', 'inherited_metadata.json',
    'correction_mode.txt', 'storage_format.txt', 'source_sha256.txt',
}


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def verify():
    selection = json.loads((BASE / 'selected_profile.json').read_text())
    profile = selection['profile']
    if not re.fullmatch('[a-z0-9_]+', profile):
        raise ValueError('Invalid selected profile name')
    source = BASE / profile
    package = json.loads((source / 'package.json').read_text())
    if set(package['files']) != REQUIRED:
        raise ValueError('Package must identify all six frozen source files')
    for name, digest in package['files'].items():
        if sha256(source / name) != digest:
            raise ValueError('Changed package file: ' + name)
    digest = package['files']['calibration.json']
    if selection['calibration_sha256'] != digest or (source / 'source_sha256.txt').read_text().strip() != digest:
        raise ValueError('Selection/source SHA does not match calibration')
    if (source / 'correction_mode.txt').read_text().strip() != 'total' or (source / 'storage_format.txt').read_text().strip() != 'json':
        raise ValueError('Expected a complete JSON profile')
    return source, selection, package


def install(source, selection, package, destination):
    """Install only into a new destination, or accept an already identical copy."""
    if destination.is_symlink():
        raise ValueError('Refuse symlink destination')
    if destination.exists():
        if not destination.is_dir() or set(p.name for p in destination.iterdir()) != REQUIRED | {'package.json'}:
            raise ValueError('Destination differs; choose a new version name')
        for name in REQUIRED | {'package.json'}:
            if sha256(destination / name) != sha256(source / name):
                raise ValueError('Destination differs: ' + name + '; choose a new version name')
        return 'already_identical'
    destination.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix='.' + destination.name + '.', dir=destination.parent))
    try:
        for name in sorted(REQUIRED | {'package.json'}):
            shutil.copy2(source / name, staging / name)
            if sha256(staging / name) != sha256(source / name):
                raise ValueError('Copy verification failed: ' + name)
        # One atomic publication; no file-by-file replacement of an old profile.
        staging.rename(destination)
    finally:
        if staging.exists():
            shutil.rmtree(staging)
    return 'installed'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=['verify', 'path', 'install'])
    parser.add_argument('--destination', type=Path, help='install directory; default comes from selected_profile.json')
    args = parser.parse_args()
    if args.destination is not None and args.command != 'install':
        parser.error('--destination is only valid for install')
    source, selection, package = verify()
    if args.command == 'path':
        print(source)
        return
    result = {'source': str(source), 'calibration_sha256': selection['calibration_sha256'],
              'status': 'verified', 'changes_active': False}
    if args.command == 'install':
        name = selection['runtime_version']
        if not re.fullmatch('[a-z0-9_]+', name):
            raise ValueError('Invalid runtime version name')
        destination = args.destination or PROJECT / 'pars/time_alignment/versions' / name
        destination = destination.absolute()
        result['status'] = install(source, selection, package, destination)
        result['destination'] = str(destination)
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
