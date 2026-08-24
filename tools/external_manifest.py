#!/usr/bin/env python3

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
EXTERNAL_ROOT = REPO_ROOT / 'external'
MANIFEST_PATH = REPO_ROOT / 'external-manifest.json'

CHUNK_SIZE = 1024 * 1024
VCS_CONTROL_NAMES = frozenset({'.git', '.hg', '.svn'})


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open('rb') as file:
        while chunk := file.read(CHUNK_SIZE):
            digest.update(chunk)

    return digest.hexdigest()


def external_files() -> list[Path]:
    if not EXTERNAL_ROOT.is_dir():
        raise RuntimeError(f'missing external directory: {EXTERNAL_ROOT}')

    files: list[Path] = []

    for dirpath, dirnames, filenames in os.walk(
        EXTERNAL_ROOT,
        followlinks=False,
    ):
        directory = Path(dirpath)

        # Submodule/VCS metadata is not vendored source and is not sealed.
        dirnames[:] = [
            dirname
            for dirname in dirnames
            if dirname not in VCS_CONTROL_NAMES
        ]

        for dirname in dirnames:
            path = directory / dirname
            if path.is_symlink():
                raise RuntimeError(
                    f'symlinked directory not allowed: {path}'
                )

        for filename in filenames:
            if filename in VCS_CONTROL_NAMES:
                continue

            path = directory / filename
            if path.is_symlink():
                raise RuntimeError(
                    f'symlinked file not allowed: {path}'
                )

            files.append(path)

    return sorted(
        files,
        key=lambda path: path.relative_to(EXTERNAL_ROOT).as_posix(),
    )


def snapshot() -> dict[str, dict[str, object]]:
    result: dict[str, dict[str, object]] = {}

    for path in external_files():
        relative = path.relative_to(EXTERNAL_ROOT).as_posix()
        result[relative] = {
            'sha256': sha256_file(path),
            'size': path.stat().st_size,
        }

    return result


def create_manifest(force: bool) -> int:
    if MANIFEST_PATH.exists() and not force:
        print(
            f'{MANIFEST_PATH.name} already exists.\n'
            'Refusing to overwrite it without --force.',
            file=sys.stderr,
        )
        return 1

    manifest = {
        'algorithm': 'sha256',
        'root': 'external',
        'schema': 1,
        'files': snapshot(),
    }

    temporary = MANIFEST_PATH.with_suffix('.json.tmp')

    with temporary.open('w', encoding='utf-8', newline='\n') as file:
        json.dump(
            manifest,
            file,
            indent=2,
            sort_keys=True,
        )
        file.write('\n')

    os.replace(temporary, MANIFEST_PATH)

    file_count = len(manifest['files'])
    print(
        f'Wrote {file_count} files to '
        f'{MANIFEST_PATH.name}.'
    )
    return 0


def load_manifest() -> dict[str, object]:
    if not MANIFEST_PATH.is_file():
        raise RuntimeError(f'missing manifest: {MANIFEST_PATH}')

    try:
        with MANIFEST_PATH.open('r', encoding='utf-8') as file:
            manifest = json.load(file)
    except (OSError, json.JSONDecodeError) as error:
        raise RuntimeError(f'unable to read manifest: {error}') from error

    if not isinstance(manifest, dict):
        raise RuntimeError('manifest root must be a JSON object')
    if manifest.get('schema') != 1:
        raise RuntimeError('unsupported manifest schema')
    if manifest.get('algorithm') != 'sha256':
        raise RuntimeError('unsupported manifest algorithm')
    if manifest.get('root') != 'external':
        raise RuntimeError('manifest root must be \'external\'')
    if not isinstance(manifest.get('files'), dict):
        raise RuntimeError('manifest files entry must be a JSON object')

    return manifest


def verify_manifest(stamp: str | None) -> int:
    manifest = load_manifest()
    expected = manifest['files']
    actual = snapshot()

    expected_paths = set(expected)
    actual_paths = set(actual)

    missing = sorted(expected_paths - actual_paths)
    added = sorted(actual_paths - expected_paths)
    modified = sorted(
        path
        for path in expected_paths & actual_paths
        if expected[path].get('sha256') != actual[path]['sha256']
        or expected[path].get('size') != actual[path]['size']
    )

    for path in missing:
        print(f'MISSING   external/{path}')
    for path in added:
        print(f'ADDED     external/{path}')
    for path in modified:
        print(f'MODIFIED  external/{path}')

    if missing or added or modified:
        print('\nExternal dependency integrity check FAILED.')
        return 1

    print(f'External dependency integrity OK ({len(actual)} files).')

    if stamp:
        stamp_path = Path(stamp)
        stamp_path.parent.mkdir(parents=True, exist_ok=True)
        stamp_path.write_text(
            'external dependency integrity verified\n',
            encoding='utf-8',
        )

    return 0


def main() -> int:
    parser = argparse.ArgumentParser(
        description='Create or verify the sealed external dependency manifest.'
    )
    subparsers = parser.add_subparsers(dest='command', required=True)

    create = subparsers.add_parser('create')
    create.add_argument(
        '--force',
        action='store_true',
        help='replace an existing manifest',
    )

    verify = subparsers.add_parser('verify')
    verify.add_argument(
        '--stamp',
        help='write a stamp file after successful verification',
    )

    args = parser.parse_args()

    try:
        if args.command == 'create':
            return create_manifest(args.force)
        if args.command == 'verify':
            return verify_manifest(args.stamp)
    except (OSError, RuntimeError) as error:
        print(
            f'External dependency integrity check failed: {error}',
            file=sys.stderr,
        )
        return 1

    return 2


if __name__ == '__main__':
    raise SystemExit(main())