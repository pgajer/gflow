#!/usr/bin/env python3
"""Check repository text and built distributions for non-public content."""

import argparse
import hashlib
import html
import json
import os
from pathlib import Path, PurePosixPath
import re
import subprocess
import sys
import tarfile
from datetime import datetime, timezone
from urllib.parse import unquote

VERSION = 1
MAX_TEXT_BYTES = 64 * 1024 * 1024
BINARY_SUFFIXES = {
    '.rds', '.rda', '.rdata', '.pdf', '.png', '.jpg', '.jpeg', '.gif', '.webp',
    '.ico', '.mp4', '.mov', '.zip', '.gz', '.bz2', '.xz', '.tar', '.o', '.so',
    '.dll', '.dylib', '.a', '.ttf', '.woff', '.woff2', '.eot', '.rdb', '.rdx',
}
RULES = {
    'private-reference': ('error', re.compile(
        r'(?:[.]codex[/\\]+(?:private|notes[/\\]+agent_instructions)(?:[/\\]|\b)'
        r'|/(?:Users|home)/[^/\s<>"\'|()\[\]{}]+(?:/|(?=[\s<>"\']|$))'
        r'|[A-Za-z]:[/\\]+Users[/\\]+[^/\\\s<>"\']+(?:[/\\]|(?=[\s<>"\']|$))'
        r'|file:/{3}[^\s<>"\']+'
        r'|/(?:private)/var/folders/)', re.I)),
    'internal-workflow': ('review', re.compile(
        r'\b(?:internal\s+(?:review|audit)\s+(?:notes?|status|checklist|handoff|summary|artifacts?|evidence|draft|only)'
        r'|(?:pending|awaiting)[-_\s]+(?:internal|agent|S3)[-_\s]+review'
        r'|review[- ]ready'
        r'|(?:agent|assistant)[- ](?:generated|created)\s+(?:report|artifact|handoff|notes?)'
        r'|(?:awaiting|pending)\s+(?:user|maintainer)\s+approval'
        r'|private\s+evidence\s+(?:is|was|lives|remains|stored)'
        r'|(?:CRAN|release)[- ]readiness\s+(?:report|assessment)'
        r'|reviewed\s+(?:functionality|experiment|means)'
        r'|implemented\s+in\s+Phase\s+[A-H]'
        r'|(?:next|another|other)\s+agent\s+(?:should|must|will))\b', re.I)),
    'prompt-remnant': ('error', re.compile(
        r'(?:\bBEGIN\s+(?:SYSTEM|USER|ASSISTANT)\s+PROMPT\b'
        r'|<\|(?:im_start|im_end|system|assistant|user)\|>'
        r'|\bAs\s+an?\s+(?:AI\s+(?:assistant|language\s+model)|large\s+language\s+model)\b'
        r'|\bYou\s+are\s+an?\s+expert\b'
        r'|\bYour\s+job\s+is\s+to\b'
        r'|\bI\s+will\s+(?:give|supply)\s+you\b'
        r'|\bIf\s+you\s+(?:want|like),[^\n]{0,180}\bI(?:[’\x27]ll|\s+will|\s+can)\b'
        r'|\[\[reply_to:[^\]]+\]\]'
        r'|^\s*#\s*\[(?:p|c)\]\s*$)', re.I | re.M)),
}
ARTIFACT = re.compile(
    r'(?:^|/)(?:[.]codex|[.]claude|cran-readiness(?:-[^/]*)?|agent-handoff[^/]*'
    r'|release-order[.]txt|citation_verification[.]html)(?:/|$)'
    r'|(?:^|/)[^/]*(?:chatgpt|claude|codex)[_-]instructions[^/]*(?:/|$)', re.I)


def normalize(text):
    # Decode escaped links without executing markup or following their targets.
    for _ in range(3):
        updated = html.unescape(unquote(text)).replace('\\/', '/')
        if updated == text:
            break
        text = updated
    return text.replace('\r\n', '\n').replace('\r', '\n')


def fingerprint(context):
    return hashlib.sha256(context.strip().encode('utf-8')).hexdigest()


def finding(scope, path, rule, level, line, context):
    return dict(scope=scope, path=path, rule=rule, level=level, line=line,
                fingerprint=fingerprint(context), context=context[:500])


def scan_text(scope, path, data):
    if Path(path).suffix.lower() in BINARY_SUFFIXES or b'\0' in data[:8192]:
        return [], False
    if len(data) > MAX_TEXT_BYTES:
        raise ValueError('Text file exceeds the scan limit: ' + path)
    text = normalize(data.decode('utf-8', errors='replace'))
    views = [text]
    if Path(path).suffix.lower() in {'.md', '.rmd', '.qmd', '.html', '.htm', '.xml', '.svg'}:
        visible = re.sub(r'<[^>]*>', lambda m: ' ' + '\n' * m[0].count('\n'), text)
        views.append(re.sub(r'[*`_]', '', visible))
    found, seen = [], set()
    lines = text.splitlines()
    for view in views:
        for rule, (level, pattern) in RULES.items():
            for match in pattern.finditer(view):
                line = view.count('\n', 0, match.start()) + 1
                last = view.count('\n', 0, match.end()) + 1
                context = '\n'.join(lines[line - 1:last]).strip()
                key = rule, line, context
                if key not in seen:
                    found.append(finding(scope, path, rule, level, line, context))
                    seen.add(key)
    return found, True


def repo_files(root):
    result = subprocess.run(
        ['git', 'ls-files', '--cached', '--others', '--exclude-standard', '-z'],
        cwd=root, capture_output=True, check=True)
    ignored = subprocess.run(
        ['git', 'ls-files', '--others', '--ignored', '--exclude-standard', '-z'],
        cwd=root, capture_output=True, check=True)
    # Ignored local documents and old source snapshots remain visible in a checkout.
    return sorted(set(os.fsdecode(p) for data in [result.stdout, ignored.stdout]
                      for p in data.split(b'\0') if p))


def source_inputs(root, tree=None):
    base = root if tree is None else tree
    if tree is None:
        paths = repo_files(root)
    else:
        if not tree.is_dir():
            raise ValueError('Rendered output directory does not exist: ' + str(tree))
        paths = sorted(p.relative_to(base).as_posix() for p in base.rglob('*')
                       if p.is_file() or p.is_symlink())
    for name in paths:
        p = base / name
        if not p.exists() and not p.is_symlink():
            continue  # Deleted tracked files are not part of the working source.
        if p.is_symlink():
            if p.readlink().is_absolute() or base not in p.resolve().parents or not p.is_file():
                raise ValueError('External, absolute or broken symbolic link: ' + name)
        if p.is_file():
            if p.suffix.lower() in BINARY_SUFFIXES:
                yield name, None
            else:
                with p.open('rb') as stream:
                    yield name, stream.read(MAX_TEXT_BYTES + 1)


def archive_inputs(archive, package):
    with tarfile.open(archive, 'r:gz') as source:
        for member in source:
            path = PurePosixPath(member.name)
            if path.is_absolute() or '..' in path.parts or not path.parts or path.parts[0] != package:
                raise ValueError('Unexpected archive member: ' + member.name)
            if member.isdir():
                continue
            if not member.isfile():
                raise ValueError('Non-regular archive member: ' + member.name)
            name = PurePosixPath(*path.parts[1:]).as_posix()
            if Path(name).suffix.lower() in BINARY_SUFFIXES:
                yield name, None
            else:
                with source.extractfile(member) as stream:
                    yield name, stream.read(MAX_TEXT_BYTES + 1)


def scan(scope, inputs):
    findings, counts = [], dict(files=0, text_files=0, binary_files=0)
    for name, data in inputs:
        counts['files'] += 1
        if ARTIFACT.search(name):
            findings.append(finding(scope, name, 'private-artifact', 'error', 0, name))
        matches, is_text = ([], False) if data is None else scan_text(scope, name, data)
        counts['text_files' if is_text else 'binary_files'] += 1
        findings.extend(matches)
    if counts['files'] == 0:
        raise ValueError('No files found for ' + scope)
    return findings, counts


def read_exceptions(root):
    config = root / 'tools/public-content-allowlist.json'
    if not config.exists():
        return []
    doc = json.loads(config.read_text())
    if doc.get('version') != VERSION or not isinstance(doc.get('exceptions'), list):
        raise ValueError('Invalid content-check exception file')
    entries, keys = [], set()
    for entry in doc['exceptions']:
        if set(entry) != {'scope', 'path', 'rule', 'fingerprint', 'reason'}:
            raise ValueError('Each exception requires scope, path, rule, fingerprint and reason')
        if not all(isinstance(v, str) for v in entry.values()):
            raise ValueError('Exception values must be strings')
        if entry['scope'] not in {'repo', 'archive', 'tree'} or entry['rule'] not in set(RULES) | {'private-artifact'}:
            raise ValueError('Unknown exception scope or rule')
        path = PurePosixPath(entry['path'])
        if path.is_absolute() or '..' in path.parts or any(c in entry['path'] for c in '*?[]'):
            raise ValueError('Exceptions require an exact relative path')
        if not re.fullmatch('[0-9a-f]{64}', entry['fingerprint']) or len(entry['reason'].strip()) < 12:
            raise ValueError('Exceptions require an exact fingerprint and a meaningful reason')
        key = tuple(entry[k] for k in ('scope', 'path', 'rule', 'fingerprint'))
        if key in keys:
            raise ValueError('Duplicate content-check exception')
        keys.add(key)
        entries.append(entry)
    return entries


def apply_exceptions(findings, entries, scopes):
    active, accepted, used = [], [], set()
    for match in findings:
        key = tuple(match[k] for k in ('scope', 'path', 'rule', 'fingerprint'))
        exception = next((e for e in entries if key == tuple(e[k] for k in
                         ('scope', 'path', 'rule', 'fingerprint'))), None)
        if exception:
            accepted.append(dict(match, reason=exception['reason']))
            used.add(key)
        else:
            active.append(match)
    stale = [e for e in entries if e['scope'] in scopes and tuple(e[k] for k in
             ('scope', 'path', 'rule', 'fingerprint')) not in used]
    return active, accepted, stale


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path, default=Path.cwd())
    target = parser.add_mutually_exclusive_group()
    target.add_argument('--archive', type=Path, help='Also inspect an exact source tarball')
    target.add_argument('--archive-dir', type=Path, help='Find the current package/version tarball in this directory')
    target.add_argument('--tree', type=Path, help='Also inspect a rendered site or documentation directory')
    parser.add_argument('--report-dir', type=Path, help='External report directory; also configurable with R_PACKAGE_CONTENT_REPORT_DIR')
    args = parser.parse_args(argv)
    root = args.repo.resolve()
    report, output = {}, None
    try:
        description = (root / 'DESCRIPTION').read_text()
        package = re.search(r'^Package:\s*(\S+)\s*$', description, re.M)[1]
        version = re.search(r'^Version:\s*(\S+)\s*$', description, re.M)[1]
        if not re.fullmatch(r'[A-Za-z][A-Za-z0-9.]*', package):
            raise ValueError('Invalid package name')
        default = Path.home() / '.codex' / 'private' / package / 'public-content-checks'
        directory = (args.report_dir or Path(os.environ.get('R_PACKAGE_CONTENT_REPORT_DIR', str(default)))).expanduser().resolve()
        if directory == root or root in directory.parents:
            raise ValueError('Reports must be written outside the package checkout')
        directory.mkdir(parents=True, exist_ok=True)
        stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S.%fZ')
        output = directory / (stamp + '.json')
        report = dict(checker_version=VERSION, package=package, version=version,
                      checked_utc=stamp, scopes={}, findings=[], exceptions=[], errors=[])
        findings, counts = scan('repo', source_inputs(root))
        report['scopes']['repo'] = counts
        archive = args.archive
        if args.archive_dir:
            archive = args.archive_dir / (package + '_' + version + '.tar.gz')
        if archive:
            matches, counts = scan('archive', archive_inputs(archive, package))
            findings.extend(matches)
            report['scopes']['archive'] = counts
            report['archive'] = dict(path=str(archive.resolve()),
                                     sha256=hashlib.sha256(archive.read_bytes()).hexdigest())
        if args.tree:
            matches, counts = scan('tree', source_inputs(root, args.tree.resolve()))
            findings.extend(matches)
            report['scopes']['tree'] = counts
        active, accepted, stale = apply_exceptions(findings, read_exceptions(root), report['scopes'])
        report.update(findings=active, exceptions=accepted, stale_exceptions=stale)
        report['status'] = 'FAIL' if active or stale else 'PASS'
        output.write_text(json.dumps(report, indent=2) + '\n')
        print(f"Public content: {report['status']} ({package}; {len(active)} findings, {len(stale)} stale exceptions)")
        for match in active[:50]:
            print(f"  {match['scope']}:{match['path']}:{match['line']}: {match['rule']} ({match['level']})")
        if len(active) > 50:
            print('  Additional findings are recorded in the external report.')
        for entry in stale:
            print(f"  {entry['scope']}:{entry['path']}: stale exception for {entry['rule']}")
        print('Report:', output)
        return 1 if active or stale else 0
    except (OSError, ValueError, TypeError, subprocess.CalledProcessError, tarfile.TarError) as error:
        if output:
            report.update(status='ERROR', errors=[str(error)])
            output.write_text(json.dumps(report, indent=2) + '\n')
        print('Public content check could not complete: ' + str(error), file=sys.stderr)
        return 2


if __name__ == '__main__':
    sys.exit(main())
