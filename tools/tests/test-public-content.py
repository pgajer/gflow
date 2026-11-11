"""Exercise release guards with isolated repositories and built-document fixtures."""
import io
import json
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile
import unittest

CHECKER = Path(__file__).resolve().parents[1] / 'check-public-content.py'
PRIVATE = '~/' + '.codex/' + 'private/example/measurements.txt'
WORKFLOW = 'Internal ' + 'review notes: awaiting sign-off.'


class ContentCheckTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.base = Path(self.temp.name)
        self.root = self.base / 'package'
        self.root.mkdir()
        self.reports = self.base / 'reports'
        self.write('DESCRIPTION', 'Package: fixture\nVersion: 1.0.0\n')
        subprocess.run(['git', 'init', '-q'], cwd=self.root, check=True)

    def write(self, name, text):
        path = self.root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
        return path

    def run_check(self, *args, status=0):
        result = subprocess.run(
            [sys.executable, str(CHECKER), '--repo', str(self.root),
             '--report-dir', str(self.reports), *map(str, args)],
            cwd=self.base, capture_output=True, text=True)
        self.assertEqual(result.returncode, status, result.stdout + result.stderr)
        reports = sorted(self.reports.glob('*.json'))
        return json.loads(reports[-1].read_text()) if reports else None

    def archive(self, files, name='fixture_1.0.0.tar.gz'):
        path = self.base / name
        with tarfile.open(path, 'w:gz') as tar:
            for member, text in files.items():
                data = text.encode()
                entry = tarfile.TarInfo(member)
                entry.size = len(data)
                tar.addfile(entry, io.BytesIO(data))
        return path

    def test_scientific_language_is_not_a_keyword_ban(self):
        self.write('methods.md', 'We review the literature and audit numerical stability.\n'
                   'Benchmark results, reviewer suggestions and regression tests are documented.\n'
                   'This limitation remains unresolved. See https://example.org/methods.\n')
        report = self.run_check()
        self.assertEqual(report['status'], 'PASS')

    def test_untracked_private_link_fails(self):
        self.write('README.md', '[Measurements](' + PRIVATE + ')\n')
        report = self.run_check(status=1)
        self.assertEqual(report['findings'][0]['path'], 'README.md')
        self.assertEqual(report['findings'][0]['line'], 1)

    def test_rbuildignore_does_not_hide_public_source(self):
        self.write('.Rbuildignore', '^notes$\n')
        self.write('notes/methods.tex', '\\href{' + PRIVATE + '}{Evidence}\n')
        self.assertEqual(self.run_check(status=1)['findings'][0]['path'], 'notes/methods.tex')

    def test_tracked_gitignored_files_are_checked(self):
        self.write('.gitignore', 'notes/\n')
        self.write('notes/results.txt', PRIVATE)
        subprocess.run(['git', 'add', '-f', 'notes/results.txt'], cwd=self.root, check=True)
        self.run_check(status=1)

    def test_ignored_builds_are_checked(self):
        self.write('.gitignore', 'site/\n')
        self.write('site/index.html', '<!-- ' + PRIVATE + ' -->')
        self.assertEqual(self.run_check(status=1)['findings'][0]['scope'], 'repo')
        report = self.run_check('--tree', self.root / 'site', status=1)
        self.assertEqual({x['scope'] for x in report['findings']}, {'repo', 'tree'})

    def test_ignored_development_prompt_is_checked(self):
        self.write('.gitignore', 'dev/\n')
        self.write('dev/rchk/instructions.md', 'You are ' + 'an expert C++ engineer.\n')
        report = self.run_check(status=1)
        self.assertEqual(report['findings'][0]['rule'], 'prompt-remnant')

    def test_conversation_markers_and_review_claims(self):
        self.write('notes.md', '# [' + 'p]\nA question\n# [' + 'c]\nAn answer\n')
        self.write('R/messages.R', 'message("The reviewed ' + 'functionality is available.")\n')
        report = self.run_check(status=1)
        self.assertEqual({x['rule'] for x in report['findings']},
                         {'prompt-remnant', 'internal-workflow'})

    def test_planning_status_in_ownership_table(self):
        self.write('api-ownership.csv', 'name,status\nmethod,KEEP-PENDING-' + 'S3-REVIEW\n')
        self.assertEqual(self.run_check(status=1)['findings'][0]['rule'], 'internal-workflow')

    def test_reply_marker_and_prompt_filename(self):
        self.write('notes.md', '[[' + 'reply_to:123]]\n')
        self.write('dev/rchk_' + 'chatGPT_instructions.md', 'Correction rules\n')
        report = self.run_check(status=1)
        self.assertEqual({x['rule'] for x in report['findings']},
                         {'prompt-remnant', 'private-artifact'})

    def test_conversational_offer_in_developer_document(self):
        self.write('tools/headers.md', 'If ' + 'you want, paste a header and I will rewrite it.\n')
        report = self.run_check(status=1)
        self.assertEqual(report['findings'][0]['rule'], 'prompt-remnant')

    def test_escaped_html_link_and_hidden_comment(self):
        encoded = PRIVATE.replace('/', '%252F').replace('.', '&#46;')
        self.write('guide.html', '<!-- <a href="' + encoded + '">data</a> -->')
        self.run_check(status=1)

    def test_windows_and_macos_paths(self):
        self.write('provenance.txt', 'C:\\Users\\someone\\Desktop\\input.csv\n'
                   + '/Users/' + 'someone/data.csv\n')
        report = self.run_check(status=1)
        self.assertEqual(len(report['findings']), 2)

    def test_home_directory_and_file_uri(self):
        self.write('methods.md', '"' + '/Users' + '/someone"\n'
                   + '[Evidence](' + 'file:' + '/' * 3 + 'tmp/report.html)')
        self.assertEqual(len(self.run_check(status=1)['findings']), 2)

    def test_wrapped_coordination_wording(self):
        self.write('methods.md', 'Internal\nreview notes: ready for the next pass.\n')
        report = self.run_check(status=1)
        self.assertEqual(report['findings'][0]['rule'], 'internal-workflow')
        self.assertEqual(report['findings'][0]['line'], 1)

    def test_prompt_marker(self):
        self.write('notes.txt', 'BEGIN ' + 'SYSTEM PROMPT\n')
        self.assertEqual(self.run_check(status=1)['findings'][0]['rule'], 'prompt-remnant')

    def test_inline_markup_does_not_hide_coordination_wording(self):
        self.write('guide.html', '<p>Internal <em>review</em> notes: draft.</p>')
        self.write('README.md', 'Internal **review** notes: draft.')
        self.assertEqual(len(self.run_check(status=1)['findings']), 2)

    def test_internal_artifact_filename(self):
        self.write('release-' + 'order.txt', 'Package dependency order\n')
        self.assertEqual(self.run_check(status=1)['findings'][0]['rule'], 'private-artifact')

    def test_generated_archive_content_is_checked(self):
        archive = self.archive({'fixture/inst/doc/guide.html': '<p>' + WORKFLOW + '</p>'})
        report = self.run_check('--archive', archive, status=1)
        self.assertEqual(report['findings'][0]['scope'], 'archive')
        self.assertEqual(len(report['archive']['sha256']), 64)

    def test_current_archive_directory_and_missing_archive(self):
        self.archive({'fixture/README.md': 'Published methods.'})
        self.run_check('--archive-dir', self.base)
        self.run_check('--archive-dir', self.base / 'absent', status=2)

    def test_archive_traversal_is_rejected_without_extracting(self):
        archive = self.archive({'fixture/../../escaped.txt': 'not extracted'})
        self.run_check('--archive', archive, status=2)
        self.assertFalse((self.base / 'escaped.txt').exists())

    def test_archive_symlink_is_rejected(self):
        archive = self.base / 'linked.tar.gz'
        with tarfile.open(archive, 'w:gz') as tar:
            member = tarfile.TarInfo('fixture/link')
            member.type = tarfile.SYMTYPE
            member.linkname = '/etc/passwd'
            tar.addfile(member)
        self.run_check('--archive', archive, status=2)

    def test_internal_relative_source_symlink_is_checked(self):
        self.write('target.md', PRIVATE)
        (self.root / 'linked.md').symlink_to('target.md')
        report = self.run_check(status=1)
        self.assertEqual({x['path'] for x in report['findings']}, {'linked.md', 'target.md'})

    def test_external_source_symlink_is_rejected(self):
        external = self.base / 'outside.md'
        external.write_text('External material')
        (self.root / 'linked.md').symlink_to('../outside.md')
        self.run_check(status=2)

    def test_exception_requires_exact_content_and_expires_on_change(self):
        self.write('tools/example.py', 'text = "' + PRIVATE + '"\n')
        first = self.run_check(status=1)['findings'][0]
        entry = {key: first[key] for key in ['scope', 'path', 'rule', 'fingerprint']}
        entry['reason'] = 'Literal input used to test the reference detector.'
        self.write('tools/public-content-allowlist.json', json.dumps({'version': 1, 'exceptions': [entry]}))
        self.assertEqual(len(self.run_check()['exceptions']), 1)
        self.write('tools/example.py', 'text = "' + PRIVATE + '/changed"\n')
        report = self.run_check(status=1)
        self.assertEqual(len(report['stale_exceptions']), 1)
        self.assertEqual(len(report['findings']), 1)

    def test_exception_does_not_apply_to_another_path(self):
        self.write('methods.md', WORKFLOW)
        first = self.run_check(status=1)['findings'][0]
        entry = {key: first[key] for key in ['scope', 'path', 'rule', 'fingerprint']}
        entry['reason'] = 'Demonstration fixture for the wording detector.'
        self.write('tools/public-content-allowlist.json', json.dumps({'version': 1, 'exceptions': [entry]}))
        self.write('README.md', WORKFLOW)
        report = self.run_check(status=1)
        self.assertEqual([x['path'] for x in report['findings']], ['README.md'])

    def test_broad_exception_is_rejected(self):
        self.write('tools/public-content-allowlist.json', json.dumps({'version': 1, 'exceptions': [{
            'scope': 'repo', 'path': '*.md', 'rule': 'internal-workflow',
            'fingerprint': '0' * 64, 'reason': 'Not sufficiently specific.'}]}))
        self.run_check(status=2)

    def test_reports_cannot_be_written_inside_checkout(self):
        result = subprocess.run([sys.executable, str(CHECKER), '--repo', str(self.root),
                                 '--report-dir', str(self.root / 'ignored-reports')],
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, 2)
        self.assertFalse((self.root / 'ignored-reports').exists())

    def test_missing_git_repository_fails_closed(self):
        import shutil
        shutil.rmtree(self.root / '.git')
        self.run_check(status=2)

    def test_deleted_tracked_files_do_not_break_scan(self):
        path = self.write('old.md', PRIVATE)
        subprocess.run(['git', 'add', 'old.md'], cwd=self.root, check=True)
        path.unlink()
        self.run_check()


if __name__ == '__main__':
    unittest.main()
