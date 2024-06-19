#!/usr/bin/env python
# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

"""
This script checks all files in the repo for the correct copyright statement.
"""
from __future__ import print_function
import argparse
import functools
import os.path
import re
import subprocess
import sys


class ErrCollector:
    def __init__(self, verbose):
        self.verbose = verbose
        # self.error_occurred = False

    def info(self, filename, msg):
        if self.verbose:
            print(filename, msg)

    def error(self, filename, msg):
        # self.error_occurred = True
        print(filename, msg)


py_coding_re = re.compile(r"#.*coding[=:]\s*([-\w.]+)")

comment_chars = {
    ".applescript": "--",
    ".cc": "//",
    ".cmake": "#",
    ".cpp": "//",
    ".sh": "#",
    ".i": "//",
    ".hpp": "//",
    ".h": "//",
    ".py": "#",
    ".pyx": "#",
    ".pxd": "#",
    "README": "",
}

excluded_files = set(
    (
        ".github/workflows/run-clang-format.py",
        "geometry-tool/HlbGmyTool/Model/GmyGeneration/Point_inside_polyhedron_3.h",
        "geometry-tool/HlbGmyTool/Model/CommonGeneration/PybindVTKTypeCaster.h",
        "geometry-tool/HlbGmyTool/Model/GmyGeneration/Triangle_3_Ray_3_do_intersect.h",
    )
)


def comment_char(filename):
    base = os.path.basename(filename)
    root, ext = os.path.splitext(base)
    if ext == ".in":
        ext = os.path.splitext(root)[1]
    if base == "CMakeLists.txt":
        ext = ".cmake"

    return comment_chars[ext]


def is_source(filename):
    try:
        cc = comment_char(filename)
        return True
    except KeyError:
        return False


def is_in_git(filename):
    git_ls_cmd = ["git", "ls-tree", "--name-only", "--full-name", "HEAD", filename]
    output = subprocess.check_output(git_ls_cmd).strip()
    if not isinstance(output, str):
        # Py3
        output = output.decode()
    return os.path.abspath(output) == os.path.abspath(filename)


_re_cache = {}


def make_emacs_varline_re(cc):
    try:
        return _re_cache[cc]
    except KeyError:
        # Pattern is: comment char; optional whitespace; literal '-*-';
        # anything; literal '-*-'; optional whitespace; end of string.
        ans = _re_cache[cc] = re.compile(cc + r"\s*-\*-.*-\*-\s*$")
        return ans


def get_initial_comment(filename, file_or_lines):
    root, ext = os.path.splitext(filename)

    cc = comment_char(filename)
    ncc = len(cc)

    # State vars
    needs_shebang_checked = ext in (".py", ".sh")
    needs_py_encoding_checked = ext == ".py"
    needs_emacs_varline_checked = True
    emacs_varline_re = make_emacs_varline_re(cc)

    # Candidate CR message
    file_msg = []

    for line in file_or_lines:
        if needs_shebang_checked:
            needs_shebang_checked = False
            if line.startswith("#!"):
                # discard
                continue

        if needs_py_encoding_checked:
            # .py may have an encoding spec in lines 1 or 2
            needs_py_encoding_checked = False
            if py_coding_re.match(line):
                # discard
                continue

        if needs_emacs_varline_checked:
            # Files may have emacs file variables set
            needs_emacs_varline_checked = False
            if emacs_varline_re.match(line):
                # discard
                continue

        # Commented message must begin immediately
        if line.startswith(cc):
            file_msg.append(line[ncc:].strip())
        else:
            break
    return file_msg


cr_msg = """This file is part of HemeLB and is Copyright (C)
the HemeLB team and/or their institutions, as detailed in the
file AUTHORS. This software is provided under the terms of the
license in the file LICENSE."""


def try_fix(log, filename, dryrun=True):
    root, ext = os.path.splitext(filename)
    cc = comment_char(filename)
    ncc = len(cc)
    commented_msg = [cc + " " + x for x in cr_msg.split("\n")]

    with open(filename) as f:
        wrong = f.readlines()

    needs_shebang_checked = ext in (".py", ".sh")
    needs_py_encoding_checked = ext == ".py"
    needs_emacs_varline_checked = True
    emacs_varline_re = make_emacs_varline_re(cc)

    discard_blank = True
    discard_comment = True

    out = []
    for line in wrong:
        if needs_shebang_checked:
            # .sh and .py files can have a #!
            needs_shebang_checked = False
            if line.startswith("#!"):
                out.append(line)
                continue

        if needs_py_encoding_checked:
            # .py may have an encoding spec in lines 1 or 2
            needs_py_encoding_checked = False
            if py_coding_re.match(line):
                out.append(line)
                continue

        if needs_emacs_varline_checked:
            needs_emacs_varline_checked = False
            if line.startswith(cc):
                if emacs_varline_re.match(line):
                    # discard
                    continue

        if discard_blank:
            if line.strip() == "":
                continue
            discard_blank = False

        if discard_comment:
            if line.startswith(cc):
                continue

            discard_comment = False
            for x in cr_msg.split("\n"):
                out.append(cc + " " + x + "\n")

        out.append(line)

    with open(filename, "w") as f:
        f.write("".join(out))
    return True


def check_cr(log, filename, fix=False):
    if not is_in_git(filename):
        log.info(filename, "Skipping non-VCed file")
        return True
    if not is_source(filename):
        log.info(filename, "Skipping non-source file")
        return True

    if filename in excluded_files:
        log.info(filename, "Skipping file as not checked")
        return True

    log.info(filename, "Checking file")

    with open(filename) as f:
        msg_lines = get_initial_comment(filename, f)

    msg = "\n".join(msg_lines)
    if msg != cr_msg:
        log.error(filename, "copyright statement doesn't read as expected")
        if fix:
            return try_fix(log, filename)
        return False
    else:
        log.info(filename, "OK")
        return True


parser = argparse.ArgumentParser()
parser.add_argument("--verbose", "-v", action="count", default=0)
parser.add_argument("--fix", action="store_true")
parser.add_argument("paths", nargs="*")
import pdb


def iter_repo_paths(paths):
    for p in paths:
        p = os.path.abspath(p)
        if os.path.isdir(p):
            for dirpath, dirnames, filenames in os.walk(p):
                if ".git" in dirnames:
                    dirnames.remove(".git")
                builds = [d for d in dirnames if d.startswith("build")]
                for b in builds:
                    dirnames.remove(b)

                for fn in filenames:
                    yield os.path.relpath(os.path.join(dirpath, fn), repo_root)
        else:
            yield p


if __name__ == "__main__":
    args = parser.parse_args()

    log = ErrCollector(args.verbose > 0)
    error_occurred = False

    self_path = os.path.abspath(__file__)
    scripts_dir = os.path.dirname(self_path)
    repo_root = os.path.dirname(scripts_dir)
    if len(args.paths) == 0:
        # No paths to check, do the whole repo
        args.paths = [repo_root]

    n_fail = 0
    for p in iter_repo_paths(args.paths):
        ok = check_cr(log, p, args.fix)
        if ok:
            continue
        n_fail += 1

    sys.exit(n_fail)
