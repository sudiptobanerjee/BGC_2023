#!/usr/bin/env python3
"""
Find files that are byte-identical duplicates elsewhere in the repo and replace
the duplicate copies with a relative symlink to one canonical copy.

Why symlinks instead of editing read.csv()/read.table() paths in the .Rmd/.R files:
Every script currently assumes its data lives in a local ./data (or similar) folder
relative to itself -- that's *why* the data got copied everywhere in the first place.
A symlink means the relative path each script already uses keeps resolving to a real
file, so zero R code has to change and Prof. Banerjee's existing local workflow is
completely unaffected. Git stores a symlink as a few bytes (the target path), not a
second copy of the file.

SAFE BY DESIGN: only files with an *exact* (hash-identical) match elsewhere are
touched. Nothing is deleted or rewritten based on filename alone.

Usage:
    python3 _site-config/dedupe_data.py            # dry run, prints the plan only
    python3 _site-config/dedupe_data.py --apply     # actually replace duplicates

Run this from your own terminal (not through the Cowork device bridge --
deleting/replacing files needs a real filesystem, which the bridge used in this
session could not do).
"""
import argparse
import hashlib
import os
import subprocess
from collections import defaultdict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent

# Don't touch anything under these -- .git internals, LFS pointers, or things that
# are duplicated on purpose (e.g. the Chapter*/ published copies, handled separately
# by the website generator, not this script).
SKIP_PREFIXES = (".git/", "website/")


def tracked_files():
    out = subprocess.run(
        ["git", "-C", str(REPO_ROOT), "ls-files"], capture_output=True, text=True
    ).stdout
    for line in out.splitlines():
        if line and not line.startswith(SKIP_PREFIXES):
            yield line


def file_hash(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def pick_canonical(paths):
    # Prefer a copy directly under the top-level data/ folder; then a copy under
    # src/ (real source, not a published Chapter*/ copy); otherwise the shortest
    # path (usually the "most central" copy).
    for p in paths:
        if p.startswith("data/"):
            return p
    src_paths = [p for p in paths if p.startswith("src/")]
    if src_paths:
        return min(src_paths, key=len)
    return min(paths, key=len)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true", help="actually replace duplicates with symlinks")
    args = ap.parse_args()

    by_hash = defaultdict(list)
    for rel in tracked_files():
        full = REPO_ROOT / rel
        if not full.is_file() or full.is_symlink():
            continue
        # Only bother hashing plausible data files -- skip source code/scripts.
        if full.suffix.lower() not in {".csv", ".txt", ".dat", ".rdata", ".shp", ".shx", ".dbf", ".sbn", ".sbx", ".img", ".rrd"}:
            continue
        by_hash[file_hash(full)].append(rel)

    dupes = {h: p for h, p in by_hash.items() if len(p) > 1}
    if not dupes:
        print("No exact-duplicate data files found.")
        return

    total_saved = 0
    for h, paths in sorted(dupes.items(), key=lambda kv: -len(kv[1])):
        canonical = pick_canonical(paths)
        others = [p for p in paths if p != canonical]
        size = (REPO_ROOT / canonical).stat().st_size
        total_saved += size * len(others)
        print(f"\ncanonical: {canonical}  ({size/1024:.0f} KB)")
        for other in others:
            other_path = REPO_ROOT / other
            rel_target = os.path.relpath(REPO_ROOT / canonical, other_path.parent)
            print(f"  -> symlink {other}  =>  {rel_target}")
            if args.apply:
                other_path.unlink()
                other_path.symlink_to(rel_target)

    print(f"\n{'Replaced' if args.apply else 'Would replace'} duplicates, "
          f"saving ~{total_saved/1024/1024:.1f} MB of tracked file content.")
    if not args.apply:
        print("Dry run only -- re-run with --apply to make the changes, "
              "then review with `git status` / `git diff` before committing.")


if __name__ == "__main__":
    main()
