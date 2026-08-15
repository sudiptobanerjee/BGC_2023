#!/usr/bin/env python3
"""
Build the Quarto website skeleton from _site-config/chapters.yaml.

This script does NOT render any R/Rmd. It only:
  1. Reads the chapter -> example manifest.
  2. Copies each already-rendered .html / static asset into website/chapterN/.
  3. Writes a chapterN/index.qmd listing links to those files.
  4. Writes the top-level index.qmd with a chapter list.

Run from the repo root:  python3 _site-config/generate_site.py
Requires: PyYAML (pip install pyyaml)
"""
import shutil
import sys
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent
MANIFEST = REPO_ROOT / "_site-config" / "chapters.yaml"
SITE_DIR = REPO_ROOT / "website"


def main():
    if not MANIFEST.exists():
        sys.exit(f"Manifest not found: {MANIFEST}")

    data = yaml.safe_load(MANIFEST.read_text())
    chapters = data.get("chapters", [])

    if SITE_DIR.exists():
        shutil.rmtree(SITE_DIR)
    SITE_DIR.mkdir(parents=True)

    missing = []
    index_lines = [
        "---",
        "title: \"Hierarchical Modeling and Analysis for Spatial Data\"",
        "subtitle: \"Companion code and worked examples\"",
        "---",
        "",
        "## Chapters",
        "",
    ]

    for chapter in chapters:
        num = chapter["number"]
        chap_dir = SITE_DIR / f"chapter{num}"
        chap_dir.mkdir(parents=True, exist_ok=True)

        index_lines.append(f"- [Chapter {num}](chapter{num}/index.qmd)")

        chap_lines = [
            "---",
            f"title: \"Chapter {num}\"",
            "---",
            "",
        ]

        for ex in chapter.get("examples", []):
            title = ex["title"]
            src_rel = ex.get("html") or ex.get("asset")
            if not src_rel:
                continue
            src_path = REPO_ROOT / src_rel

            if not src_path.exists() or src_path.suffix.lower() == ".rmd":
                # Either genuinely missing, or the manifest still points at an .Rmd
                # because no rendered output has been pushed yet.
                missing.append((num, title, src_rel))
                chap_lines.append(
                    f"- {title} -- *not yet available (source: `{src_rel}`)*"
                )
                continue

            dest_name = src_path.name
            dest_path = chap_dir / dest_name
            shutil.copy2(src_path, dest_path)
            chap_lines.append(f"- [{title}]({dest_name})")

        (chap_dir / "index.qmd").write_text("\n".join(chap_lines) + "\n")

    (SITE_DIR / "index.qmd").write_text("\n".join(index_lines) + "\n")

    # Copy the Quarto project config into the site working directory Quarto expects
    shutil.copy2(REPO_ROOT / "_quarto.yml", SITE_DIR / "_quarto.yml")

    print(f"Generated {len(chapters)} chapter pages under {SITE_DIR}")
    if missing:
        print("\nWARNING -- these examples have no rendered output to publish yet:")
        for num, title, path in missing:
            print(f"  Chapter {num}: {title}  (expected at {path})")
        print(
            "\nThese will show as 'not yet available' on the site until the "
            ".html is knitted locally and pushed."
        )


if __name__ == "__main__":
    main()
