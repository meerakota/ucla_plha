"""Updates the versioned documentation in a checkout of the gh-pages branch.

The documentation of each release is kept in its own folder, named by the major and minor
version (e.g. 2.1/, 3.0/). This script

1. moves a documentation site found at the root of gh-pages (the unversioned documentation
   published before versioning) into the folder of its version, so it is kept,
2. copies a new build, if given, into the folder of its version, replacing the earlier build,
3. writes versions.json, the list of versions read by the version switcher of the docs, and
4. writes redirect pages at the root: index.html and every page of the documentation redirect
   to the same page of the latest version, so links to the unversioned documentation keep working.

Usage:
    python update_docs_versions.py SITE [--build BUILD --version VERSION] [--base-url URL]
"""

import argparse
import json
import pathlib
import re
import shutil

VERSION_DIR = re.compile(r"^\d+\.\d+$")
REDIRECT_MARKER = "<!-- ucla_plha docs version redirect -->"
# files at the root of gh-pages that are not part of a documentation site
KEEP_AT_ROOT = {".git", ".nojekyll", "CNAME", "versions.json"}


def version_key(name):
    return tuple(int(part) for part in name.split("."))


def is_redirect(path):
    if path.suffix != ".html":
        return False
    with open(path, encoding="utf-8", errors="ignore") as f:
        return REDIRECT_MARKER in f.read(500)


def migrate_root_site(site):
    """Moves an unversioned documentation site at the root into the folder of its version."""
    options = site / "_static" / "documentation_options.js"
    if not (site / "objects.inv").exists() or not options.exists():
        return None
    match = re.search(r"VERSION:\s*'(\d+\.\d+)", options.read_text(encoding="utf-8"))
    if match is None:
        raise RuntimeError(f"no version found in {options}")
    version = match.group(1)
    target = site / version
    if target.exists():
        raise RuntimeError(f"{target} already exists; not moving the root site into it")
    target.mkdir()
    for path in list(site.iterdir()):
        if path.name in KEEP_AT_ROOT or path == target or VERSION_DIR.match(path.name):
            continue
        if path.is_file() and is_redirect(path):
            continue
        shutil.move(str(path), str(target / path.name))
    return version


def copy_build(site, build, version):
    if not VERSION_DIR.match(version):
        raise ValueError(f'version must be "major.minor", not "{version}"')
    target = site / version
    if target.exists():
        shutil.rmtree(target)
    shutil.copytree(build, target, ignore=shutil.ignore_patterns(".buildinfo", ".doctrees"))


def remove_redirects(site):
    for path in site.rglob("*.html"):
        rel = path.relative_to(site)
        if rel.parts[0] in KEEP_AT_ROOT or VERSION_DIR.match(rel.parts[0]):
            continue
        if is_redirect(path):
            path.unlink()
    # remove folders left empty by the redirects
    for path in sorted(site.rglob("*"), key=lambda p: len(p.parts), reverse=True):
        rel = path.relative_to(site)
        if rel.parts[0] in KEEP_AT_ROOT or VERSION_DIR.match(rel.parts[0]):
            continue
        if path.is_dir() and not any(path.iterdir()):
            path.rmdir()


def write_redirect(path, target):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        f"""<!DOCTYPE html>
{REDIRECT_MARKER}
<html>
  <head>
    <meta charset="utf-8">
    <title>Redirecting to the latest ucla_plha documentation</title>
    <meta http-equiv="refresh" content="0; url={target}">
    <link rel="canonical" href="{target}">
    <script>window.location.replace("{target}" + window.location.hash);</script>
  </head>
  <body>
    <p>The documentation has moved to <a href="{target}">{target}</a>.</p>
  </body>
</html>
""",
        encoding="utf-8",
    )


def documentation_pages(folder):
    """Relative paths of the html pages of a documentation folder, except Sphinx's own folders."""
    for path in folder.rglob("*.html"):
        rel = path.relative_to(folder)
        if not rel.parts[0].startswith("_") or rel.parts[0] == "_modules":
            yield rel


def update(site, base_url):
    versions = sorted(
        (p.name for p in site.iterdir() if p.is_dir() and VERSION_DIR.match(p.name)),
        key=version_key,
        reverse=True,
    )
    if not versions:
        raise RuntimeError(f"no version folders in {site}")
    latest = versions[0]
    base = base_url.rstrip("/")
    switcher = [
        {
            "name": f"{v} (latest)" if v == latest else v,
            "version": v,
            "url": f"{base}/{v}/",
            **({"preferred": True} if v == latest else {}),
        }
        for v in versions
    ]
    (site / "versions.json").write_text(json.dumps(switcher, indent=2) + "\n", encoding="utf-8")

    # redirect pages at the root: each page goes to the newest version that has it
    remove_redirects(site)
    pages = {}
    for v in reversed(versions):
        for rel in documentation_pages(site / v):
            pages[rel] = v
    for rel, v in pages.items():
        target = "../" * (len(rel.parts) - 1) + f"{v}/{rel.as_posix()}"
        write_redirect(site / rel, target)
    write_redirect(site / "index.html", f"{latest}/")
    (site / ".nojekyll").touch()
    return versions


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("site", type=pathlib.Path, help="checkout of the gh-pages branch")
    parser.add_argument("--build", type=pathlib.Path, help="Sphinx html build to publish")
    parser.add_argument("--version", help='version of the build, "major.minor"')
    parser.add_argument("--base-url", default="https://meerakota.github.io/ucla_plha")
    args = parser.parse_args()
    if (args.build is None) != (args.version is None):
        parser.error("--build and --version must be given together")

    migrated = migrate_root_site(args.site)
    if migrated:
        print(f"moved the unversioned documentation at the root into {migrated}/")
    if args.build is not None:
        copy_build(args.site, args.build, args.version)
        print(f"published {args.build} as {args.version}/")
    versions = update(args.site, args.base_url)
    print(f"versions: {', '.join(versions)} (latest {versions[0]})")


if __name__ == "__main__":
    main()
