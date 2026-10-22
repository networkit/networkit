#!/usr/bin/env python3

import argparse
import json
import subprocess
import urllib.request

from packaging.version import InvalidVersion, Version

API_URL = "https://api.anaconda.org/package/{owner}/{package}"


def published_versions(owner: str, package: str) -> set[str]:
    url = API_URL.format(owner=owner, package=package)
    with urllib.request.urlopen(url, timeout=60) as response:
        package_info = json.load(response)
    return {release["version"] for release in package_info["releases"]}


def split_versions(versions: set[str], keep: int) -> tuple[list[str], list[str]]:
    """Return the newest `keep` versions and the older ones, both newest first."""
    ordered = []
    for version in versions:
        try:
            ordered.append((Version(version), version))
        except InvalidVersion:
            # We never remove what we cannot order.
            print(f"Ignoring {version}: not a valid version")
    ordered.sort(reverse=True)
    names = [version for _, version in ordered]
    return names[:keep], names[keep:]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Remove all but the newest versions of a package on anaconda.org."
    )
    parser.add_argument("--owner", required=True, help="User or organization")
    parser.add_argument("--package", required=True)
    parser.add_argument("--keep", type=int, required=True)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    if args.keep < 1:
        parser.error("--keep must be at least 1")

    versions = published_versions(args.owner, args.package)
    kept, removed = split_versions(versions, args.keep)
    for version in kept:
        print(f"Keeping {args.owner}/{args.package}/{version}")
    for version in removed:
        spec = f"{args.owner}/{args.package}/{version}"
        print(f"Removing {spec}")
        if not args.dry_run:
            # anaconda-client takes the token from ANACONDA_API_TOKEN.
            subprocess.run(["anaconda", "remove", "--force", spec], check=True)


if __name__ == "__main__":
    main()
