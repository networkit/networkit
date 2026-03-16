#!/usr/bin/env python3

import argparse
import datetime
import json
from pathlib import Path


def normalize_base_path(base_path: str) -> str:
    cleaned = base_path.strip()
    if not cleaned:
        return "/"
    if not cleaned.startswith("/"):
        cleaned = "/" + cleaned
    return cleaned.rstrip("/") or "/"


def build_channel_url(base_path: str, channel: str) -> str:
    if channel == "stable":
        return f"{base_path}/"
    return f"{base_path}/{channel}/"


def main() -> None:
    parser = argparse.ArgumentParser(description="Update docs version metadata.")
    parser.add_argument("--file", required=True, help="Path to versions.json")
    parser.add_argument("--channel", required=True, choices=("stable", "nightly"))
    parser.add_argument("--label", required=True)
    parser.add_argument("--base-path", required=True)
    parser.add_argument("--nightly-date", default="")
    args = parser.parse_args()

    metadata_path = Path(args.file)
    base_path = normalize_base_path(args.base_path)

    data = {}
    if metadata_path.exists():
        with metadata_path.open("r", encoding="utf-8") as infile:
            data = json.load(infile)

    channels = data.get("channels")
    if not isinstance(channels, dict):
        channels = {}

    stable_entry = channels.get("stable")
    if not isinstance(stable_entry, dict):
        stable_entry = {}

    nightly_entry = channels.get("nightly")
    if not isinstance(nightly_entry, dict):
        nightly_entry = {}

    stable_entry["url"] = build_channel_url(base_path, "stable")
    nightly_entry["url"] = build_channel_url(base_path, "nightly")

    channels["stable"] = stable_entry
    channels["nightly"] = nightly_entry

    channel_entry = channels[args.channel]
    channel_entry["label"] = args.label
    channel_entry["updated_at"] = datetime.datetime.now(
        datetime.timezone.utc
    ).isoformat(timespec="seconds")

    if args.channel == "nightly" and args.nightly_date:
        channel_entry["date"] = args.nightly_date

    data["base_path"] = base_path
    data["default_channel"] = "stable"
    data["channels"] = channels
    data["generated_at"] = datetime.datetime.now(datetime.timezone.utc).isoformat(
        timespec="seconds"
    )

    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    with metadata_path.open("w", encoding="utf-8") as outfile:
        json.dump(data, outfile, indent=2, sort_keys=True)
        outfile.write("\n")


if __name__ == "__main__":
    main()
