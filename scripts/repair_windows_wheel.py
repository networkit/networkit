import argparse
import pathlib
import subprocess
import sys


def find_windows_package_dir():
	build_root = pathlib.Path("build")
	candidates = sorted(build_root.glob("lib.win-*/networkit"))
	for candidate in candidates:
		if (candidate / "networkit_state.dll").is_file():
			return candidate
	raise FileNotFoundError("Unable to find build/lib.win-*/networkit/networkit_state.dll")


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("--wheel", required=True)
	parser.add_argument("--dest-dir", required=True)
	args = parser.parse_args()

	package_dir = find_windows_package_dir()
	subprocess.check_call([
		sys.executable,
		"-m",
		"delvewheel",
		"repair",
		"-w",
		args.dest_dir,
		"-v",
		"--add-path",
		str(package_dir),
		args.wheel,
	])


if __name__ == "__main__":
	main()
