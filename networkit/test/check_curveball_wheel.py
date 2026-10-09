"""Check an installed repaired wheel in an isolated environment, without sampling.

New verification code under NetworKit's existing MIT terms. It does not copy
the planner implementation or change the Apache component's license.
"""
import argparse
from dataclasses import FrozenInstanceError
from fractions import Fraction
import hashlib
from importlib import metadata, resources
from importlib.machinery import EXTENSION_SUFFIXES
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import venv
import zipfile


NOTICE_HASHES = {
    "_licenses/curveball-planner/LICENSE": "c71d239df91726fc519c6eb72d318ec65820627232b2f796219e87dcf35d0ab4",
    "_licenses/curveball-planner/NOTICE": "2a2708df253f2c2cae0341a469ed60a4acfe4499a76e2b71e77139f79655a6ac",
}


def installedCheck(sourceRoot):
    import networkit as nk
    import networkit.randomization as randomization

    packagePath = Path(nk.__file__).resolve()
    extensionPath = Path(randomization.__file__).resolve()
    prefix = Path(sys.prefix).resolve()
    if packagePath.is_relative_to(sourceRoot) or extensionPath.is_relative_to(sourceRoot):
        raise RuntimeError("source-checkout import cannot validate an installed wheel")
    if not packagePath.is_relative_to(prefix) or not extensionPath.is_relative_to(prefix):
        raise RuntimeError("imports must originate in the fresh wheel environment")
    if not any(str(extensionPath).endswith(suffix) for suffix in EXTENSION_SUFFIXES):
        raise RuntimeError("randomization must be the compiled extension")

    plan = randomization.curveballTradePlan([2] * 6, Fraction(1, 19900), maxTrades=314)
    if type(plan) is not randomization.CurveballTradePlan:
        raise AssertionError("unexpected public record type")
    if plan.attemptedTrades != 315 or plan.withinTradeLimit is not False:
        raise AssertionError("installed planner count/cap mismatch")
    if plan.degrees != (2,) * 6 or plan.epsilon != Fraction(1, 19900):
        raise AssertionError("installed planner input record mismatch")
    try:
        plan.attemptedTrades = 0
    except FrozenInstanceError:
        pass
    else:
        raise AssertionError("installed planner record is mutable")

    observed = {}
    for name, expected in NOTICE_HASHES.items():
        observed[name] = hashlib.sha256(resources.files("networkit").joinpath(name).read_bytes()).hexdigest()
        if observed[name] != expected:
            raise AssertionError("installed component notice mismatch: " + name)
    distributions = metadata.packages_distributions().get("networkit", [])
    matches = [metadata.distribution(name) for name in distributions
               if name.replace("_", "-").lower() in ("networkit", "networkit-nightly")]
    if len(matches) != 1:
        raise AssertionError("expected one installed NetworKit distribution")
    distribution = matches[0]
    classifiers = distribution.metadata.get_all("Classifier", [])
    for classifier in ("License :: OSI Approved :: MIT License",
                       "License :: OSI Approved :: Apache Software License"):
        if classifier not in classifiers:
            raise AssertionError("missing component license classifier")
    licenseText = " ".join(distribution.metadata.get_all("License", [])
                           + distribution.metadata.get_all("License-Expression", []))
    if "MIT" not in licenseText or "Apache" not in licenseText:
        raise AssertionError("distribution metadata must identify both component licenses")
    print(json.dumps({"check": "installed_curveball_wheel", "passed": True,
                      "networkit_origin": str(packagePath), "randomization_origin": str(extensionPath),
                      "installed_notice_sha256": observed, "license_metadata_checked": True,
                      "public_api_count_cap_type_immutability_checked": True,
                      "sampler_invocations": 0}, sort_keys=True))


def wheelCheck(wheelDirectory, sourceRoot):
    wheels = sorted(wheelDirectory.glob("*.whl"))
    if len(wheels) != 1:
        raise RuntimeError("expected exactly one repaired wheel")
    wheel = wheels[0].resolve()
    with zipfile.ZipFile(wheel) as archive:
        if archive.testzip() is not None:
            raise AssertionError("corrupt wheel archive")
        for name, expected in NOTICE_HASHES.items():
            if hashlib.sha256(archive.read("networkit/" + name)).hexdigest() != expected:
                raise AssertionError("wheel component notice mismatch: " + name)
        archive.getinfo("networkit/_curveball_planning.py")

    with tempfile.TemporaryDirectory(prefix="curveball-wheel-check-") as temporary:
        temporaryPath = Path(temporary).resolve()
        if temporaryPath.is_relative_to(sourceRoot):
            raise RuntimeError("wheel check requires a directory outside the checkout")
        if not temporaryPath.is_relative_to(Path(tempfile.gettempdir()).resolve()):
            raise RuntimeError("unexpected temporary directory")
        environmentPath = temporaryPath / "venv"
        venv.EnvBuilder(with_pip=True, system_site_packages=False).create(environmentPath)
        python = environmentPath / ("Scripts/python.exe" if os.name == "nt" else "bin/python")
        environment = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
                           MKL_NUM_THREADS="1", PYTHONNOUSERSITE="1")
        environment.pop("PYTHONPATH", None)
        subprocess.run([str(python), "-I", "-B", "-m", "pip", "--isolated", "install",
                        "--no-cache-dir", "--no-input",
                        "--disable-pip-version-check", "--only-binary=:all:", str(wheel)],
                       cwd=temporaryPath, env=environment, check=True)
        subprocess.run([str(python), "-I", "-B", str(Path(__file__).resolve()), "--installed",
                        "--source-root", str(sourceRoot)], cwd=temporaryPath, env=environment, check=True)
    print(json.dumps({"check": "repaired_curveball_wheel", "passed": True,
                      "wheel_name": wheel.name, "wheel_sha256": hashlib.sha256(wheel.read_bytes()).hexdigest(),
                      "archive_notice_bytes_checked": True, "sampler_invocations": 0}, sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--wheel-dir", type=Path)
    parser.add_argument("--installed", action="store_true")
    parser.add_argument("--source-root", type=Path, default=Path(__file__).resolve().parents[2])
    arguments = parser.parse_args()
    if arguments.installed == (arguments.wheel_dir is not None):
        parser.error("choose exactly one of --installed or --wheel-dir")
    if arguments.installed:
        installedCheck(arguments.source_root.resolve())
    else:
        wheelCheck(arguments.wheel_dir.resolve(), arguments.source_root.resolve())
