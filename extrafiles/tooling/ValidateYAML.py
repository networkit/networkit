#!/usr/bin/env python3
"""
This tool validates NetworKit YAML configuration files.
"""
import nktooling as nkt
import os
import yaml

nkt.setup()
os.chdir(nkt.getNetworKitRoot())

for configFile in [".clang-format", ".clang-tidy"]:
    try:
        with open(configFile, "r") as f:
            yaml.load(f, yaml.SafeLoader)
    except yaml.YAMLError as err:
        nkt.failIfMalformed(configFile, err)
