#!/usr/bin/env python3

import io
from os.path import dirname, join
import extern

def get_version(relpath):
  """Read version info from a file without importing it"""
  for line in io.open(join(dirname(__file__), relpath), encoding="utf-8"):
    if "version" in line:
      if '"' in line:
        return line.split('"')[1]
      elif "'" in line:
        return line.split("'")[1]

if __name__ == "__main__":
    print("running release.sh")
    extern.run("./release.sh")

    version = get_version('Cargo.toml')
    print("version is {}".format(version))

    print("Checking if repo is clean ..")
    extern.run('if [[ $(git diff --shortstat 2> /dev/null | tail -n1) != "" ]]; then exit 1; fi')

    extern.run('git tag v{}'.format(version))
    print("Now run 'cargo publish && git push && git push --tags'".format(version))
