"""This generates exe_metadata.txt. It reads in the existing template and
replaces the placeholders with the hs version number and git commit.
The output is saved to the build folder and used by PyInstaller to make the windows
executables. This script is intended to be run on the remote when the Github
workflow build_windows_exe.yml is triggered.
"""

import sys
import subprocess
from importlib.metadata import version
from pathlib import Path
from packaging.version import Version

version_string = version("heatsource9")
major, minor, micro = Version(version_string).release
version_tuple = (major, minor, micro, 0)

# Get the seven character commit that triggered the build
git_commit = subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()[:7]

template = Path(__file__).with_name("exe_metadata.txt")
exe_metadata = template.read_text(encoding="utf-8")
exe_metadata = exe_metadata.replace("__VERSION_TUPLE__", str(version_tuple))
exe_metadata = exe_metadata.replace("__VERSION_STRING__", version_string)
exe_metadata = exe_metadata.replace("__GIT_COMMIT__", git_commit)

output = Path(sys.argv[1])
output.parent.mkdir(parents=True, exist_ok=True)
output.write_text(exe_metadata, encoding="utf-8")