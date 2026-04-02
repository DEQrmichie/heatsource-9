
"""
Command line entrypoints: `python -m heatsource9` or `hs`

Equivalent command forms:

  hs [-md MODEL_DIR] run -t | -s | -hy
  hs [-md MODEL_DIR] setup (-cf | -mi) [setup options]
  python -m heatsource9  [-md MODEL_DIR] run -t | -s | -hy
  python -m heatsource9 [-md MODEL_DIR] setup (-cf | -mi) [setup options]

Global options:

  -h / --help              : show help
  -md / --model-dir        : path to the model directory
  -v / --version           : print Heat Source version

Run types:

  -t / --temperature        : run the full temperature model
  -s / --solar              : run solar routines only
  -hy / --hydraulics        : run hydraulics only

Run options:
  (none)

Setup types:

  -cf / --control-file      : write a blank control file template
  -mi / --model-inputs      : write blank input files from a parameterized control file

Setup options:

  -csv / --csv-mode         : with -cf, write a CSV control file instead of XLSX
  -t / --timestamp          : add a timestamp to the template file name
  -o / --overwrite          : overwrite existing files (default is to keep existing files)

Usage examples:
`hs setup -cf -md /path/to/model_dir` : writes a blank control file in the model directory.
`hs setup -mi -o` : writes blank input files from a parameterized control file and overwrites any existing files.
`hs run -t` : runs a temperature model using `HeatSource_Control.*` found in the selected model directory.
"""


import argparse
import sys
from os.path import abspath

import heatsource9.run as hs_run
from heatsource9.setup import (
    write_cf,
    write_mi,
)


def _get_model_dir(args):
    return abspath(getattr(args, "model_dir", None) or ".")


def _parse_args(argv):
    # Build command line options for global settings and subcommands.
    parser = argparse.ArgumentParser(prog="hs", usage="hs <command> [options]")
    parser.add_argument(
        "-md",
        "--model-dir",
        nargs="?",
        default=abspath("."),
        help="Path to the model directory. Default is current working directory.",
    )
    parser.add_argument(
        "-v",
        "--version",
        action="store_true",
        help="Print installed Heat Source version.",
    )

    subparsers = parser.add_subparsers(dest="command")

    # Add command options for model runs.
    run_parser = subparsers.add_parser("run", help="Run a model with -t | -s | -hy")
    rt = run_parser.add_mutually_exclusive_group(required=True)
    rt.add_argument("-t", "--temperature", action="store_true", help="Runs a temperature model.")
    rt.add_argument("-s", "--solar", action="store_true", help="Runs solar routines only.")
    rt.add_argument("-hy", "--hydraulics", action="store_true", help="Runs hydraulics only.")
    run_parser.add_argument(
        "-md",
        "--model-dir",
        nargs="?",
        default=None,
        help="Path to the model directory. Overrides global --model-dir when provided after 'run'.",
    )
    # Add command options for setup.
    setup_parser = subparsers.add_parser("setup", help="Setup model templates with -cf | -mi")
    setup_group = setup_parser.add_mutually_exclusive_group(required=True)
    setup_group.add_argument("-cf", "--control-file", action="store_true", help="Writes a blank control file.")
    setup_parser.add_argument(
        "-md",
        "--model-dir",
        nargs="?",
        default=None,
        help="Path to the model directory. Overrides global --model-dir when provided after 'setup'.",
    )
    setup_group.add_argument(
        "-mi",
        "--model-inputs",
        action="store_true",
        help="Writes blank input files from a parameterized control file.",
    )
    setup_parser.add_argument(
        "-csv",
        "--csv-mode",
        action="store_true",
        help="Write CSV templates instead of XLSX for control file setup.",
    )
    setup_parser.add_argument(
        "-t",
        "--timestamp",
        action="store_true",
        help="Prefix output template filenames with timestamp.",
    )
    setup_parser.add_argument(
        "-o",
        "--overwrite",
        action="store_true",
        help="Overwrite existing files.",
    )
    args = parser.parse_args(argv)
    if not args.version and args.command is None:
        parser.error("the following arguments are required: command")
    return args


def main(argv = None):
    argv = list(sys.argv[1:] if argv is None else argv)
    args = _parse_args(argv)

    # model version
    if getattr(args, "version", False):
        from heatsource9.__version__ import __version__
        print(f"hs {__version__}")
        return 0

    # Setup template writing
    if args.command == "setup":
        model_dir = _get_model_dir(args)
        if args.control_file:
            write_cf(
                model_dir=model_dir,
                use_timestamp=bool(args.timestamp),
                overwrite=bool(args.overwrite),
                csv_mode=bool(args.csv_mode),
            )
        if args.model_inputs:
            write_mi(
                model_dir=model_dir,
                use_timestamp=bool(args.timestamp),
                overwrite=bool(args.overwrite),
            )
        return 0

    # Run command starts here.
    model_dir = _get_model_dir(args)

    if args.temperature:
        hs_run.temperature(model_dir)
    elif args.solar:
        hs_run.solar(model_dir)
    else:
        hs_run.hydraulics(model_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
