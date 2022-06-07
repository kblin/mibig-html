# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.
"""The MIBiG HTML generation pipeline"""

import logging
import os
import sys
from typing import List, Optional

import antismash
from antismash.common.path import (
    changed_directory,
    get_full_path,
    locate_file,
)
from antismash.common.subprocessing import execute

import mibig_html

GIT_VERSION_FALLBACK_FILENAME = get_full_path(__file__, "git_hash")


def get_git_version(fallback_filename: Optional[str] = GIT_VERSION_FALLBACK_FILENAME) -> str:
    """Get the sha1 of the current git version"""
    git_version = ""
    try:
        with changed_directory(os.path.dirname(__file__)):
            version_cmd = execute(['git', 'rev-parse', '--short', 'HEAD'])
            status_cmd = execute(['git', 'status', '--porcelain'])
        if version_cmd.successful() and status_cmd.successful():
            git_version = version_cmd.stdout.strip()
            changes = status_cmd.stdout.splitlines()
            if changes:
                git_version += "(changed)"
    except OSError:
        pass
    if git_version == "" and fallback_filename:
        if locate_file(fallback_filename, silent=True):
            with open(fallback_filename, 'rt') as handle:
                git_version = handle.read().strip()
    return git_version


def get_version() -> str:
    """Get the current version string"""
    version = antismash.__version__
    git_version = get_git_version()
    if git_version:
        version += "-%s" % git_version

    return version


# set arg version to avoid cyclic imports
antismash.config.args.ANTISMASH_VERSION = get_version()


def main(args: List[str]) -> int:
    """ The entrypoint of MIBiG-html as if it was on the command line

        Arguments:
            args: a list of args as would be given on the command line
                    e.g. ["inputfile", "--minimal", "--enable-nrps_pks"]

        Returns:
            zero if successful, non-zero otherwise

    """
    all_modules = mibig_html.get_all_modules()
    parser = antismash.config.args.build_parser(from_config_file=True, modules=all_modules)

    # if --help, show help texts and exit
    if set(args).intersection({"-h", "--help", "--help-showall"}):
        parser.print_help(None, "--help-showall" in args)
        return 0

    options = antismash.config.build_config(args, parser=parser)
    # pylint has some issues properly handling the config members, where mypy does not
    assert isinstance(options, antismash.config.ConfigType)
    # pylint: disable=no-member

    if options.write_config_file:
        parser.write_to_config_file(options.write_config_file, options.__dict__)
        return 0

    # if -V, show version text and exit
    if options.version:
        print("MIBiG-HTML %s" % get_version())
        return 0

    if len(options.sequences) > 1:
        parser.error("Only one sequence file should be provided")
    if len(options.sequences) < 1 and not options.reuse_results \
            and not options.check_prereqs_only and not options.list_plugins:
        parser.error("One of an input file or --reuse-results must be specified")
    if options.sequences and options.reuse_results:
        parser.error("Provide a sequence file or results to reuse, not both.")
    if options.sequences:
        sequence = options.sequences[0]
        if not os.path.exists(sequence):
            parser.error("Input file does not exist: %s" % sequence)
        if not os.path.isfile(sequence):
            parser.error("input %s is not a file" % sequence)
    else:
        sequence = ""

    if options.reuse_results and not os.path.exists(options.reuse_results):
        parser.error("Input file does not exist: %s" % options.reuse_results)

    if os.sep in options.output_basename:
        parser.error("Output basename cannot contain a path separator character {}".format(os.sep))

    options.version = get_version()

    try:
        mibig_html.run_mibig(sequence, options)
    except antismash.common.errors.AntismashInputError as err:
        if not str(err):
            raise
        logging.error(str(err))
        return 1

    return 0
    # pylint: enable=no-member


def entrypoint() -> None:
    """This is needed for the script generated by setuptools."""
    with antismash.common.logs.changed_logging():
        sys.exit(main(sys.argv[1:]))


if __name__ == '__main__':
    entrypoint()
