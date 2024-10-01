import os
import sys


def parse_args():
    import argparse
    from deli.scripts import get_args

    parser = argparse.ArgumentParser()

    parser.add_argument('input', type=str,
                        help='called sequence file to convert into cube format')
    get_args(parser, "cube_gen")
    return parser.parse_args(sys.argv[1:])


if __name__ == '__main__':
    from deli.scripts import run_cube_gen, check_args
    from deli.logger import setup_logger
    from deli.html_report import HTMLReport

    args = parse_args()
    check_args(args)

    logger = setup_logger("cube-generation", debug=args.debug)
    report = HTMLReport()  # for compatability, no report will be made

    os.makedirs(args.outdir, exist_ok=True)

    run_cube_gen(args, logger, report)
