import argparse
import sys
import os


def parse_args():
    from deli.scripts import get_args

    parser = argparse.ArgumentParser()

    parser.add_argument('input', type=str,
                        help='input fastq file with reads to be searched for matches')

    get_args(parser, "matching")

    return parser.parse_args(sys.argv[1:])


if __name__ == '__main__':
    from deli.scripts import run_matching, check_args
    from deli.logger import setup_logger
    from deli.html_report import HTMLReport

    args = parse_args()
    check_args(args, check_lib=False, check_index=False)
    logger = setup_logger("barcode-matching", debug=args.debug)
    report = HTMLReport()  # for compatability, no report will be made

    os.makedirs(args.outdir, exist_ok=True)

    run_matching(args, logger, report)
