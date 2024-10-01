import sys
import os


def parse_args():
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('input', nargs='+', type=str,
                        help='matched files to call DEL IDs for')

    get_args(parser, "calling")
    return parser.parse_args(sys.argv[1:])


if __name__ == '__main__':
    from deli.scripts import run_calling, get_args, check_args
    from deli.logger import setup_logger
    from deli.html_report import HTMLReport

    args = parse_args()
    check_args(args)

    logger = setup_logger("sequence-caller", debug=args.debug)
    report = HTMLReport()  # for compatability, no report will be made

    os.makedirs(args.outdir, exist_ok=True)

    run_calling(args, logger, report)
