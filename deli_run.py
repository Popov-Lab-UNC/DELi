import argparse
import sys
import os
from multiprocessing import cpu_count


def prep_args():
    from deli.scripts import get_args

    parser = argparse.ArgumentParser()

    parser.add_argument('input', type=str,
                        help='input fastq file with reads to be searched for matches')

    get_args(parser, "all")
    return parser.parse_args(sys.argv[1:])


if __name__ == '__main__':
    from deli.logger import setup_logger
    from deli.html_report import HTMLReport
    from deli.scripts import run_cube_gen, run_matching, run_calling, check_args

    report = HTMLReport()

    args = prep_args()
    check_args(args)

    logger = setup_logger('deli', debug=args.debug)
    logger.debug("using args {}".format(vars(args)))

    os.makedirs(args.outdir, exist_ok=True)
    logger.info(f"saving files to {args.outdir}")

    _pool_size = args.n_workers - 1 if args.n_workers != -1 else cpu_count() - 1
    args.n_workers = _pool_size
    logger.info("number of worker is {}".format(_pool_size))

    _match_files = run_matching(args, logger, report)
    _call_file = run_calling(args, logger, report, _match_files)
    run_cube_gen(args, logger, report, _call_file)

    report.get_html_report(out_path=os.path.join(args.outdir, args.prefix + "deli_report.html"))
    logger.info(f"DELi run for {args.input} completed")
