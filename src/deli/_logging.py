"""for basic python-based logging"""

import logging
import os
import sys


def get_dummy_logger() -> logging.Logger:
    """Get a dummy logger; for compatability"""
    return logging.getLogger("dummy")


def get_logger(name: str = "root", debug: bool = False) -> logging.Logger:
    """
    Create a logger with the passed settings

    Notes
    -----
    Will always save logs to the CWD as 'deli.log'

    Parameters
    ----------
    name: str
        name to use for logging
    debug: bool, Default False
        if True, turn on debug logging

    Returns
    -------
    logging.Logger
    """
    logging.captureWarnings(True)  # hook python warning to the logger
    logger = logging.getLogger(name=name)

    # setup handlers
    file_handler = logging.FileHandler(os.path.join(os.getcwd(), "deli.log"))
    file_handler.setFormatter(
        logging.Formatter("%(asctime)s | %(name)s | %(levelname)s | %(message)s")
    )
    # set levels
    if debug:
        logger.setLevel("DEBUG")
        file_handler.setLevel("DEBUG")
    else:
        logger.setLevel("INFO")
        file_handler.setLevel("INFO")

    # add handlers to logger
    logger.addHandler(file_handler)

    # add an exception hook to the logger to allow for uncaught exceptions to be logged
    def handle_exception(exc_type, exc_value, exc_traceback):
        if issubclass(exc_type, KeyboardInterrupt):
            sys.__excepthook__(exc_type, exc_value, exc_traceback)
            return
        logger.critical("Uncaught exception", exc_info=(exc_type, exc_value, exc_traceback))
        sys.__excepthook__(exc_type, exc_value, exc_traceback)

    sys.excepthook = handle_exception

    return logger
