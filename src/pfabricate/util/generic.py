import os
import sys
import datetime


def print_header(title):
    """
    Print header for script

    params
        title : str
            Title of header.
    returns
        t0 : datetime.datetime
            Time when `print_header()` was run.

    """

    t0 = datetime.datetime.now().replace(microsecond=0)
    print("=" * 80)
    print(title)
    print("-" * 80)
    print("Command: %s" % " ".join(sys.argv))
    print("Run on host: %s" % os.uname().nodename)
    print("Operating system: %s" % os.uname().sysname)
    print("Machine: %s" % os.uname().machine)
    print("Started at: %s" % t0.strftime("%Y-%m-%d %H:%M:%S"))
    print("=" * 80)

    return t0


def print_footer(t0):
    """
    Print footer for a script

    params
        t0 : datetime object
            Time script was initialised, passed from
            `print_header`. Used to compute runtime.
    returns
        None

    """

    t1 = datetime.datetime.now().replace(microsecond=0)
    print("-" * 80)
    print("Finished at: %s" % datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    print("Time Elapsed: %s" % (t1 - t0))
    print("=" * 80)

    return None


def produce_dir(*args):
    """
    Produce a new directory by concatenating `args`,
    if it does not already exist

    params
        *args: str1, str2, str3 ...
            Comma-separated strings which will
            be combined to produce the directory,
            e.g. str1/str2/str3

    returns
        dir_name: str
            Directory name created from *args.

    """

    # Define directory path
    dir_name = os.path.join(*args)

    # Create if doesn't exist
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)

    return dir_name


# def create_logger(log_fn):

#     log = logging.getLogger(__name__)
#     log.setLevel(logging.INFO)
#     file_handler = logging.FileHandler(filename=log_fn)
#     stream_handler = logging.StreamHandler()
#     log.addHandler(file_handler)
#     log.addHandler(stream_handler)

#     return log
