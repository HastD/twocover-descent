"""This is a simple multi-process manager for running two-cover descent on many curves.
The program maintains a queue of curves, starts processes for multiple curves at once,
and monitors the processes for excess memory usage or overly long runtimes."""

import resource
import time
import os

def run_curve(label):
    print("Starting {}...".format(label))
    os.system("magma -b LABEL:={} runcurve.m".format(label))
    print("Label {} complete.".format(label))
    return None

def set_memory_limit(limit):
    """Set maximum amount of memory each worker process can allocate."""
    soft, hard = resource.getrlimit(resource.RLIMIT_AS)
    resource.setrlimit(resource.RLIMIT_AS, (limit, hard))
    resource.setrlimit(resource.RLIMIT_RSS, (limit, hard))


