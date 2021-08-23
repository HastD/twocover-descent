#!/usr/bin/env python
"""This is a simple multi-process manager for running two-cover descent on many curves.
The program maintains a queue of curves, starts processes for multiple curves at once,
and monitors the processes for excess memory usage or overly long runtimes."""

import resource
import time

def set_memory_limit(limit):
    """Set maximum amount of memory each worker process can allocate."""
    soft, hard = resource.getrlimit(resource.RLIMIT_AS)
    resource.setrlimit(resource.RLIMIT_AS, (limit, hard))
    resource.setrlimit(resource.RLIMIT_RSS, (limit, hard))

def run_curve(label):
    print("Starting {}...".format(label))
    string = ""
    for _ in range(1024):
        string += 1024 * label
    time.sleep(1)
    print("Memory limits: {}".format(resource.getrlimit(resource.RLIMIT_AS)))
    time.sleep(1)
    print("Memory usage: {}".format(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
    time.sleep(1)
    print("Label {} complete.".format(label))
    return string[200]

