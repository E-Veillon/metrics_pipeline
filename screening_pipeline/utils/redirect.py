#!/usr/bin/python
"""Redirect unnecessary warnings output."""


import io
import os
import sys
import ctypes
import tempfile
from typing import TextIO
from contextlib import contextmanager


libc = ctypes.CDLL(None)
c_stdout = ctypes.c_void_p.in_dll(libc, "stdout")


@contextmanager
def redirect_c_stdout(stream: TextIO = None):
    """
    A context when redirecting stdout to a given stream.

    Args:
        stream (TextIO): The stream to redirect stdout to. Can be None.
    """

    # The original fd stdout points to. Usually 1 on POSIX systems.
    original_stdout_fd = sys.stdout.fileno()

    def _redirect_stdout(to_fd):
        """Redirect stdout to the given file descriptor."""
        # Flush the C-level buffer stdout
        libc.fflush(c_stdout)
        # Flush and close sys.stdout - also closes the file descriptor (fd)
        sys.stdout.close()
        # Make original_stdout_fd point to the same file as to_fd
        os.dup2(to_fd, original_stdout_fd)
        # Create a new sys.stdout that points to the redirected fd
        sys.stdout = io.TextIOWrapper(os.fdopen(original_stdout_fd, "wb"))

    # Save a copy of the original stdout fd in saved_stdout_fd
    saved_stdout_fd = os.dup(original_stdout_fd)
    try:
        # Create a temporary file and redirect stdout to it
        tfile = tempfile.TemporaryFile(mode="w+b")
        _redirect_stdout(tfile.fileno())
        # Yield to caller, then redirect stdout back to the saved fd
        yield
        _redirect_stdout(saved_stdout_fd)
        # Copy contents of temporary file to the given stream
        tfile.flush()
        tfile.seek(0, io.SEEK_SET)
        if stream is not None:
            stream.write(tfile.read())
    finally:
        tfile.close()
        os.close(saved_stdout_fd)


c_stderr = ctypes.c_void_p.in_dll(libc, "stderr")


@contextmanager
def redirect_c_stderr(stream: TextIO = None):
    """
    A context when redirecting stderr to a given stream.

    Args:
        stream (TextIO): The stream to redirect stderr to. Can be None.
    """

    # The original fd stderr points to. Usually 1 on POSIX systems.
    original_stderr_fd = sys.stderr.fileno()

    def _redirect_stderr(to_fd):
        """Redirect stderr to the given file descriptor."""
        # Flush the C-level buffer stderr
        libc.fflush(c_stderr)
        # Flush and close sys.stderr - also closes the file descriptor (fd)
        sys.stderr.close()
        # Make original_stderr_fd point to the same file as to_fd
        os.dup2(to_fd, original_stderr_fd)
        # Create a new sys.stderr that points to the redirected fd
        sys.stderr = io.TextIOWrapper(os.fdopen(original_stderr_fd, "wb"))

    # Save a copy of the original stderr fd in saved_stderr_fd
    saved_stderr_fd = os.dup(original_stderr_fd)
    try:
        # Create a temporary file and redirect stderr to it
        tfile = tempfile.TemporaryFile(mode="w+b")
        _redirect_stderr(tfile.fileno())
        # Yield to caller, then redirect stderr back to the saved fd
        yield
        _redirect_stderr(saved_stderr_fd)
        # Copy contents of temporary file to the given stream
        tfile.flush()
        tfile.seek(0, io.SEEK_SET)
        if stream is not None:
            stream.write(tfile.read())
    finally:
        tfile.close()
        os.close(saved_stderr_fd)
