from __future__ import annotations
import sys
import typing as tp
import argparse
import getpass

from utag.types import Path


if sys.platform.startswith("linux") or sys.platform.startswith("freebsd"):
    try:
        APP_DIR = Path("~/.utag").expanduser().mkdir()
    except Exception:
        APP_DIR = Path("./.utag").absolute().mkdir()
        print(
            f"Could not create directory '~/.utag', will use '{APP_DIR.as_posix()}'. "
            "Using a resource directory in the user's home directory improves "
            "management of resources. Make sure your home directory is writable."
        )
elif sys.platform.startswith("darwin"):
    # home = Path(f"/User/{getpass.getuser()}")
    try:
        APP_DIR = Path("/Applications/.utag").mkdir()
    except Exception:
        APP_DIR = Path("./.utag").absolute().mkdir()
        print(
            f"Could not create directory '/User/{getpass.getuser()}/.utag', "
            "will use '{APP_DIR.as_posix()}'. "
            "Using a resource directory in the user's home directory improves "
            "management of resources. Make sure your home directory is writable."
        )
#elif sys.platform.startswith("win") or sys.platform.startswith("cygwin"):
#    raise NotImplementedError("Windows support is not yet available!")
else:
    print(
        "Warning, OS could not be easily identified. Using default dir ~/.utag to store "
        "resources but that might not work!"
    )
    APP_DIR = Path("~/.utag").expanduser().mkdir()


def clear_package_data():
    """Remove static resources of package."""
    from shutil import rmtree

    rmtree(APP_DIR)
