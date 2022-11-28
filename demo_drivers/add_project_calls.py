#!/usr/bin/env python3

import os
import re
from pathlib import Path

demo_driver_dir = Path.cwd()


def load(path: Path) -> str:
    with open(path, "r") as f:
        text = f.read()
    return text


def filter_cmake_file(path: Path) -> bool:
    return True if (re.search("project", load(path)) is None) else False


def find_cmake_files():
    return list(filter(filter_cmake_file, demo_driver_dir.rglob("CMakeLists.txt")))


def get_text(folder, rpath):
    return f"""\
cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project({folder} C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "{rpath}")
endif()
include(CTest)

set(SUBDIRS"""


files = list(sorted(find_cmake_files()))

install_dir = demo_driver_dir.parent / "install"

for file in files:
    rpath = os.path.relpath(install_dir, file.parent)
    text_to_insert = get_text(file.parent.name, rpath)
    file_text = load(file)
    new_text = re.sub(r"set\(SUBDIRS", text_to_insert, file_text)
    with open(file, "w") as f:
        f.write(new_text)
