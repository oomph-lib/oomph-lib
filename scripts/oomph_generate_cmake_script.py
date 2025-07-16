#!/usr/bin/env python3

# ------------------------------------------------------------------------------------------

from __future__ import annotations

import difflib
import re
import logging
from argparse import ArgumentParser, Namespace
from pathlib import Path
from typing import Any, Callable, List, Optional, Tuple, Union

# ------------------------------------------------------------------------------------------


class CustomFormatter(logging.Formatter):
    GREY = "\x1b[38;20m"
    YELLOW = "\x1b[33;20m"
    RED = "\x1b[31;20m"
    BOLD_RED = "\x1b[31;1m"
    RESET = "\x1b[0m"
    format = "%(asctime)s - %(levelname)s - (%(filename)s:%(lineno)d)\n%(message)s"

    FORMATS = {
        logging.DEBUG: GREY + format + RESET,
        logging.INFO: GREY + format + RESET,
        logging.WARNING: YELLOW + format + RESET,
        logging.ERROR: RED + format + RESET,
        logging.CRITICAL: BOLD_RED + format + RESET
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)


logger = logging.getLogger("")
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
ch.setFormatter(CustomFormatter())
logger.addHandler(ch)

# ------------------------------------------------------------------------------------------


# Up-to-date as of 26/11/22.
# A list of libraries built by oomph-lib. If a target defined in a Makefile.am is linking
# against a library that isn't in this list then I'm going to assume that you're creating
# it yourself.
OOMPHLIB_LIBRARIES = [
    "oomph_zlib",
    "oomph_hsl",
    "oomph_crbond_bessel",
    "oomph_triangle",
    "oomph_tetgen",
    "oomph_superlu_6.0.1",
    "oomph_metis_from_parmetis_4.0.3",
    "generic",
    "meshes",
    "poisson",
    "unsteady_heat",
    "advection_diffusion",
    "linear_wave",
    "linear_elasticity",
    "navier_stokes",
    "axisym_navier_stokes",
    "polar_navier_stokes",
    "spherical_navier_stokes",
    "fluid_interface",
    "constitutive",
    "solid",
    "axisym_spherical_solid",
    "beam",
    "shell",
    "multi_physics",
    "womersley",
    "advection_diffusion_reaction",
    "spherical_advection_diffusion",
    "axisym_advection_diffusion",
    "biharmonic",
    "flux_transport",
    "young_laplace",
    "steady_axisym_advection_diffusion",
    "mesh_smoothing",
    "helmholtz",
    "time_harmonic_linear_elasticity",
    "rigid_body",
    "time_harmonic_fourier_decomposed_linear_elasticity",
    "fourier_decomposed_helmholtz",
    "foeppl_von_karman",
    "axisym_foeppl_von_karman",
    "pml_helmholtz",
    "axisym_linear_elasticity",
    "pml_fourier_decomposed_helmholtz",
    "darcy",
    "ode",
    "poroelasticity",
    "axisym_poroelasticity",
    "pml_time_harmonic_linear_elasticity",
    "axisym_displ_based_foeppl_von_karman",
    "linearised_axisym_navier_stokes",
    "linearised_navier_stokes",
    "generalised_newtonian_navier_stokes",
    "generalised_newtonian_axisym_navier_stokes",
    "space_time_block_preconditioner",
    "space_time_navier_stokes",
    "space_time_unsteady_heat_equal_order_galerkin",
    "space_time_unsteady_heat_equal_order_galerkin_petrov",
    "space_time_unsteady_heat_mixed_order_galerkin_petrov",
]


class MakefilePatterns:
    AM_CPPFLAGS = r"AM_CPPFLAGS.+\n"
    EXTRA_DIST = r"EXTRA_DIST.+\n"
    CLEAN_AND_DISTCLEAN_LOCAL = r"\b(dist|)clean-local:(\s+)?\n(\t.+?(\n|$))*"
    DOC_DIRS = r"\bdoc_dirs:(\s+)?\n(\t.+?(\n|$))+"
    NOINST_PROGRAMS = r"noinst_PROGRAMS.*(\n|$)"
    CHECK_PROGRAMS = r"check_PROGRAMS.*(\n|$)"
    SOURCES = r"([a-zA-Z0-9_]+)_SOURCES\s*=\s*(.+)\n"
    REMAINING_SOURCES = r"(.+)_SOURCES\s*=\s*(.+)\n"
    def MATCHING_CXXFLAGS(name: str): return name + r"_CXXFLAGS\s*=\s*?(.+)?\n"
    def MATCHING_LDADD(name: str): return name + r"_LDADD\s*=\s*(.+)\n"
    def MATCHING_LDFLAGS(name: str): return name + r"_LDFLAGS\s*=\s*(.+)\n"
    SUBDIRS = r"SUBDIRS\s?=\s*?(.+?)?\n"
    INCLUDES = r"INCLUDES.*\n"
    INCLUDE_TOP_SRCDIR = r"include \$\(top_srcdir\).+\n"
    ANY_RECIPE = r"\b.+?\s*:(\s+)?(.+)?\n(\t.+?(\n|$))+"
    IF_BLOCKS = r"if\s+([^\(].+?[^\)])(\s+?(.|\s)+?)endif"
    HAVE_PYTHON_IF_BLOCK = r"if HAVE_PYTHON\s+(TESTS)?(.+)\n(else(.|\s)+?(TESTS)?(.+)\s+)?endif"


def remove_prefix(text: str, prefix: st) -> str:
    return text[len(prefix):] if text.startswith(prefix) else prefix


class MakefileToCMakeListsConverter:
    def construct_header(self, project_name: str) -> str:
        return "\n".join([
            "# ------------------------------------------------------------------------------",
            "list(APPEND CMAKE_MESSAGE_INDENT \" \")",
            f"message(VERBOSE \"Entered {project_name} subdirectory\")\n",
            "cmake_minimum_required(VERSION 3.22 FATAL_ERROR)",
            f"project({project_name} C CXX Fortran)",
            "find_package(oomphlib CONFIG REQUIRED)",
            "include(CTest)\n\n",
        ])

    def construct_footer(self, project_name: str) -> str:
        return "\n".join([
            f"\nmessage(VERBOSE \"Leaving {project_name} subdirectory\")",
            "# ------------------------------------------------------------------------------\n",
        ])

    def find_makefiles(self, base_dir: Path) -> MakefileToCMakeListsConverter:
        if not base_dir.is_dir():
            raise ValueError(f"Expected 'base_dir' argument:\n\t{base_dir}\nto be a directory!")
        return base_dir.rglob("Makefile.am")

    def raw_load(self, fpath: Path) -> str:
        text = self._load(fpath)
        return text

    def load_makefile(self, fpath: Path) -> str:
        text = self._load(fpath)
        text = self._remove_comments(text)
        text = self._remove_line_breaks(text)
        text = self._remove_tabs(text)
        text = self._squeeze_whitespaces(text)
        text = self._squeeze_newlines(text)
        return text

    def _load(self, fpath: Path) -> str:
        if (not fpath.is_file()) or (fpath.name != "Makefile.am"):
            raise ValueError(f"Expected 'fpath' argument:\n\t{fpath}\nto be a Makefile.am!")
        text: str = None
        with open(fpath, "r") as f:
            text = f.read()
        return text

    def _remove_comments(self, text: str) -> str:
        clean_text = []
        for line in text.split("\n"):
            index = line.lstrip().find("#")
            if index == -1:
                clean_text.append(line)
            elif index > 0:
                clean_text.append(line[:index])
        clean_text = "\n".join(clean_text)
        return clean_text

    def _remove_line_breaks(self, text: str) -> str:
        return text.replace('\\\n', " ")

    def _remove_tabs(self, text: str) -> str:
        return text.replace(r"\t", " ")

    def _squeeze_whitespaces(self, text: str) -> str:
        return re.sub(r" +", " ", text)

    def _squeeze_newlines(self, text: str, max_newline: int = 1) -> str:
        return re.sub(r"\n{" + str(max_newline) + ",}", "\n" * max_newline, text)

    def has_pattern(self, text: str, pattern: MakefilePatterns) -> bool:
        return bool(re.search(pattern, text))

    def find_all(self, text: str, pattern: MakefilePatterns) -> Optional[List[re.Match]]:
        return list(re.finditer(pattern, text))

    def find(self, text: str, pattern: MakefilePatterns) -> Optional[re.Match]:
        return re.search(pattern, text)

    def remove_pattern(self, text: str, pattern: MakefilePatterns, count: Optional[int] = 0) -> str:
        return re.sub(pattern, "", text, count=count)

    def replace_pattern(self, text: str, pattern: MakefilePatterns, substitution: Union[str, Callable], count: Optional[int] = 0) -> str:
        return re.sub(pattern, substitution, text, count=count)

    def convert_if_blocks(self, text: str) -> str:
        def _substitute(match: re.Match) -> str:
            (block_key, block_text) = (match.group(1), match.group(2).strip())
            if len(block_text) == 0:
                return ""
            if block_key == "SYMBOLIC_LINKS_FOR_HEADERS":
                return ""
            substitution = "\n".join([
                f"if ({block_key})",
                block_text.strip(),
                "endif()"
            ])
            return substitution
        text = self.replace_pattern(text, MakefilePatterns.IF_BLOCKS, _substitute)
        return text

    def convert_subdirs(self, text: str) -> str:
        def _substitute(match: re.Match) -> str:
            n_subdir = len(match.groups())
            if (len(match.groups()) == 0) or (match.group(1) is None):
                return ""
            subdirs = match.group(1).strip()
            substitution = "\n".join([
                f"set(SUBDIRS {subdirs})\n",
                "foreach(SUBDIR IN LISTS SUBDIRS)",
                "  add_subdirectory(${SUBDIR})",
                "endforeach()\n\n",
            ])
            return substitution
        text = self.replace_pattern(text, MakefilePatterns.SUBDIRS, _substitute)
        return text

    def convert_targets(self, text: str) -> str:
        def _convert_ldadd_libs(match: re.Match) -> Optional[str]:
            target_libs = match.group(1).strip().split()
            names_to_remove = ["-L@libdir@", "$(EXTERNAL_LIBS)", "$(FLIBS)"]
            target_libs = list(
                filter(lambda x: True if x not in names_to_remove else False, target_libs)
            )
            target_libs = [remove_prefix(x, "-l") for x in target_libs]
            target_libs_with_namespace = []
            for lib in target_libs:
                if lib in OOMPHLIB_LIBRARIES:
                    target_libs_with_namespace += [f"oomph::{lib}"]
                else:
                    target_libs_with_namespace += [lib]
                    logger.warning(
                        f"Target in Makefile.am is linking against a library '{lib}' that oomph-lib does not define. You will likely need to define it yourself."
                    )
            if "oomph::generic" not in target_libs_with_namespace:
                target_libs_with_namespace += ["oomph::generic"]
            # insert oomph::meshes before oomph::generic
            if "oomph::meshes" not in target_libs_with_namespace:
                generic_pos = target_libs_with_namespace.index("oomph::generic")
                target_libs_with_namespace = [
                    *target_libs_with_namespace[:generic_pos],
                    "oomph::meshes",
                    *target_libs_with_namespace[generic_pos:],
                ]
            if len(target_libs_with_namespace) == 0:
                return None
            target_libs_with_namespace = " ".join(target_libs_with_namespace)
            return target_libs_with_namespace

        def _convert_ldflags(match: re.Match) -> Optional[str]:
            target_name = str(match[0]).split("_LDFLAGS")[0]
            warning_message = f"""
****************************************************************************************
Linker flags have been specified for the target '{target_name}' in this Makefile.am
with _LDFLAGS. We do not currently support the conversion of the _LDFLAGS as it is a
non-trivial task. However, you should question whether arguments of the form '-l<libname>'
can be specified as libraries to link against via the 'LIBRARIES' argument to
oomph_add_executable(...).
****************************************************************************************
            """
            logger.warning(warning_message)
            return ""

        def _convert_cxxflags(match: re.Match) -> Optional[str]:
            all_flags = match.group(1)
            if all_flags is None:
                return None
            all_flags = all_flags.strip().strip("\"").split()
            target_flags = [x for x in all_flags if x.startswith("-D")]
            if len(target_flags) == 0:
                return None
            target_flags = " ".join(target_flags)
            return target_flags

        match = self.find(text, MakefilePatterns.SOURCES)
        while match is not None:
            (target_name, target_sources) = (match.group(1), match.group(2).strip())

            # Assumption: there is only one occurrence of each unique <target>_SOURCES,
            # <target>_LDADD, <target>_CXXFLAGS, <target>_LDFLAGS
            ld_match = self.find(text, MakefilePatterns.MATCHING_LDADD(target_name))
            cxxflags_match = self.find(text, MakefilePatterns.MATCHING_CXXFLAGS(target_name))
            ldflags_match = self.find(text, MakefilePatterns.MATCHING_LDFLAGS(target_name))

            (target_libs, target_cxxflags, target_ldflags) = (None, None, None)
            if ld_match is not None:
                target_libs = _convert_ldadd_libs(ld_match)
                if target_libs.strip() == "":
                    logger.warning(
                        f"Target '{target_name}' does not link to any libraries. I'll write the conversion for you but it won't work if you forgot to specify libraries to link to!")
            if ldflags_match is not None:
                ldflags_match = _convert_ldflags(ldflags_match)
            if cxxflags_match is not None:
                target_cxxflags = _convert_cxxflags(cxxflags_match)

            cmake_text = [
                f"\noomph_add_executable(",
                f"  NAME {target_name}",
                f"  SOURCES {target_sources}",
            ]
            if target_libs is not None:
                cmake_text += [f"  LIBRARIES {target_libs}"]
            if target_cxxflags is not None:
                cmake_text += [f"  CXX_DEFINITIONS {target_cxxflags}"]
            if target_libs is None:
                cmake_text += [f"  SILENCE_NO_LIBS_SUPPLIED_WARNING"]
            cmake_text = "\n".join(cmake_text) + ")\n"

            # Replace the _SOURCES target with the CMake text
            (start, end) = (match.start(), match.end())
            text = text[:start] + cmake_text + text[end:]

            # Remove the matching _LDADD and, if there was one, _CXXFLAGS line
            text = self.remove_pattern(text, MakefilePatterns.MATCHING_LDADD(target_name), count=1)
            if cxxflags_match:
                text = self.remove_pattern(
                    text, MakefilePatterns.MATCHING_CXXFLAGS(target_name), count=1)
            if ldflags_match:
                text = self.remove_pattern(
                    text, MakefilePatterns.MATCHING_LDFLAGS(target_name), count=1)

            # Move onto next match
            match = self.find(text, MakefilePatterns.SOURCES)
        return text

    def convert_recipes(self, text: str) -> str:
        def _substitute(match: re.Match) -> str:
            warning_message = """
****************************************************************************************
This Makefile.am contains a Makefile recipe but the conversion script cannot convert
recipes reliably yet. Converting recipes that just run a command is simple; you just
need to wrap the recipe in an 'add_custom_target' command. For example, given the recipe

    long_0.8:
        ./make_symbolic_links.bash prepared_mesh_files/long_0.8

you would write

    add_custom_target(long_0.8 COMMAND ./make_symbolic_links.bash prepared_mesh_files/long_0.8)

For more complex recipes where you wish to precompile some code into an object file,
you can use 'add_library(<NAME> OBJECT <SOURCES>)'. For example, the recipe

    triangle.o : triangle.c triangle.h
        $(CC) $(CFLAGS) $(AM_CFLAGS) -DANSI_DECLARATORS -DLINUX -DTRILIBRARY -MT $@ -MD -MP -MF $(DEPDIR)/$*.Tpo -c -o $@ $<
        mv -f $(DEPDIR)/$*.Tpo $(DEPDIR)/$*.Po

would be written in your CMakeLists.txt as

    add_library(triangle_lib OBJECT triangle.c triangle.h)
    target_compile_definitions(triangle_lib PRIVATE -DANSI_DECLARATORS -DLINUX -DTRILIBRARY)

****************************************************************************************
            """
            logger.warning(warning_message)
            return ""
        text = self.replace_pattern(text, MakefilePatterns.ANY_RECIPE, _substitute)
        return text

    def print_delta(self, original_text: str, new_text: str) -> None:
        print("----------------------------------------------- START OF DELTA")
        diff = difflib.ndiff(original_text.splitlines(), text.splitlines())
        delta = [l for l in diff if l.startswith('+ ') or l.startswith('- ')]
        print("\n".join(delta))
        print("-------------------------------------------------- END OF DELTA")

    def _convert_with_catch(self, fn: Callable, *args: Tuple[Any], fpath: Path) -> str:
        try:
            text = fn(*args)
        except Exception as e:
            print(f"************************************************************************")
            print(f"Failed to convert file:\n\t{fpath}")
            print(f"------------------------------------------------------------------------")
            print("ERROR MESSAGE:", str(e))
            print(f"************************************************************************")
            raise
        return text

    def save(self, outpath: Path, text: str, overwrite: bool = False) -> str:
        if (not outpath.exists()) or (overwrite):
            with open(outpath, "w") as f:
                f.write(text)
        else:
            logger.error(f"File\n\t'{outpath}'\nalready exists!")
        print(f"Saved to: {outpath}")

    def convert(self, fpath: Path, save: bool = False, overwrite: bool = False) -> str:
        if (not fpath.is_file()) or (fpath.name != "Makefile.am"):
            raise ValueError(f"Expected argument 'fpath':\n\t{fpath}\nto be a Makefile.am!")
        text = self.load_makefile(fpath)
        if re.search(r"libname\s*?=", text):
            raise NotImplementedError(
                "Unable to convert Makefile.am files that define libraries. You're going to have to do this one by hand...")
        header = self.construct_header(project_name=fpath.parent.name)
        footer = self.construct_footer(project_name=fpath.parent.name)
        if text.strip() == "":
            print("...input Makefile.am is empty.")
            if save:
                print("...creating a CMakeLists.txt file with no targets.")
                text = header + footer
                outpath = fpath.parent / "CMakeLists.txt"
                self.save(outpath, text, overwrite)
                return
        text = self.remove_pattern(text, MakefilePatterns.CHECK_PROGRAMS)
        text = self.remove_pattern(text, MakefilePatterns.NOINST_PROGRAMS)
        text = self.remove_pattern(text, MakefilePatterns.AM_CPPFLAGS)
        text = self.remove_pattern(text, MakefilePatterns.EXTRA_DIST)
        text = self.remove_pattern(text, MakefilePatterns.HAVE_PYTHON_IF_BLOCK)
        text = self.remove_pattern(text, MakefilePatterns.CLEAN_AND_DISTCLEAN_LOCAL)
        text = self.remove_pattern(text, MakefilePatterns.DOC_DIRS)
        text = self.remove_pattern(text, MakefilePatterns.INCLUDES)
        text = self.remove_pattern(text, MakefilePatterns.INCLUDE_TOP_SRCDIR)
        text = self.convert_if_blocks(text)
        text = self.convert_subdirs(text)
        text = self.convert_targets(text)
        text = self.convert_recipes(text)
        text = header + text + footer
        text = self._squeeze_newlines(text, max_newline=2)
        if save:
            outpath = fpath.parent / "CMakeLists.txt"
            self.save(outpath, text, overwrite)
        return text


def parse_args() -> Namespace:
    # fmt: off
    parser = ArgumentParser("Makefile.am to CMakeLists.txt Converter")
    mutex = parser.add_mutually_exclusive_group(required=True)
    mutex.add_argument("-d", "--directory", type=Path, required=False, help="The directory to process")
    mutex.add_argument("-f", "--file", type=Path, required=False, help="Single file to process")
    mutex = parser.add_mutually_exclusive_group(required=True)
    mutex.add_argument("-s", "--save", action="store_true", help="Write it to disk")
    mutex.add_argument("-p", "--preview", action="store_true", help="Preview the converted CMakeLists.txt file")
    parser.add_argument("-o", "--overwrite", action="store_true", help="If a previously created CMakeLists.txt file is found, overwrite it")
    args = parser.parse_args()
    # fmt: on
    return args


def main():
    args = parse_args()
    c = MakefileToCMakeListsConverter()
    files = list(c.find_makefiles(args.directory)) if args.file is None else [args.file]
    print(f"{len(files)} Makefile.am files found:")
    for file in files:
        print(f"  * {file}")
    print("")

    for file in files:
        file = file.absolute()
        print(f"Currently processing... {file}")
        try:
            cmake_text = c.convert(file, save=args.save, overwrite=args.overwrite)
        except Exception as e:
            print(f"************************************************************************")
            print(f"Failed to convert file:\n\t./{file.relative_to(Path.cwd())}")
            print(f"------------------------------------------------------------------------")
            logger.error(str(e))
            print(f"************************************************************************")
            continue
        if args.preview:
            print()
            print(">" * 90)
            print(f"{file}:")
            print("*" * (len(str(file)) + 1))
            print(c.raw_load(file))
            print("=" * 90)
            print("CMakeLists.txt:")
            print("***************")
            print(cmake_text)
            print("<" * 90)
            print("")


if __name__ == "__main__":
    main()
