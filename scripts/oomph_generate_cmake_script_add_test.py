#!/usr/bin/env python3

import os

ENABLE_DEBUGGING = False


def debug(func: "Callable"):
    def debugged_function(*args, **kwargs):
        if ENABLE_DEBUGGING:
            print(f"Entering function: {func.__name__}")
        result = func(*args, **kwargs)
        if ENABLE_DEBUGGING:
            print(f"Leaving function: {func.__name__}")
        return result
    return debugged_function


@debug
def get_project_cmake_script_paths(base_dir: str) -> list:
    """Walks through all of the subdirectories below 'base_dir' and appends the
    path to each CMakeLists.txt.
    """
    path_list = []
    for root, dirs, files in os.walk(base_dir, topdown=True):
        if "CMakeLists.txt" in files:
            with open(root + "/CMakeLists.txt") as f:
                file_contents = f.read()
                if "oomph_add_executable(" in file_contents:
                    path_list.append(root)
    return path_list


def already_has_test(cmake_script_path: str) -> bool:
    """Returns True if the CMake script contains "oomph_add_test"."""
    file_contents = open(cmake_script_path + "/CMakeLists.txt", "r").read()
    return file_contents.find("oomph_add_test") >= 0


def is_simple_cmake_script(cmake_script_path: str) -> bool:
    """Returns True if a CMake script contains no conditional statements, as
    conditional statements require special attention. Returns False otherwise.
    """
    file_contents = open(cmake_script_path + "/CMakeLists.txt", "r").read()

    # Remove comments
    contents_as_list = []
    for line in file_contents.splitlines():
        index = line.find("#")
        if index != -1:
            contents_as_list.append(line[:index])
        else:
            contents_as_list.append(line)
    file_contents = "\n".join(contents_as_list)

    # Determine whether the file contains a conditional statement
    is_simple = True
    if file_contents.find("if(") >= 0:
        is_simple = False

    return is_simple


@debug
def create_test_string(project_name: str,
                       executables: list,
                       equations_type: str) -> str:
    """Combine the inputs to create a test command with the appropriate inputs.
    """
    test_string = f"oomph_add_test(\n"
    test_string += f"  TEST_NAME {equations_type}.{project_name}\n"
    test_string += f"  DEPENDS_ON {' '.join(executables)}\n"
    test_string += f"  TEST_FILES validate.sh validata\n"
    test_string += f"  LABELS {equations_type} {project_name}"

    # Add any additional labels we might want
    if "one_d_" in project_name:
        test_string += f" one_d_{project_name}"
    if "two_d_" in project_name:
        test_string += f" two_d_{project_name}"
    if "three_d_" in project_name:
        test_string += f" three_d_{project_name}"
    if "adapt" in project_name:
        test_string += f" adapt"
    test_string += ")\n\n"

    return test_string


@debug
def get_equations_type(cmake_script_path: str) -> str:
    """Returns the name of the immediate subfolder of demo_drivers that contains
    this CMakeLists.txt file, which is the equation type.
    """
    return cmake_script_path.split("demo_drivers/")[1].split("/")[0]


@debug
def get_list_of_targets_in_cmake_script(cmake_script_path: str):
    """Returns all names following the TARGET keyword in a CMake script."""
    with open(cmake_script_path + "/CMakeLists.txt", "r") as f:
        file_contents_as_list = f.read().splitlines()

    executables_list = []
    for line in file_contents_as_list:
        if "TARGET " in line:
            target_args = line.strip().split(" ")[1:]
            if len(target_args) > 1:
                raise ValueError(f"""
                Expected a single argument to the TARGET keyword, but file at:
                    {cmake_script_path}
                contains the {len(target_args)} arguments: {line}
                """)
            executables_list.append(target_args[0])
    return executables_list


@debug
def get_header_and_footer_comment_line_nums(cmake_script_path: str) -> list:
    # Read the file
    with open(cmake_script_path + "/CMakeLists.txt", "r") as f:
        file_contents_as_list = f.read().splitlines()

    line_nums = []
    for i, line in enumerate(file_contents_as_list):
        if "# --------" in line:
            line_nums.append(i)
    return line_nums


@debug
def create_cmake_script_with_test(cmake_script_path: str) -> str:
    """Process a CMake project script and add a add_test(...) command to it if
    it doesn't already contain one.

    Canonical example of what we want to add:

    oomph_add_test(
        TEST_NAME one_d_poisson
        DEPENDS_ON one_d_poisson
        TEST_FILES validate.sh validata
        LABELS poisson one_d_poisson)
    """

    project_name = os.path.basename(cmake_script_path)
    executables_list = get_list_of_targets_in_cmake_script(cmake_script_path)
    equations_type = get_equations_type(cmake_script_path)
    test_string = create_test_string(
        project_name, executables_list, equations_type)
    hf_line_nums = get_header_and_footer_comment_line_nums(cmake_script_path)

    # Sanity check
    if len(hf_line_nums) != 2:
        raise ValueError(f"""
            Unexpected number of header/footer comments in file:
                {cmake_script_path}
            """)

    # Convert both multiline text strings into lists of strings
    with open(cmake_script_path + "/CMakeLists.txt", "r") as f:
        file_contents_as_list = f.read().splitlines()
    test_string_as_list = test_string.splitlines()

    # Insert the test string just before the footer comment line
    cmake_script_content = "\n".join(file_contents_as_list[:hf_line_nums[1]] +
                                     test_string_as_list[:] +
                                     file_contents_as_list[hf_line_nums[1]:])

    # print("/*******************************************************************/")
    # print(cmake_script_content)
    # print("/*******************************************************************/")
    return cmake_script_content


def main():
    base_dir = os.getcwd()
    if (os.path.basename(base_dir) != "demo_drivers"):
        raise ValueError(
            "Script needs to be run from the demo_drivers directory!")

    # Get the paths to all the CMake scripts we need to edit
    cmake_script_paths = get_project_cmake_script_paths(base_dir)

    # Differentiate between simple scripts and scripts with conditional statements
    simple_makefiles = [
        p for p in cmake_script_paths if is_simple_cmake_script(p)]
    complex_makefiles = [
        p for p in cmake_script_paths if not is_simple_cmake_script(p)]

    # Get the CMake scripts we haven't edited yet
    simple_makefiles = [
        p for p in cmake_script_paths if not already_has_test(p)]

    print("--------------- Handling unedited simple makefiles ----------------")

    for cmake_script_path in simple_makefiles:
        cmake_script_content = create_cmake_script_with_test(cmake_script_path)
        print(f"DIRECTORY: {cmake_script_path}")

        # Output it...
        # cmake_file_path = cmake_script_path + "/CMakeLists.txt"
        # with open(cmake_file_path, "w") as f:
        #     print(cmake_script_content, file=f)

    print("\n\n")
    print("------------ List of makefiles I can't handle right now -----------")
    for cmake_script_path in complex_makefiles:
        rel_path = os.path.relpath(cmake_script_path, base_dir)
        if rel_path != ".":
            print(rel_path)


if __name__ == "__main__":
    main()
