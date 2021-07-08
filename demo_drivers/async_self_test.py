#!/usr/bin/env python3

import asyncio, os, re, shutil, subprocess

async def async_run_self_test(cmd):
    # See https://docs.python.org/3/library/asyncio-subprocess.html#asyncio-subprocess
    proc = await asyncio.create_subprocess_shell(
        cmd,
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE)

    stdout, stderr = await proc.communicate()

    print(f'[{cmd!r} exited with {proc.returncode}]')
    if stdout:
        print(f'[stdout]\n{stdout.decode()}')
    if stderr:
        print(f'[stderr]\n{stderr.decode()}')
    return None


def run_self_test(path: str, just_build=False) -> None:
    os.chdir(path)

    build_dir_path = os.path.join(path, "build")

    # Delete the build directory if it exists and create a new one
    if os.path.isdir(build_dir_path):
        shutil.rmtree(build_dir_path)
    os.mkdir("build")
    os.chdir("build")

    print(f"Current directory: {os.getcwd()}")
    subprocess.run(["cmake -G Ninja .."], shell=True)
    subprocess.run(["ninja"], shell=True)
    print(f"\n==> DIRECTORY: {os.getcwd()}\n\n")


    return None


def clean_up(path: str) -> None:
    build_dir_path = os.path.join(path, "build")
    os.chdir(path)
    if os.path.isdir(build_dir_path):
        shutil.rmtree(build_dir_path)



def remove_comments(file_contents: str) -> str:
    # Remove empty lines and strip leading and trailing whitespaces
    file_contents = [line.strip() for line in file_contents.split('\n') if line.strip() != ""]

    # Remove comments
    contents_as_list = []
    for line in file_contents:
        index = line.find("#")
        if index == -1:
            # No comments => take whole line
            contents_as_list.append(line)
        elif index > 0:
            # The whole line isn't a comment so take the bit that isn't
            contents_as_list.append(line[:index])

    # Return the result
    return "\n".join(contents_as_list)


def is_project_script(path: str) -> bool:
    """If the CMakeLists.txt file contains a project(...) call then it
    constitutes a project that we can build."""
    file_contents = open(path + "/CMakeLists.txt", "r").read()
    file_contents = remove_comments(file_contents)
    found_project_command = (file_contents.find("project") >= 0)
    return found_project_command


def get_cmakelists_paths(base_dir: str) -> list:
    """Walks through all of the subdirectories below base_dir and appends the
    path to each CMakeLists.txt file."""
    path_list = []
    for root, dirs, files in os.walk(base_dir, topdown=True):
        if "CMakeLists.txt" in files:
            path_list.append(root)
    return path_list


def sorting_hook_fcn(paths: list) -> list:
    return paths

def main():
    base_dir = "/Users/PuneetMatharu/Dropbox/programming/cmake/oomph-lib/oomph-lib-cmake/demo_drivers"
    # base_dir = os.getcwd()
    if (os.path.basename(base_dir) != "demo_drivers"):
        raise ValueError("Script needs to be run from the demo_drivers directory!")

    all_cmakelists_paths = get_cmakelists_paths(base_dir)

    paths_to_projects_to_build = [x for x in all_cmakelists_paths if is_project_script(x)]
    paths_to_projects_to_build = sorting_hook_fcn(paths_to_projects_to_build)
    for directory in paths_to_projects_to_build[:5]:
        # run_self_test(directory, just_build=True)
        clean_up(directory)



# async def main():
#     print('Hello ...')
#     await asyncio.sleep(1)
#     print('... World!')

if __name__ == "__main__":
    # asyncio.run(main())
    main()