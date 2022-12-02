import re
from pathlib import Path

PATTERN = r"LIBRARIES((?!.*\boomph::meshes oomph::generic\b).*)"
COUNT = 0


def load(p):
    assert p.exists()
    with open(p, "r") as f:
        return f.read()


def filterer(p):
    return False if re.search(PATTERN, load(p)) is None else True


def update_file(p):
    def substitution(match) -> str:
        libs = match.group(1)
        endc = ""
        if libs.endswith(")"):
            libs = libs.rstrip(")")
            endc += ")"
        return f"LIBRARIES{libs} oomph::meshes oomph::generic{endc}"
    text = load(p)
    new_text = re.sub(PATTERN, substitution, text)
    with open(p, "w") as f:
        f.write(new_text)


files = Path.cwd().rglob("CMakeLists.txt")
files = list(filter(filterer, files))
list(map(update_file, files))
