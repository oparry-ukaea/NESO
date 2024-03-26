import glob
import hpp2plantuml
import os
import plantuml


def is_access(str, access_type):
    access_strs = dict(private="-", protected="#", public="+")
    assert access_type in access_strs
    access_str = access_strs[access_type]
    tmp = str.strip()
    return tmp.startswith(access_str)


def edit_diag(diag_str, skip_access=[]):
    sep = "\n"
    lines = diag_str.split(sep)
    # Filter out one or more access modes
    for access in skip_access:
        lines = [line for line in lines if not is_access(line, access)]
    diag_str = sep.join(lines)
    return diag_str


def get_all_header_paths(solver):
    repo_root = os.path.dirname(os.path.dirname(get_this_dir_path()))
    solver_dir = os.path.join(repo_root, "solvers", solver)
    return glob.glob(solver_dir + "/**/*.hpp", recursive=True)


def gen_diag_str(header_paths, **diagram_kwargs):
    # Construct diagram object
    diagram_kwargs = {}
    diag = hpp2plantuml.Diagram(**diagram_kwargs)
    # Populate diagram with headers
    diag.create_from_file_list(header_paths)
    # Get diagram string
    return diag.render()


def get_this_dir_path():
    return os.path.dirname(__file__)


def write_puml(diag_str, output_path):
    with open(output_path, "wt") as fh:
        fh.write(diag_str)


def main(solver="H3LAPD", **edit_kws):
    # Generate diagram string from header files
    header_paths = get_all_header_paths(solver)
    if not header_paths:
        raise RuntimeError(f"No header files found for solver '{solver}'")
    diagram_kwargs = dict()
    diag_str = gen_diag_str(header_paths, **diagram_kwargs)
    diag_str = edit_diag(diag_str, **edit_kws)

    # Generate the png from the diagram file
    pl = plantuml.PlantUML("http://www.plantuml.com/plantuml/img/")
    png_path = os.path.join(get_this_dir_path(), f"{solver}.png")

    try:
        png_content = pl.processes(diag_str)
        with open(png_path, "wb") as fh:
            fh.write(png_content)
    except plantuml.PlantUMLHTTPError as e:
        print("Failed to generate png:")
        print(e)


main(skip_access=["private"])
