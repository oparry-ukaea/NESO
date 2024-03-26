import glob
import hpp2plantuml
import os
import plantuml

access_strs = dict(private="-", protected="#", public="+")


def is_class_var(s):
    tmp = s.strip()
    for access_str in access_strs.values():
        if tmp.startswith(access_str) and not "(" in tmp:
            return True
    return False


def is_access(s, access_type):
    assert access_type in access_strs
    access_str = access_strs[access_type]
    tmp = s.strip()
    return tmp.startswith(access_str)


def edit_diag(
    diag_str,
    exclude_elements=[],
    funcs_only=False,
    landscape=False,
    replace_elements=[],
    skip_access=[],
    wrap_col=None,
):
    sep = "\n"
    lines = diag_str.split(sep)
    # Filter out one or more access modes
    for access in skip_access:
        lines = [line for line in lines if not is_access(line, access)]

    if funcs_only:
        lines = [line for line in lines if not is_class_var(line)]

    for exclude_element in exclude_elements:
        lines = [line for line in lines if not exclude_element in line]

    for replace_element in replace_elements:
        lines = [line.replace(replace_element[0], replace_element[1]) for line in lines]

    # orientation
    if landscape:
        lines.insert(1, "left to right direction")

    if wrap_col is not None:
        lines.insert(1, f"skinparam wrapWidth {wrap_col}")

    bg_color = "LightYellow"
    arrow_color = "Grey"
    border_color = "Black"

    class_color_block = [
        "skinparam class {",
        f"BackgroundColor {bg_color}",
        f"ArrowColor {arrow_color}",
        f"BorderColor {border_color}" + "}",
    ]
    lines[1:1] = class_color_block

    diag_str = sep.join(lines)

    return diag_str


def get_all_header_paths(solver, exclude_headers=[]):
    repo_root = os.path.dirname(os.path.dirname(get_this_dir_path()))
    solver_dir = os.path.join(repo_root, "solvers", solver)
    headers = glob.glob(solver_dir + "/**/*.hpp", recursive=True)
    for exclude_h in exclude_headers:
        headers = [h for h in headers if not exclude_h in h]
    return headers


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


def main(
    solver="H3LAPD",
    exclude_headers=[],
    out_dir=get_this_dir_path(),
    out_fname=None,
    write_int_file=False,
    **edit_kws,
):
    # Generate diagram string from header files
    header_paths = get_all_header_paths(solver, exclude_headers=exclude_headers)
    if not header_paths:
        raise RuntimeError(f"No header files found for solver '{solver}'")
    diagram_kwargs = dict()
    diag_str = gen_diag_str(header_paths, **diagram_kwargs)
    diag_str = edit_diag(diag_str, **edit_kws)

    # Write the diagram string to file (in addition to generating the png directly)
    if write_int_file:
        puml_path = os.path.join(get_this_dir_path(), f"{solver}.puml")
        write_puml(diag_str, puml_path)

    # Generate the png from the diagram file
    pl = plantuml.PlantUML("http://www.plantuml.com/plantuml/img/")
    if out_fname is None:
        out_fname = f"{solver}_eqn_systems_puml.png"
    png_path = os.path.join(out_dir, out_fname)

    try:
        png_content = pl.processes(diag_str)
        with open(png_path, "wb") as fh:
            fh.write(png_content)
        print(f"Wrote png to {png_path}")
    except plantuml.PlantUMLHTTPError as e:
        print("Failed to generate png:")
        print(e)


# Original version For M4c.4, replaced by below
# main(
#     funcs_only=True,
#     exclude_headers=["GrowthRatesRecorder", "MassRecorder"],
#     exclude_elements=["print_arr"],
#     landscape=False,
#     out_dir="/home/oparry/code/excalibur-wa/tex/pics/",
#     out_fname="H3LAPD_eqn_sys_classes.png",
#     write_int_file=True,
#     wrap_col=1000,
# )

# For NESO demo talk
main(
    funcs_only=True,
    exclude_headers=["GrowthRatesRecorder", "MassRecorder", "NeutralParticleSystem"],
    exclude_elements=[
        "print_arr",
        "create",
        "NESO.Solvers.H3LAPD.HWSystem *-- NESO.Solvers.H3LAPD.HWSystem",
    ],
    landscape=False,
    out_dir="/home/oparry/code/excalibur-wa/tex/pics",
    out_fname="H3LAPD_puml_eqn_sys.png",
    skip_access=["private"],
    write_int_file=True,
    wrap_col=800,
)

# main(
#     funcs_only=True,
#     exclude_headers=[
#         "DriftReducedSystem",
#         "HW2Din3DSystem",
#         "HW3DSystem",
#         "LAPDSystem",
#         "NeutralParticleSystem",
#     ],
#     exclude_elements=[
#         "print_arr",
#         "create",
#         "v_",
#     ],
#     landscape=False,
#     out_dir="/home/oparry/code/excalibur-wa/tex/pics",
#     out_fname="H3LAPD_puml_diagnostic_eg.png",
#     skip_access=["private"],
#     wrap_col=500,
#     write_int_file=True,
# )


# main(
#     funcs_only=True,
#     exclude_headers=[
#         "HWSystem",
#         "HW2Din3DSystem",
#         "HW3DSystem",
#         "LAPDSystem",
#         "GrowthRatesRecorder",
#         "MassRecorder",
#     ],
#     exclude_elements=[
#         "print_arr",
#         "v_",
#         "create",
#     ],
#     landscape=False,
#     out_dir="/home/oparry/code/excalibur-wa/tex/pics",
#     out_fname="H3LAPD_puml_coupling.png",
#     skip_access=["private"],
#     wrap_col=1000,
#     write_int_file=True,
# )
