import glob
import hpp2plantuml
import os
import plantuml

from hpp2plantuml.hpp2plantuml import MEMBER_PROP_MAP as access_type_map
from hpp2plantuml.hpp2plantuml import LINK_TYPE_MAP as link_type_map


def is_mem_var(l, mem_var_name):
    tmp = l.strip()
    if tmp and tmp[0] in access_type_map.values():
        return tmp[1:].startswith(mem_var_name)
    return False


def is_static(l):
    tmp = l.strip()
    static_pref = "{static}"
    return tmp and tmp[1 : 1 + len(static_pref)] == static_pref


def is_access(l, access_type):
    try:
        access_str = access_type_map[access_type]
    except KeyError:
        print(f"'{access_type}' not recognised as a valid access type")
        raise
    tmp = l.strip()
    return tmp.startswith(access_str)


def is_skipped_relationship(l, rel):
    tmp = l.split()
    if len(tmp) == 3 and tmp[1] in link_type_map.values():
        return tmp[0].endswith(rel) or tmp[2].endswith(rel)
    return False


def edit_diag(
    diag_str,
    skip_access=[],
    skip_member_vars=[],
    skip_relationships=[],
    skip_static=False,
    replace_strs={},
):
    sep = "\n"
    lines = diag_str.split(sep)
    # Filter out some lines
    access_excluded = lambda l: True in [is_access(l, access) for access in skip_access]
    mem_var_excluded = lambda l: True in [
        is_mem_var(l, mem_var) for mem_var in skip_member_vars
    ]
    relationship_excluded = lambda l: True in [
        is_skipped_relationship(l, rel) for rel in skip_relationships
    ]
    static_excluded = lambda l: is_static(l) if skip_static else False
    lines = [
        replace_strs.get(l, l)
        for l in lines
        if not access_excluded(l)
        and not mem_var_excluded(l)
        and not relationship_excluded(l)
        and not static_excluded(l)
    ]
    # for access in skip_access:
    #     lines = [line for line in lines if not is_access(line, access)]
    diag_str = sep.join(lines)
    for find, rep in replace_strs.items():
        diag_str = diag_str.replace(find, rep)
    return diag_str


def get_solver_header_paths(solver):
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
    return os.path.realpath(os.path.dirname(__file__))


def write_puml(diag_str, output_path):
    with open(output_path, "wt") as fh:
        fh.write(diag_str)


def gen_base_classes_diagram(solver="Diffusion", skip_headers=[], **edit_kws):
    repo_root = os.path.dirname(os.path.dirname(get_this_dir_path()))
    p = os.path.join(repo_root, "include/nektar_interface/solver_base/*.hpp")
    header_paths = glob.glob(
        p,
        recursive=True,
    )
    # Add one solver which inherits from TimeEvo
    header_paths.append(
        os.path.join(repo_root, f"solvers/{solver}/EquationSystems/{solver}System.hpp")
    )

    # Filter out excluded headers
    header_excluded = lambda hp: True in [hp.endswith(f"{s}.hpp") for s in skip_headers]
    header_paths = [hp for hp in header_paths if not header_excluded(hp)]

    if not header_paths:
        raise RuntimeError(f"No header files found for base classes")
    gen_diagram_from_headers(header_paths, "base_classes", **edit_kws)


def gen_solver_diagram(solver="H3LAPD", skip_headers=[], **edit_kws):
    header_paths = get_solver_header_paths(solver)

    # Filter out excluded headers
    header_excluded = lambda hp: True in [hp.endswith(f"{s}.hpp") for s in skip_headers]
    header_paths = [hp for hp in header_paths if not header_excluded(hp)]

    if not header_paths:
        raise RuntimeError(f"No header files found for solver '{solver}'")
    gen_diagram_from_headers(header_paths, solver, **edit_kws)


def gen_diagram_from_headers(header_paths, output_base, **edit_kws):
    """Generate diagram string from header files"""
    diagram_kwargs = dict()
    diag_str = gen_diag_str(header_paths, **diagram_kwargs)
    diag_str = edit_diag(diag_str, **edit_kws)
    # Write the diagram string to file
    # Using the string directly with pl.processes doesn't seem to work...
    puml_path = os.path.join(get_this_dir_path(), f"{output_base}.puml")
    write_puml(diag_str, puml_path)

    gen_diagram_from_puml(puml_path, output_base)


def gen_diagram_from_puml(puml_path, output_base):
    # Generate the png from the diagram file
    pl = plantuml.PlantUML("http://www.plantuml.com/plantuml/img/")
    png_path = os.path.join(get_this_dir_path(), f"{output_base}.png")
    pl.processes_file(puml_path, outfile=png_path)
    print(f"Wrote {png_path}")


# gen_solver_diagram(skip_access=["private"])
gen_base_classes_diagram(
    skip_access=["private"],
    skip_member_vars=[
        "options",
        "zero_array_of_arrays",
        "v_DoInitialise",
        "v_DoSolve",
        "v_GenerateSummary",
        "v_InitObject",
        "cell_id_translation",
        "h5part",
        "nektar_graph_local_mapper",
    ],
    skip_headers=["empty_partsys", "particle_reader"],
    skip_relationships=[],
    skip_static=True,
    replace_strs={
        "SharedPtr": "",
        "NESO.Particles.PartSysBase *-- NESO.Particles.PartSysOptions": "NESO.Solvers.EqnSysBase o-- NESO.Particles.PartSysBase",
        "Array<OneD, Array<OneD, NekDouble>>": "NekArrArr",
        "Array<OneD, const Array<OneD, NekDouble>>": "NekArrConstArr",
        "LU::": "",
        "NESO::": "",
        "SD::": "",
        "std::": "",
        "\t\tclass PartSysOptions {\n\t\t\t+extend_halos_offset : int\n\t\t}": "",
        "NekDouble": "double",
    },
)

# gen_diagram_from_puml(
#     "/home/oparry/code/NESO/scripts/plantuml/base_classes.puml", "base_classes"
# )
