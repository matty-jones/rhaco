import glob
import os
import pytest
import subprocess as sp

@pytest.fixture(
    scope="module",
    params=[
        {},
        {"--stoichiometry": "{'V': 1}"},
        {"--dimensions": "5x5x3"},
        {"--crystal_separation": "10.0"},
        {
            "--template": "corundum.pdb", "--crystal_x": "0.4759",
            "--crystal_y": "0.4759", "--crystal_z": "1.299", "--dimensions": "10x10x2",
        },
        {
            "--template": "NiMnGa_UnitCell.pdb", "--crystal_x": "0.546",
            "--crystal_y": "0.546", "--crystal_z": "0.658", "--dimensions": "10x10x2",
            "--integrate_crystal": "", "--forcefield": "NiMnGa_FF.xml",
        },
        {"--crystal_bonds": "", "--dimensions": "3x3x1"},
        {"--z_reactor_size": "10.0"},
    ],
    ids=[
        "M1_defaults", "stoichiometry_check", "dimension_check",
        "plane_separation_check", "corundum_unit_cell_size_check",
        "NiMnGa_forcefield_check", "crystal_bond_visualisation_check", "z_size_check",
    ],
)
def create_morph(request):
    cmd = ["rhaco-create-morph"]
    # Add a couple of reactant molecules to the flags list so that we can run some
    # basic simulations and check that everything works
    flags_dict = request.param

    # Unless we're testing integrate_crystal, then we need to have some reactants in
    # there otherwise nothing will be integrated and rhaco-run-hoomd will fail.
    if "--integrate_crystal" not in flags_dict:
        flags_dict["--reactant_composition"] = "{'C2H6': 1}"
        flags_dict["--reactant_num_mol"] = "10"

    for flag, val in sorted(request.param.items()):
        cmd.append(flag)
        if val is not "":
            cmd.append(val)
    create_morph_job = sp.Popen(cmd, stdout=sp.PIPE)
    shell_output = create_morph_job.communicate()
    print(shell_output)
    output_file_name_line = shell_output[0].decode("utf-8").split("\n")[-3]
    output_file = output_file_name_line.split("XML file written to ")[1][:-1]
    return os.path.abspath(output_file)


class TestCreateMorphology():
    @pytest.mark.order1
    def test_check_xml_created(self, create_morph):
        (directory, file_name) = os.path.split(create_morph)
        files = os.listdir(directory)
        assert file_name in files, "".join(
            [
                "Expected the file ",
                str(file_name),
                " to exist in ",
                str(directory),
                ", but it doesn't.",
            ]
        )

    @pytest.mark.order2
    def test_check_run_hoomd(self, create_morph):
        run_hoomd_output = sp.Popen(
            ["rhaco-run-hoomd", "-r", "10", create_morph], stdout=sp.PIPE
        ).communicate()
        print(run_hoomd_output[0].decode("utf-8"))
        (directory, file_name) = os.path.split(create_morph)
        morphology_name = os.path.splitext(file_name)[0]
        hoomd_run_files = [
            create_morph,
            "".join([os.path.abspath(morphology_name), ".log"]),
            "".join([os.path.abspath(morphology_name), "_final.gsd"]),
            "".join([os.path.abspath(morphology_name), "_traj.gsd"]),
        ]
        for file_name in hoomd_run_files:
            assert os.path.isfile(file_name), "".join(
                [
                    "Expected the file ",
                    str(file_name),
                    " to exist, but it doesn't.",
                ]
            )
            os.remove(file_name)
