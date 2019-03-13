import glob
import os
import pytest
import subprocess as sp

@pytest.fixture(
    scope="module",
    params=[
        {"--dimensions": "3x3x1", "--z_reactor_size": "15.0", "--crystal_separation": "9.0", "--reactant_composition": "{'C2H6': 1.0}", "--reactant_num_mol": "10", "--forcefield": "FF_opls_uff"},
        {"--dimensions": "3x3x1", "--z_reactor_size": "15.0", "--crystal_separation": "9.0", "--reactant_composition": "{'C2H6': 1.0}", "--reactant_density": "0.005", "--forcefield": "FF_opls_uff"},
        {"--dimensions": "3x3x1", "--z_reactor_size": "15.0", "--crystal_separation": "9.0", "--reactant_composition": "{'perylene': 1.0}", "--reactant_density": "0.05", "--reactant_rigid": "['perylene']", "--forcefield": "FF_opls_uff_aromatic"},
        {"--dimensions": "3x3x1", "--z_reactor_size": "15.0", "--crystal_separation": "9.0", "--reactant_composition": "{'Ag_5nm': 'pos'}", "--reactant_position": "[[0.0, 0.0, 3.5]]", "--forcefield": "FF_opls_uff_aromatic"},
        {"--dimensions": "3x3x1", "--z_reactor_size": "15.0", "--crystal_separation": "9.0", "--reactant_composition": "{'C2H6': 1.0, 'perylene': 1.0}", "--reactant_density": "0.05", "--reactant_rigid": "['perylene']", "--forcefield": "FF_opls_uff_aromatic"},
        {"--dimensions": "3x3x1", "--z_reactor_size": "15.0", "--crystal_separation": "9.0", "--reactant_composition": "{'C2H6': 1.0, 'Ag_5nm': 'pos'}", "--reactant_density": "0.05", "--reactant_position": "[[0.0, 0.0, 3.5]]", "--forcefield": "FF_opls_uff_aromatic"},
        {"--dimensions": "3x3x1", "--z_reactor_size": "15.0", "--crystal_separation": "9.0", "--reactant_composition": "{'perylene': 1.0, 'Ag_5nm': 'pos'}", "--reactant_density": "0.05", "--reactant_rigid": "['perylene']", "--reactant_position": "[[0.0, 0.0, 3.5]]", "--forcefield": "FF_opls_uff_aromatic"},
        {"--dimensions": "3x3x1", "--z_reactor_size": "15.0", "--crystal_separation": "9.0", "--reactant_composition": "{'C2H6': 1.0, 'perylene': 1.0, 'Ag_5nm': 'pos'}", "--reactant_density": "0.05", "--reactant_rigid": "['perylene']", "--reactant_position": "[[0.0, 0.0, 3.5]]", "--forcefield": "FF_opls_uff_aromatic"},
        {"--dimensions": "7x3x1", "--z_reactor_size": "15.0", "--crystal_separation": "9.0", "--reactant_composition": "{'C2H6': 1.0, 'perylene': 1.0, 'Ag_5nm': 'pos'}", "--reactant_density": "0.05", "--reactant_rigid": "['perylene']", "--reactant_position": "[[-3.0, 0.0, 3.5], [3.0, 0.0, 3.5]]", "--forcefield": "FF_opls_uff_aromatic"},
        {"--template": "Ag_surface.pdb", "--dimensions": "10x10x1", "--z_reactor_size": "15.0", "--crystal_separation": "9.0", "--reactant_composition": "{'Ag_5nm': 'pos'}", "--reactant_position": "[[0.0, 0.0, 4.5]]", "--forcefield": "Ag_Zhou04.eam.alloy", "--crystal_x": "0.781", "--crystal_y": "0.781", "--crystal_z": "0.644"},
        {"--template": "corundum.pdb", "--dimensions": "15x15x1", "--z_reactor_size": "15.0", "--crystal_separation": "15.0", "--reactant_composition": "{'Ag_5nm': 'pos'}", "--reactant_position": "[[0.0, 0.0, 4.5]]", "--forcefield": "Ag_Zhou04.eam.alloy", "--crystal_x": "0.4759", "--crystal_y": "0.4759", "--crystal_z": "1.299"},
        {"--template": "corundum.pdb", "--dimensions": "15x15x1", "--z_reactor_size": "15.0", "--crystal_separation": "15.0", "--reactant_composition": "{'Ag_5nm': 'pos', 'C2H6': 1}", "--reactant_position": "[[0.0, 0.0, 4.5]]", "--reactant_density": "0.05", "--forcefield": "['Ag_Zhou04.eam.alloy', 'FF_opls_elliott.xml']", "--crystal_x": "0.4759", "--crystal_y": "0.4759", "--crystal_z": "1.299", "--omit_lj": "['Ag-Ag']"},
        #{"--template": "corundum.pdb", "--dimensions": "15x15x1", "--z_reactor_size": "15.0", "--crystal_separation": "15.0", "--reactant_composition": "{'Ag_5nm': 'pos', 'perylene': 1}", "--reactant_position": "[[0.0, 0.0, 4.5]]", "--reactant_density": "0.05", "--reactant_rigid": "['perylene']", "--forcefield": "['Ag_Zhou04.eam.alloy', 'FF_opls_elliott_aromatic.xml']", "--crystal_x": "0.4759", "--crystal_y": "0.4759", "--crystal_z": "1.299", "--omit_lj": "['Ag-Ag']"},
    ],
    ids=[
        "M1_flex_num", "M1_flex_dens", "M1_rig", "M1_pos", "M1_flex_rig",
        "M1_flex_pos", "M1_rigid_pos", "M1_flex_rig_pos", "M1_flex_rig_pos_large",
        "EAM_complete", "EAM_incomplete", "EAM_UFF_flex",
        #"EAM_UFF_rig",
    ],
)
def create_morph(request):
    cmd = ["rhaco-create-morph"]
    # Add a couple of reactant molecules to the flags list so that we can run some
    # basic simulations and check that everything works
    flags_dict = request.param

    for flag, val in sorted(request.param.items()):
        cmd.append(flag)
        if val is not "":
            cmd.append(val)
    create_morph_job = sp.Popen(cmd, stdout=sp.PIPE)
    shell_output = create_morph_job.communicate()
    cmd_string = ""
    for element in cmd:
        if " " in element:
            cmd_string = " ".join([cmd_string, repr(element)])
        else:
            cmd_string = " ".join([cmd_string, element])
    print(cmd_string)
    print(shell_output[0].decode("utf-8"))
    output_file_name_line = shell_output[0].decode("utf-8").split("\n")[-3]
    output_file = output_file_name_line.split("XML file written to ")[1][:-1]
    simulate_arguments = None
    if "--omit_lj" in flags_dict.keys():
        simulate_arguments = {"--omit_lj": flags_dict["--omit_lj"]}
    return os.path.abspath(output_file), simulate_arguments


class TestCreateMorphology():
    @pytest.mark.order1
    def test_check_xml_created(self, create_morph):
        (directory, file_name) = os.path.split(create_morph[0])
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
        run_hoomd_cmd = ["rhaco-run-hoomd", "-r", "10"]
        if create_morph[1] is not None:
            for key, val in create_morph[1].items():
                run_hoomd_cmd.append(key)
                run_hoomd_cmd.append(val)
        run_hoomd_cmd.append(create_morph[0])
        run_hoomd_output = sp.Popen(run_hoomd_cmd, stdout=sp.PIPE).communicate()
        print(run_hoomd_output[0].decode("utf-8"))
        (directory, file_name) = os.path.split(create_morph[0])
        morphology_name = os.path.splitext(file_name)[0]
        hoomd_run_files = [
            create_morph[0],
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
                    " Current directory file list: ",
                    repr(os.listdir(directory)),
                ]
            )
            os.remove(file_name)
        optional_cleanup_files = ["err.txt", "log.txt"]
        for file_name in optional_cleanup_files:
            if os.path.isfile(file_name):
                os.remove(file_name)
