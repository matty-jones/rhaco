import glob
import os
import pytest
import subprocess as sp

@pytest.fixture(
    scope="module",
    params=[
        #{"--dimensions": "3x3x1", "--z_reactor_size": "15.0", "--crystal_separation": "9.0", "--reactant_composition": "{'C2H6': 1.0}", "--reactant_num_mol": "10", "--forcefield": "FF_opls_uff"},
        #{"--dimensions": "3x3x1", "--z_reactor_size": "15.0", "--crystal_separation": "9.0", "--reactant_composition": "{'C2H6': 1.0}", "--reactant_density": "0.005", "--forcefield": "FF_opls_uff"}
        #{"--dimensions": "3x3x1", "--z_reactor_size": "15.0", "--crystal_separation": "9.0", "--reactant_composition": "{'perylene': 1.0}", "--reactant_density": "0.05", "--reactant_rigid": "['perylene']", "--forcefield": "FF_opls_uff_aromatic"}
        #{"--dimensions": "3x3x1", "--z_reactor_size": "15.0", "--crystal_separation": "9.0", "--reactant_composition": "{'Ag_5nm': 'pos'}", "--reactant_position": "[[0.0, 0.0, 3.5]]", "--forcefield": "FF_opls_uff_aromatic"}
        #{"--dimensions": "3x3x1", "--z_reactor_size": "15.0", "--crystal_separation": "9.0", "--reactant_composition": "{'C2H6': 1.0, 'perylene': 1.0}", "--reactant_density": "0.05", "--reactant_rigid": "['perylene']", "--forcefield": "FF_opls_uff_aromatic"}
        #{"--dimensions": "3x3x1", "--z_reactor_size": "15.0", "--crystal_separation": "9.0", "--reactant_composition": "{'C2H6': 1.0, 'Ag_5nm': 'pos'}", "--reactant_density": "0.05", "--reactant_position": "[[0.0, 0.0, 3.5]]", "--forcefield": "FF_opls_uff_aromatic"}
        #{"--dimensions": "3x3x1", "--z_reactor_size": "15.0", "--crystal_separation": "9.0", "--reactant_composition": "{'perylene': 1.0, 'Ag_5nm': 'pos'}", "--reactant_density": "0.05", "--reactant_rigid": "['perylene']", "--reactant_position": "[[0.0, 0.0, 3.5]]", "--forcefield": "FF_opls_uff_aromatic"}
        #{"--dimensions": "3x3x1", "--z_reactor_size": "15.0", "--crystal_separation": "9.0", "--reactant_composition": "{'C2H6': 1.0, 'perylene': 1.0, 'Ag_5nm': 'pos'}", "--reactant_density": "0.05", "--reactant_rigid": "['perylene']", "--reactant_position": "[[0.0, 0.0, 3.5]]", "--forcefield": "FF_opls_uff_aromatic"}
        #{"--dimensions": "7x3x1", "--z_reactor_size": "15.0", "--crystal_separation": "9.0", "--reactant_composition": "{'C2H6': 1.0, 'perylene': 1.0, 'Ag_5nm': 'pos'}", "--reactant_density": "0.05", "--reactant_rigid": "['perylene']", "--reactant_position": "[[-3.0, 0.0, 3.5], [3.0, 0.0, 3.5]]", "--forcefield": "FF_opls_uff_aromatic"}
        {"--template": "Ag_surface.pdb", "--dimensions": "10x10x1", "--z_reactor_size": "15.0", "--crystal_separation": "9.0", "--reactant_composition": "{'Ag_5nm': 'pos'}", "--reactant_position": "[[0.0, 0.0, 4.5]]", "--forcefield": "Ag_Zhou04.eam.alloy", "--crystal_x": "0.781", "--crystal_y": "0.781", "--crystal_z": "0.644"}
    ],
    # ids=[
    #     "M1_flex_num", "M1_flex_dens", "M1_rig", "M1_pos", "M1_flex_rig",
    #     "M1_flex_pos", "M1_rigid_pos", "M1_flex_rig_pos", "M1_flex_rig_pos_large",
    #"EAM_complete",
    #     "EAM_incomplete", "EAM_UFF_flex", "EAM_UFF_rig",
    # ],
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
    print(cmd)
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
            #os.remove(file_name)
