import glob
import os
import pytest
import subprocess as sp

@pytest.fixture(
    scope="module",
    params=[
        {},
        # {"--stoichiometry": "{'V': 1}"},
        # {"--dimensions": "5x5x3"},
        # {"--crystal_separation": "10.0"},
        # {
        #     "--template": "corundum", "--crystal_x": "0.4759", "--crystal_y": "0.4759",
        #     "--crystal_z": "1.299", "-d": "10x10x2"
        # },
        # {"--crystal_bonds": ""},
        # {"--z_reactor_size": "10.0"},
    ],
)
def create_morph(request):
    cmd = ["rhaco-create-morph"]
    # Add a couple of reactant molecules to the flags list so that we can run some
    # basic simulations and check that everything works
    flags_dict = request.param
    flags_dict["--reactant_composition"] = "{'C2H6': 1}"
    flags_dict["--reactant_num_mol"] = "10"

    for flag, val in sorted(request.param.items()):
        cmd.append(flag)
        if val is not "":
            cmd.append(val)
    create_morph_job = sp.Popen(cmd, stdout=sp.PIPE)
    shell_output = create_morph_job.communicate()
    output_file_name_line = shell_output[0].decode("utf-8").split("\n")[-3]
    print(output_file_name_line)
    output_file = output_file_name_line.split("XML file written to ")[1][:-1]
    return os.path.abspath(output_file)


class TestCreateMorphology():
    def test_check_xml_created(self, create_morph):
        (self.directory, self.file_name) = os.path.split(create_morph)
        files = os.listdir(self.directory)
        assert self.file_name in files, "".join(
            [
                "Expected the file ",
                str(self.file_name),
                " to exist in ",
                str(self.directory),
                ", but it doesn't.",
            ]
        )

    def test_check_run_hoomd(self, create_morph):
        sp.Popen(["rhaco-run-hoomd", "-r", "10", create_morph])
        morphology_name = os.path.splitext(os.path.split(create_morph)[1])[0]
        os.remove(create_morph)
        os.remove("".join([os.path.abspath(morphology_name), ".log"]))
        os.remove("".join([os.path.abspath(morphology_name), "_final.gsd"]))
        os.remove("".join([os.path.abspath(morphology_name), "_traj.gsd"]))
