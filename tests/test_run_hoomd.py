import glob
import os
import pytest
import shutil
import subprocess as sp

@pytest.fixture(
    scope="module",
    params=[
        {"--temperature": "300", "--run_time": "10", "--timestep": "0.001", "--tau": "0.01", "--omit_lj": "[]", "--energy_scale_unit": "1.0", "--distance_scale_unit": "1.0", "--nl_type": "tree"},
        {"--temperature": "500", "--run_time": "10", "--timestep": "0.001", "--tau": "0.01", "--omit_lj": "[]", "--energy_scale_unit": "1.0", "--distance_scale_unit": "1.0", "--nl_type": "tree"},
        {"--temperature": "300", "--run_time": "100", "--timestep": "0.001", "--tau": "0.01", "--omit_lj": "[]", "--energy_scale_unit": "1.0", "--distance_scale_unit": "1.0", "--nl_type": "tree"},
        {"--temperature": "300", "--run_time": "10", "--timestep": "0.01", "--tau": "0.01", "--omit_lj": "[]", "--energy_scale_unit": "1.0", "--distance_scale_unit": "1.0", "--nl_type": "tree"},
        {"--temperature": "300", "--run_time": "10", "--timestep": "0.001", "--tau": "1.0", "--omit_lj": "[]", "--energy_scale_unit": "1.0", "--distance_scale_unit": "1.0", "--nl_type": "tree"},
        {"--temperature": "300", "--run_time": "10", "--timestep": "0.001", "--tau": "0.01", "--omit_lj": "['C-Ag']", "--energy_scale_unit": "1.0", "--distance_scale_unit": "1.0", "--nl_type": "tree"},
        {"--temperature": "300", "--run_time": "10", "--timestep": "0.001", "--tau": "0.01", "--omit_lj": "[]", "--energy_scale_unit": "10.0", "--distance_scale_unit": "1.0", "--nl_type": "tree"},
        {"--temperature": "300", "--run_time": "10", "--timestep": "0.001", "--tau": "0.01", "--omit_lj": "[]", "--energy_scale_unit": "1.0", "--distance_scale_unit": "10.0", "--nl_type": "tree"},
        {"--temperature": "300", "--run_time": "10", "--timestep": "0.001", "--tau": "0.01", "--omit_lj": "[]", "--energy_scale_unit": "1.0", "--distance_scale_unit": "1.0", "--nl_type": "cell"},
        {"--temperature": "300", "--run_time": "10", "--timestep": "0.001", "--tau": "0.01", "--omit_lj": "[]", "--energy_scale_unit": "1.0", "--distance_scale_unit": "1.0", "--nl_type": "tree", "--r_cut": "5.0"},
    ],
    ids=[
        "defaults", "temperature", "run_time", "timestep", "tau", "omit_lj",
        "energy_scale_unit", "distance_scale_unit", "nl_type", "r_cut",
    ],
)
def run_hoomd(request):
    # Copy the test asset to the current directory first
    shutil.copy(
        os.path.join(os.getcwd(), "tests", "assets", "run_hoomd_test.hoomdxml"),
        os.getcwd()
    )
    cmd = ["rhaco-run-hoomd"]
    flags_dict = request.param
    for flag, val in sorted(request.param.items()):
        cmd.append(flag)
        if val is not "":
            cmd.append(val)
    cmd.append("run_hoomd_test.hoomdxml")
    run_hoomd_job = sp.Popen(cmd, stdout=sp.PIPE)
    shell_output = run_hoomd_job.communicate()
    cmd_string = ""
    for element in cmd:
        if " " in element:
            cmd_string = " ".join([cmd_string, repr(element)])
        else:
            cmd_string = " ".join([cmd_string, element])
    print(cmd_string)
    print(shell_output[0].decode("utf-8"))
    return "run_hoomd_test.hoomdxml"


class TestCreateMorphology():
    def test_check_run_hoomd(self, run_hoomd):
        morphology_name = os.path.splitext(run_hoomd)[0]
        directory = os.path.abspath(os.getcwd())
        hoomd_run_files = [
            run_hoomd,
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
