import subprocess

def test_imports():
    import rhaco
    from rhaco import generate
    from rhaco import simulate
    from rhaco import definitions
    from rhaco import modify_EAM_forcefield

def test_entry_points():
    subprocess.Popen(["rhaco-create-morph", "-h"])
    subprocess.Popen(["rhaco-run-hoomd", "-h"])
