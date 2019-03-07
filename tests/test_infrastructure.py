import subprocess as sp

def test_imports():
    import rhaco
    from rhaco import generate
    from rhaco import simulate
    from rhaco import definitions
    from rhaco import modify_EAM_forcefield

def test_entry_points():
    sp.Popen(["rhaco-create-morph", "-h"])
    sp.Popen(["rhaco-run-hoomd", "-h"])
