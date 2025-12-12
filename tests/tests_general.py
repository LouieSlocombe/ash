import os
import shutil

def test_import_ash():
    import ash
    assert ash is not None


def test_simple_optimization():
    from ash.modules.module_coords import Fragment
    from ash.modules.module_singlepoint import ZeroTheory
    from ash.interfaces.interface_geometric_new import geomeTRICOptimizer

    # Define coordinate string
    coords = """
    O       -1.377626260      0.000000000     -1.740199718
    H       -1.377626260      0.759337000     -1.144156718
    H       -1.377626260     -0.759337000     -1.144156718
    """
    # Defining fragment
    H2Ofragment = Fragment(coordsstring=coords, charge=0, mult=1)

    # Defining dummy theory
    zerotheorycalc = ZeroTheory()

    # Optimize with dummy theory
    result = geomeTRICOptimizer(fragment=H2Ofragment, theory=zerotheorycalc)
    # remove directory created by geomeTRIC
    if os.path.exists("geometric_tmpdir"):
        shutil.rmtree("geometric_OPTtraj.tmp")

    os.remove("ASH_Optimizer.result")
    os.remove("Fragment-currentgeo.xyz")
    os.remove("Fragment-optimized.xyz")
    os.remove("Fragment-optimized.ygg")
    os.remove("geometric_OPTtraj.log")
    os.remove("geometric_OPTtraj_optim.xyz")
    os.remove("initialxyzfiletric.xyz")
