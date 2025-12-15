from ash import *
import os

def test_openmm_md():
    # Read raw PDB-file, fix using pdbfixer, setup using Modeller and run MD
    pdb_file = f"./tests/pdbfiles/1aki.pdb"

    # Setting up new system, adding hydrogens, solvent, ions and defining forcefield, topology
    omm, fragment = OpenMM_Modeller(pdbfile=pdb_file,
                                    forcefield='Amber14',
                                    watermodel="tip3p-fb",
                                    pH=7.0,
                                    solvent_padding=10.0,
                                    ionicstrength=0.1,
                                    platform='CPU')


    # openmmobject = OpenMMTheory(platform=platform, forcefield=forcefield_obj, topoforce=True,
    #                                 topology=modeller.topology, pdbfile=None, periodic=periodic,
    #                                 autoconstraints='None', rigidwater=False, printlevel=0,
    #                                 residuetemplate_choice=residuetemplate_choice)


    OpenMM_MD(fragment=fragment,
              theory=omm,
              simulation_steps=100,
              traj_frequency=2,
              trajectory_file_option='PDB',
              integrator='RPMDIntegrator',
              )