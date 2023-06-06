import multiprocessing
import subprocess

import pyrosetta
import os
from subprocess import Popen, PIPE


def make_chunks(data: list, thread_count) -> dict:
    """
    Takes a list and splits it into parts based on the thread count
    :param data: a list that needs to be split up amongst the threads
    :param thread_count: the number of threads to use
    :return: None
    """
    threads = {}

    for x in range(0, thread_count):
        threads[x] = []

    thread = 0
    for x in range(0, len(data)):
        threads[thread].append(data[x])
        thread += 1
        if thread == thread_count:
            thread = 0
    return threads


def makeFaceFiles(pose: pyrosetta.Pose, faceA: pyrosetta.rosetta.core.select.residue_selector.ResidueSelector,
                  faceB: pyrosetta.rosetta.core.select.residue_selector.ResidueSelector, outputFaceA: os.path,
                  outputFaceB: os.path) -> None:

    with open(outputFaceA, "w") as output:
        for index, boolValue in enumerate(faceA.apply(pose)):
            if boolValue:
                resInfo = pose.pdb_info().pose2pdb(index + 1)
                output.write(resInfo.split(" ")[1] + " " + resInfo.split(" ")[0] + " _\n")

    with open(outputFaceB, "w") as output:
        for index, boolValue in enumerate(faceB.apply(pose)):
            if boolValue:
                resInfo = pose.pdb_info().pose2pdb(index + 1)
                output.write(resInfo.split(" ")[1] + " " + resInfo.split(" ")[0] + " _\n")


def move_list_file(files: list, destination_folder: os.path) -> None:
    for file in files:
        destinantion_path = os.path.join(destination_folder, str(file).split("/")[-1])
        os.replace(file, destinantion_path)


def runInterfaceScript(faceA: os.path, faceB: os.path, pdbPath: os.path, outputPath: os.path) -> float:
    pipe = subprocess.Popen(
        f"/home/dan/rosetta/main/source/bin/interface_energy.linuxgccrelease -s {pdbPath} -face1 {faceA} -face2 {faceB}",
        shell=True, stdout=subprocess.PIPE)
    text = str(pipe.communicate()[0])
    try:
        return float(
            text.split("##### TOTAL INTERFACE ENERGY: ")[-1].replace("n", "").replace("\\", "").replace('"', ""))
    except ValueError:
        return 10000


def makeResidueSelectors(pose: pyrosetta.Pose, chain1: str, chain2: str, radius: float)->[pyrosetta.rosetta.core.select.residue_selector.ResidueSelector,pyrosetta.rosetta.core.select.residue_selector.ResidueSelector]:
    chainA = pyrosetta.rosetta.core.select.residue_selector.ChainSelector(pyrosetta.rosetta.core.pose.get_chain_id_from_chain(chain1, pose))
    chainB = pyrosetta.rosetta.core.select.residue_selector.ChainSelector(pyrosetta.rosetta.core.pose.get_chain_id_from_chain(chain2, pose))

    nearChainA = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector()
    nearChainA.set_focus_selector(chainA)
    nearChainA.set_distance(radius)
    nearChainA.set_include_focus_in_subset(True)

    nearChainB = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector()
    nearChainB.set_focus_selector(chainB)
    nearChainB.set_distance(radius)
    nearChainB.set_include_focus_in_subset(True)
    interface = pyrosetta.rosetta.core.select.residue_selector.AndResidueSelector(selector1=nearChainA,
                                                                                  selector2=nearChainB)

    chainAList = pyrosetta.rosetta.utility.vector1_unsigned_long()
    chainBList = pyrosetta.rosetta.utility.vector1_unsigned_long()
    for index, boolean in enumerate(interface.apply(pose)):
        if boolean and pose.pdb_info().pose2pdb(index + 1).split(" ")[1] == chain1:
            chainAList.append(index + 1)

    for index, boolean in enumerate(interface.apply(pose)):
        if boolean and pose.pdb_info().pose2pdb(index + 1).split(" ")[1] == chain2:
            chainBList.append(index + 1)

    chainAInterfaceResidues = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector(chainAList)
    chainBInterfaceResidues = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector(chainBList)

    return [chainAInterfaceResidues,chainBInterfaceResidues]



def runScript( thread_num : int,chain1 : str, chain2: str, radius: float,pdb: os.path):

    pose: pyrosetta.Pose = pyrosetta.pose_from_pdb(pdb)
    faceAOutputPath: str = "FaceFiles/chainA_" + str(thread_num) + ".txt"
    faceBOutputPath = "FaceFiles/chainB_" + str(thread_num) + ".txt"
    residues = makeResidueSelectors(pose, chain1, chain2, radius)
    try:
        makeFaceFiles(pose, residues[0], residues[1], faceAOutputPath, faceBOutputPath)
    except:
        return +9999999
    score = runInterfaceScript(faceAOutputPath, faceBOutputPath, pdb, ":")
    return score
