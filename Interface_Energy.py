import pyrosetta
import os


# By Michael Pacella

def find_energy(energy_graph: pyrosetta.rosetta.core.scoring.EnergyGraph,
                weights: pyrosetta.rosetta.core.scoring.EMapVector, res1: int, res2: int) -> float:
    edge: pyrosetta.rosetta.core.scoring.EnergyEdge = energy_graph.find_energy_edge(res1, res2)
    if not edge:
        return 0
    edge.fill_energy_map()
    score = edge.dot(weights)
    return score


def find_long_range(resi1: int, resi2: int, pose: pyrosetta.Pose, score: pyrosetta.ScoreFunction) -> float:
    res1: pyrosetta.rosetta.core.conformation.Residue = pose.residue(resi1)
    res2: pyrosetta.rosetta.core.conformation.Residue = pose.residue(resi2)

    res1_atoms = res1.natoms()
    res2_atoms = res2.natoms()

    atom_energy = 0
    for atom1 in range(1, res1_atoms + 1):
        for atom2 in range(1, res2_atoms + 1):
            energy = pyrosetta.etable_atom_pair_energies(res1, atom1, res2, atom2, score)
            for term in energy:
                atom_energy += term
    return atom_energy


def get_interface_energy(pose: pyrosetta.Pose, face1: list, face2: list, output_path: os.path) -> dict:
    working_pose: pyrosetta.Pose = pose.clone()
    scoreLeg: pyrosetta.ScoreFunction = pyrosetta.rosetta.core.scoring.get_score_function()
    scoreLeg.score(pose)
    energies_sr: pyrosetta.rosetta.core.scoring.Energies = working_pose.energies()
    energy_graph_sr: pyrosetta.rosetta.core.scoring.EnergyGraph = energies_sr.energy_graph()
    SR_weights: pyrosetta.rosetta.core.scoring.EMapVector = energies_sr.weights()
    total_sr = 0
    energies_dict = {}
    for res_face1 in face1:
        for res_face2 in face2:
            res_key = str(res_face1) + "_" + str(res_face2)
            energy = find_energy(energy_graph_sr, SR_weights, res_face1, res_face2)
            energies_dict[res_key] = energy
            total_sr += energy


    list_file_name = os.path.join( output_path,str(output_path).split("/")[-1] +"_short_range_pair_wise_interaction_energies.csv")
    with open(list_file_name, "w") as output_file:
        output_file.write("residue_pair,pairwise_short_range_energies\n")
        for key in energies_dict:
            energy = energies_dict[key]
            line = str(key +",")
            line += str(energy)
            output_file.write(line+"\n")
        output_file.write("total,"+str(total_sr))
    output_file.close()
    return energies_dict
