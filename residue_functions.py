import pyrosetta
import os


def calculate_residue_energy(pose: pyrosetta.Pose, residues: list, output_folder: os.path, face: str) -> dict:
    list_file_name = os.path.join(output_folder,
                                  str(output_folder).split("/")[-1] + "_face_{}_per_residue_total_energies.csv".format(
                                      face))

    with open(list_file_name, "w") as energy_list:
        working_pose: pyrosetta.Pose = pose.clone()
        score_dict = {}
        pose_energies: pyrosetta.rosetta.core.scoring.Energies = working_pose.energies()
        energy_list.write("Residues,Total Residue Energy\n")
        total = 0
        for x in residues:
            residue_energy = pose_energies.residue_total_energy(x)
            total += residue_energy
            score_dict[x] = residue_energy
            energy_list.write(str(x) + "," + str(total) + "\n")

        score_dict["total"] = total
        energy_list.write(str("total") + "," + str(total))

    energy_list.close()

    return score_dict


def energy_by_residue_type(pose: pyrosetta.Pose, residues: list, score_types: list, output_folder: os.path) -> dict:
    working_pose: pyrosetta.Pose = pose.clone()
    score_dict = {}
    pose_energies: pyrosetta.rosetta.core.scoring.Energies = working_pose.energies()
    score_fxn_score_terms: pyrosetta.ScoreFunction = pyrosetta.rosetta.core.scoring.ScoreFunction()
    for score_type in score_types:
        score_fxn_score_terms.set_weight(score_type, 1)
    score_fxn_score_terms.score(working_pose)
    for residue in residues:
        residue_dict = {}
        for score_type in score_types:
            residue_dict[score_type] = pose_energies.residue_total_energies(residue)[score_type]
        score_dict[residue] = residue_dict
    for residue_key in score_dict:
        list_file_name = os.path.join(output_folder, str(output_folder).split("/")[-1] + "_" + str(working_pose.pdb_info().pose2pdb(residue_key)).split(" ")[1] + "_"+ str(working_pose.pdb_info().pose2pdb(residue_key)).split(" ")[0] + "_energy_by_score_term.csv")
        with open(list_file_name, "w") as residue_output_file:
            residue_output_file.write("score_type,score\n")
            for term_key in score_dict[residue_key]:
                input_line = str(term_key) + "," + str(score_dict[residue_key][term_key]) + "\n"
                residue_output_file.write(input_line)
            residue_output_file.close()
    return score_dict
