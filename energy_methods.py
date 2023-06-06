import pyrosetta

def get_chain_energy(pose: pyrosetta.Pose, chain: int) -> float:
    total_energy = 0
    eng: pyrosetta.rosetta.core.scoring.Energies = pose.energies()
    chain_start = pose.chain_begin(chain)
    chain_end = pose.chain_end(chain)
    print(chain_start, chain_end, pose.total_residue())
    for x in range(chain_start, chain_end + 1):
        total_energy += eng.residue_total_energies(x)[pyrosetta.rosetta.core.scoring.ScoreType.total_score]
    return total_energy


def residue_energies_by_residue_type_list(pose: pyrosetta.Pose, score_fxn: pyrosetta.ScoreFunction,
                                          score_type: str, rosetta_index: int) -> dict:
    working_pose: pyrosetta.Pose = pose.clone()
    score_fxn.score(working_pose)
    eng: pyrosetta.rosetta.core.scoring.Energies = working_pose.energies()
    energies_dict = {}
    score_type_list = score_type.split(" ")
    if "total_score" not in score_type_list:
        score_type_list.append("total_score")
    for score_type_name in score_type_list:
        energy_type = pyrosetta.rosetta.core.scoring.score_type_from_name(score_type_name)
        energies_dict[score_type_name] = eng.residue_total_energies(rosetta_index)[energy_type]
    return energies_dict

