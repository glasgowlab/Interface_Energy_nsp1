from typing import List, Any

import pyrosetta
import argparse
import energy_methods
import os
import pandas as pd

import interface_kyle


def make_dir(dir_path: os.path) -> None:
    if not os.path.isdir(dir_path):
        if not os.path.isfile(dir_path):
            os.makedirs(dir_path)
        else:
            print("Location not found)")


def main(args):
    pose: pyrosetta.Pose = pyrosetta.pose_from_pdb(args.pdb_path)
    target_chain_id = pyrosetta.rosetta.core.pose.get_chain_id_from_chain(str(args.target_chain), pose)
    target_and_target_neighbor_selector = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector()
    target_chain_selector = pyrosetta.rosetta.core.select.residue_selector.ChainSelector(target_chain_id)
    target_and_target_neighbor_selector.set_focus_selector(target_chain_selector)
    target_and_target_neighbor_selector.set_distance(8)
    target_and_target_neighbor_selector.set_include_focus_in_subset(True)
    target_and_target_neighbor_residues = target_and_target_neighbor_selector.apply(pose)
    residue_energies = []
    neighbor_chains = []
    residue_number_row = []
    score_fxn: pyrosetta.ScoreFunction = pyrosetta.rosetta.core.scoring.ScoreFunctionFactory.create_score_function(
        args.score_fxn)
    emopts = pyrosetta.rosetta.core.scoring.methods.EnergyMethodOptions(score_fxn.energy_method_options())
    emopts.hbond_options().decompose_bb_hb_into_pair_energies(True)
    score_fxn.set_energy_method_options(emopts)
    for i in range(1, len(target_and_target_neighbor_residues) + 1):
        res_id_pymol = pose.pdb_info().pose2pdb(i).split(" ")
        if target_and_target_neighbor_residues[i]:
            residue_number_row.append('_'.join(res_id_pymol))
            residue_energies.append(
                energy_methods.residue_energies_by_residue_type_list(pose, score_fxn, args.score_types, i))
        chain = pose.pdb_info().pose2pdb(i).split(" ")[1]
        if chain not in neighbor_chains and chain != args.target_chain:
            neighbor_chains.append(chain)
    interface_energies = []
    for c in neighbor_chains:
        interface_energies.append(interface_kyle.runScript(0, args.target_chain, c, 8, args.pdb_path))

    neighbor_chains.append("total score")
    score_fxn.score(pose)
    interface_energies.append(energy_methods.get_chain_energy(pose, pyrosetta.rosetta.core.pose.get_chain_id_from_chain
    (args.target_chain, pose)))

    data_frame = pd.DataFrame(residue_energies)
    data_frame.index = residue_number_row
    print(data_frame)
    print(neighbor_chains)
    print(interface_energies)

    file_name = args.pdb_path.split("/")[-1][0:-4]
    containing_folder = os.path.join(args.output_dir, file_name)
    make_dir(os.path.join(containing_folder))
    data_frame.to_csv(os.path.join(containing_folder, "per_residue_per_score_type_energies" + "_" + file_name + ".csv"))
    with open(
            os.path.join(containing_folder, "chain_and_total_energy" + "_" + file_name + ".csv"),'w') as writting_file:
        writting_file.write(",".join(neighbor_chains) + "\n")
        writting_file.write(",".join(interface_energies))
        writting_file.close()

    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--pdb_path", type=str, default="", help="Path to pdb file; .../.../.../inputs/pdb_file.pdb")
    parser.add_argument("--output_dir", type=str, default="",
                        help="Path to output directory; .../.../.../project/outputs")
    parser.add_argument("--target_chain", type=str, default="A", help="PDB chain of interest. Default; A")
    parser.add_argument("--score_fxn", type=str, default="rna_res_level_energy4.wts",
                        help="Score function name. Default; rna_res_level_energy4.wts")
    parser.add_argument("--score_types", type=str,
                        default="fa_elec_rna_phos_phos rna_torsion fa_stack rna_sugar_close hbond_sr_bb_sc",
                        help="Rosetta score function types. Default; fa_elec_rna_phos_phos rna_torsion fa_stack "
                             "rna_sugar_close hbond_sr_bb_sc")
    pyrosetta.init('-mute all')
    for file in os.listdir("inputs"):
        args = parser.parse_args()
        setattr(args, "pdb_path", os.path.join("inputs", file))
        setattr(args, "output_dir", "output")
        main(args)

