import pyrosetta
import argparse
import residue_functions
import Interface_Energy
import os


def make_dir(path: os.path) -> None:
    if not os.path.exists(path):
        os.makedirs(path)

def main(args):
    pose :pyrosetta.Pose= pyrosetta.pose_from_pdb(args.pdb_path)
    chain_face_1 = pyrosetta.rosetta.core.pose.get_chain_id_from_chain(str(args.face_1_chain),pose)
    chain_face_2 = pyrosetta.rosetta.core.pose.get_chain_id_from_chain(str(args.face_2_chain),pose)
    if args.face_1_res_list:

        face_1_res = str(args.face_1_res_list).split(",")
        pose_face_1_res = []
        for x in face_1_res:
            pose_face_1_res.append(int(pose.pdb_info().pdb2pose(str(args.face_1_chain),int(x))))
    else:
        pose_face_1_res = []
        for residue in range(pose.chain_begin(chain_face_1), pose.chain_end(chain_face_1)+1):
            pose_face_1_res.append(residue)

    if args.face_2_res_list:
        face_2_res = str(args.face_2_res_list).split(",")
        pose_face_2_res = []
        for x in face_2_res:
            pose_face_2_res.append(int(pose.pdb_info().pdb2pose(str(args.face_2_chain), int(x))))
    else:
        pose_face_2_res = []
        for residue in range(pose.chain_begin(chain_face_2), pose.chain_end(chain_face_2)+1):
            pose_face_2_res.append(residue)

    score_types_str = str(args.score_types).split(" ")
    score_types = []
    residues_of_interest = []
    for x in score_types_str:
        score_types.append(pyrosetta.rosetta.core.scoring.score_type_from_name(x))

    output_path: os.path = os.path.join(args.output_dir, str(args.pdb_path).split("/")[-1][0:-4])
    make_dir(output_path)
    score_fxn = pyrosetta.create_score_function(args.score_fxn)
    score_fxn.score(pose)

    face_1_energy_dict = residue_functions.calculate_residue_energy(pose, pose_face_1_res, output_path, "1")
    face_2_energy_dict = residue_functions.calculate_residue_energy(pose, pose_face_2_res, output_path, "2")
    interface_energy_dict = Interface_Energy.get_interface_energy(pose, pose_face_1_res, pose_face_2_res, output_path)


    for face_1_key in pose_face_1_res:
        if -args.score_threshold > face_1_energy_dict[face_1_key] or face_1_energy_dict[face_1_key] > args.score_threshold:
            residues_of_interest.append(face_1_key)
    for face_2_key in pose_face_2_res:
        if -args.score_threshold > face_2_energy_dict[face_2_key] or face_2_energy_dict[face_2_key] > args.score_threshold:
            residues_of_interest.append(face_2_key)
    for interface_key in interface_energy_dict:
        key_1 = int(str(interface_key).split("_")[0])
        key_2 = int(str(interface_key).split("_")[1])
        if key_1 in residues_of_interest and key_2 in residues_of_interest:
            continue
        if -args.score_threshold > interface_energy_dict[interface_key] or interface_energy_dict[interface_key] > args.score_threshold:
            if not key_1 in residues_of_interest:
                residues_of_interest.append(key_1)
            if not key_2 in residues_of_interest:
                residues_of_interest.append(key_2)

    score_type_interest_res = residue_functions.energy_by_residue_type(pose,residues_of_interest, score_types ,output_path)

    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb_path", type=str, default="", help="Path to pdb file; .../.../.../inputs/pdb_file.pdb")
    parser.add_argument("--output_dir", type=str, default="",
                        help="Path to output directory; .../.../.../project/outputs")
    parser.add_argument("--face_1_res_list", type=str, default="", help="Face 1 residue PDB index;1,2,3,4...,N")
    parser.add_argument("--face_1_chain", type=str, default="", help="Face 1 PDB chain; A")
    parser.add_argument("--face_2_res_list", type=str, default="", help="Face 2 PDB index;1,2,3,4...,N")
    parser.add_argument("--face_2_chain", type=str, default="", help="Face 2 chain; B")
    parser.add_argument("--score_fxn", type=str, default="rna_res_level_energy4.wts",
                        help="Score function name. Default; rna_res_level_energy4.wts")
    parser.add_argument("--score_threshold", type=float, default=0.5,
                        help="Absolute value of Rosetta Energy Units(REU) threhold for important residue. Default; 0.5")
    parser.add_argument("--score_types", type=str,
                        default="fa_elec_rna_phos_phos rna_torsion fa_stack rna_sugar_close hbond_sr_bb_sc",
                        help="Rosetta score function types. Default; fa_elec_rna_phos_phos rna_torsion fa_stack "
                             "rna_sugar_close hbond_sr_bb_sc")
    args = parser.parse_args()
    print(args)
    pyrosetta.init()
    main(args)
