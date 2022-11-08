Interface Energy: Take two list of residues alongisde its PDB chain to compare the two faces: 1) it finds residue total energy of all residues saves, it find the short range interface energy defined by rosetta and it defines critical residues based on a score threshold and it breaks down the desire score types 

**Input Flags for interface analyser:**
'
> - parser.add_argument("--pdb_path", type=str, default="", help="Path to pdb file; .../.../.../inputs/pdb_file.pdb") <br>
> - parser.add_argument("--output_dir", type=str, default="", help="Path to output directory; .../.../.../project/outputs")  <br>
> - parser.add_argument("--face_1_res_list", type=str, default="", help="Face 1 residue PDB index;1,2,3,4...,N")<br>
> - parser.add_argument("--face_1_chain", type=str, default="", help="Face 1 PDB chain; A")<br>
> - parser.add_argument("--face_2_res_list", type=str, default="", help="Face 2 PDB index;1,2,3,4...,N")<br>
> - parser.add_argument("--face_2_chain", type=str, default="", help="Face 2 chain; B")<br>
> - rser.add_argument("--score_fxn", type=str, default="rna_res_level_energy4.wts", help="Score function name. <br> Default; rna_res_level_energy4.wts")<br>
> - parser.add_argument("--score_threshold", type=float, default=0.5, help="Absolute value of Rosetta Energy Units(REU) threhold for important residue. Default; 0.5")<br>
> - parser.add_argument("--score_types", type=str, default="fa_elec_rna_phos_phos rna_torsion fa_stack rna_sugar_close hbond_sr_bb_sc", help="Rosetta score function types. Default; fa_elec_rna_phos_phos rna_torsion fa_stack rna_sugar_close hbond_sr_bb_sc")<br>

**Examples:** 


>COV1_01 Model:
>>python protein_rna_interface_analysis_run.py --pdb_path inputs/COV1_01.pdb --output_dir output --face_1_chain A --face_2_res_list   603,604,605,606,607,608,609,628,628,630,631 --face_2_chain 2 


