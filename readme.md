Interface Analysis: Finds Residues within 8 Å of a target chain, scores them using the Rosetta Score Function, and output per residue energy breakdown by residue type and total energy, interface energy, and total chain energy as csv files. 

-----------------------------------------------------------------------------------------------------
**Input Flags for InterfaceAnalysis_Run.py:**:
```
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
```




**Examples:** 


>COV1_01 Model:
>>python InterfaceAnalaysis_run.py --pdb_path inputs/COV1_01.pdb --output_dir output --target_chain A 


