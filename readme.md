Interface Analysis: Finds Residues within 8 Å of a target chain, scores them using the Rosetta Score Function, and output per residue energy breakdown by residue type and total energy, interface energy, and total chain energy as csv files. 

Inputs: 
```
    --pdb_path: path to desired pdb file ; /PATH/TO/FILE/input.pdb; required and must be a pdb
    --output_dir: path to diretory to save generated files; /PATH/TO/OUTPUT/ ; required abd must be a pdb file
    --taget_chain: string of chain to focus; A ; required and must be a valid string found in the pdb file.
    --score_fxn: string to reference rosetta score fuction; "ref2015"; not required, defines the score fuction used in scoring models. 
    --score_types: string of valid rosetta scoring types: "fa_elec_rna_phos_phos rna_torsion fa_stack";not required, score types must be separated by a space
```

The script generates two files:

1. chain_interface_and_chain_total_score_FILE_NAME.csv contains interface score of all chains within 8 Å of the target chain.  
2. per_residue_score_type_energies_FILE_NAME.csv 
-----------------------------------------------------------------------------------------------------

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


