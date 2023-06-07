Interface Analysis: Finds Residues within 8 Å of a target chain, scores them using the Rosetta Score Function, and output per residue energy breakdown by residue type and total energy, interface energy, and total chain energy as csv files. 

Inputs: 
```
    --pdb_path: path to desired pdb file ; /PATH/TO/FILE/input.pdb; required and must be a pdb
    --output_dir: path to diretory to save generated files; /PATH/TO/OUTPUT/ ; required abd must be a pdb file
    --taget_chain: string of chain to focus; A ; required and must be a valid string found in the pdb file.
    --score_fxn: string to reference rosetta score fuction; "ref2015"; not required, defines the score fuction used in scoring models. 
    --score_types: string of valid rosetta scoring types: "fa_elec_rna_phos_phos rna_torsion fa_stack";not required, score types must be separated by a space
```
----------------------------------------------------------------------------------------------------- 
The script generates a directory containing two files:

1. chain_interface_and_chain_total_score_FILE_NAME.csv contains interface score of all chains within 8 Å of the target chain. first row shows the chain name based on the pdb,total interface score, and total score of the target chain. Second shows the values of interface scores between 
2. per_residue_score_type_energies_FILE_NAME.csv contains per residue energies based on given score types. First column shows pdb index followed by the chain number. First row shows the scoring category includiong all score types, total per residue energy from score types, and total energy of residue. 

Code Organization:
1. interface_scorer.py helper functions to calculate the interface energy between chains. Call interface energy Rosetta Method. 
2. energy_methods.py helper functions to calculate per residue energies. 
3. InterfaceAnalysis_Run.py runs interface analysis, requires the follow. 
-----------------------------------------------------------------------------------------------------
 **Example:** 
```
    >>python InterfaceAnalaysis_run.py --pdb_path inputs/COV1_01.pdb --output_dir output --target_chain A 
    >>cd output/COV_01
    >>output/COV_01 ls
    >>chain_interface_and_chain_total_score_COV1_01.csv per_residue_per_score_type_energies_COV1_02.csv
```





