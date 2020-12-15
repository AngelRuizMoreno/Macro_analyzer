import argparse
import pandas as pd
from pymol import cmd
from rdkit import Chem, RDConfig
from rdkit.Chem import AllChem, TorsionFingerprints, Descriptors3D, rdMolDescriptors,rdFMCS,rdMolAlign
from rdkit.Chem.AllChem import GetBestRMS

def run_comparison (references=None,conformers=None):

	references_file =args.references
	conformers_file =args.conformers

	print (' -- -- -- -- -- LOADING REFERENCES IN {} -- -- -- -- -- '.format(references_file))
	print (' -- -- -- -- -- LOADING CONFORMERS IN {} -- -- -- -- -- '.format(conformers_file))

	lowest_rmsd =[]
	#out_macro_refs		= Chem.SDWriter('Macrocycles_references.sdf')
	out_RMSD			= Chem.SDWriter('Aligment_RMSD.sdf')
	#out_macro_confs		= Chem.SDWriter('Macrocycles_conformers.sdf')

	print (' -- -- -- -- -- STARTING ANALYSIS -- -- -- -- -- ')
	ref_index=1
  
	for ref in Chem.SDMolSupplier(references_file):

		print (ref_index,')',ref.GetProp('_Name').split('_')[0])
		
		macrocycle_atoms=[]
		all_cycles = Chem.GetSymmSSSR(ref)
		cycles_size=[i for i in all_cycles if len(i)>=8]
		for element in cycles_size:
			macrocycle_atoms+=list(element)
		all_atoms=[i.GetIdx() for i in ref.GetAtoms()]
		atoms_to_remove=(list(set(all_atoms) - set(macrocycle_atoms)))

		macrocycle = Chem.RWMol(ref)
		for i in sorted(atoms_to_remove,reverse=True):
			macrocycle.RemoveAtom (i)
	
		m_ref=macrocycle.GetMol()
		m_ref.UpdatePropertyCache()
		
		macrocycle_atoms=sorted(list(set(macrocycle_atoms)))
		print ('Initial Num Atoms:',len (all_atoms))
		print ('Macrocycle length:',len (macrocycle_atoms))
		
		m_ref_smiles=Chem.MolFragmentToSmiles(ref,macrocycle_atoms,kekuleSmiles=True)
		m_ref_smiles=Chem.MolFromSmiles(m_ref_smiles,sanitize=False)
		  
		ref_index=ref_index +1
		mol_index=0
		table=pd.DataFrame()
		
		for mol in Chem.SDMolSupplier(conformers_file):
			if ref.GetProp('_Name').split('_')[0]==mol.GetProp('_Name').split('_')[0]:
				
				table.loc[mol_index,'Conformer'] = [mol.GetProp('_Name')]

				ref_atoms=ref.GetSubstructMatch(m_ref_smiles)
				mol_atoms=mol.GetSubstructMatch(m_ref_smiles)
				amap=zip(mol_atoms,ref_atoms)
				rms_macrocycle=AllChem.GetBestRMS (mol,ref,map=[list(amap)])
				
				mol.SetProp('RMSD_macrocycle',str(rms_macrocycle)) 
				table.loc[mol_index,'RMSD_macrocycle'] = [rms_macrocycle]

				macrocycle_atoms=[]
				all_cycles = Chem.GetSymmSSSR(mol)
				cycles_size=[i for i in all_cycles if len(i)>=8]
				for element in cycles_size:
					macrocycle_atoms+=list(element)
				all_atoms=[i.GetIdx() for i in mol.GetAtoms()]
				atoms_to_remove=(list(set(all_atoms) - set(macrocycle_atoms)))

				macrocycle = Chem.RWMol(mol)
				for i in sorted(atoms_to_remove,reverse=True):
					macrocycle.RemoveAtom (i)
			
				m_mol=macrocycle.GetMol()
				m_mol.UpdatePropertyCache()

				#m_mol=Chem.MolFragmentToSmiles(mol,macrocycle_atoms,kekuleSmiles=True)
				#m_mol=Chem.MolFromSmiles(m_mol,sanitize=False)
				
				radious_macro = Descriptors3D.RadiusOfGyration (m_mol)
				table.loc[mol_index,'RoG_macrocycle']=radious_macro
				
				tt_macro=rdMolDescriptors.GetTopologicalTorsionFingerprint (m_mol)
				table.loc[mol_index,'TF_macrocycle']= [tt_macro.GetTotalVal()]
				
				r_list=Chem.TorsionFingerprints.CalculateTorsionLists (m_ref)
				r_angles=Chem.TorsionFingerprints.CalculateTorsionAngles (m_ref,r_list[0],r_list[1])
				c_list=Chem.TorsionFingerprints.CalculateTorsionLists (m_mol)
				c_angles=Chem.TorsionFingerprints.CalculateTorsionAngles (m_mol,c_list[0],c_list[1])
			
				if len(r_angles) == len(c_angles):
					torsion_macro =Chem.TorsionFingerprints.CalculateTFD (r_angles,c_angles)
					table.loc[mol_index,'TFD_macrocycle']=[torsion_macro]
				else:
					table.loc[mol_index,'TFD_macrocycle']=['NA']  
				
		
				cmd.read_molstr(Chem.MolToMolBlock(ref),'ref')
				cmd.read_molstr(Chem.MolToMolBlock(mol),'mol')
				rmsd=cmd.rms_cur('ref','mol')
				cmd.deselect()
				cmd.delete('all')
				
				mol.SetProp('RMSD_heavy_atoms',str(rmsd))
				table.loc[mol_index,'RMSD_heavy_atoms'] = [rmsd]
				
				out_RMSD.write(mol)
				
				
				radious = Descriptors3D.RadiusOfGyration (mol)
				table.loc[mol_index,'RoG_heavy_atoms']=radious

				tt=rdMolDescriptors.GetTopologicalTorsionFingerprint (mol)
				table.loc[mol_index,'TF_heavy_atoms']= [tt.GetTotalVal()]
				
				r_list=Chem.TorsionFingerprints.CalculateTorsionLists (ref)
				r_angles=Chem.TorsionFingerprints.CalculateTorsionAngles (ref,r_list[0],r_list[1])
				c_list=Chem.TorsionFingerprints.CalculateTorsionLists (mol)
				c_angles=Chem.TorsionFingerprints.CalculateTorsionAngles (mol,c_list[0],c_list[1])
			
				if len(r_angles) == len(c_angles):
					torsion =Chem.TorsionFingerprints.CalculateTFD (r_angles,c_angles)
					table.loc[mol_index,'TFD_heavy_atoms']=[torsion]
				else:
					table.loc[mol_index,'TFD_heavy_atoms']=['NA']
				
				mol_index=mol_index+1
				
		if len(table.index) > 0:
			sort =table.sort_values ('RMSD_macrocycle',ascending=True)
			sort = sort.reset_index(drop=True)
			sort.to_csv (ref.GetProp('_Name')+'.csv')
			
			sort['Nconf']=len(sort.index)
			print ('Number of conformers analyzed:',len(sort.index))
			print ('data in file:',ref.GetProp('_Name')+'.csv')
			print ('-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- \n')

			sort['Span_Rog_macrocycle'] = float(max(sort['RoG_macrocycle'])-min(sort['RoG_macrocycle']))
			sort['Span_Rog_heavy_atoms']= float(max(sort['RoG_heavy_atoms'])-min(sort['RoG_heavy_atoms']))

			lowest_rmsd.append(sort.loc[0])
		else:
			print ('No reference or conformers found in input files for {}'.format(ref.GetProp('_Name')))
			print (' ************************************ \n')

	#out_macro_refs.close()
	out_RMSD.close()
	#out_macro_confs.close()

	print ('SAVING DATA OF LOWEST RMSD OF CONFORMERS')

	summary=pd.DataFrame(lowest_rmsd)
	summary=summary.reset_index(drop=True)
	summary.to_csv ('Lowest_RMSD_Data.csv')

	print ('Lowest RMSD Data in file: Lowest_RMSD_Data.csv')
	print ('***************************************************\n')
	print ('Structures in files: Alignment_RMSD.sdf')
	print ('***************************************************\n')
	print ('CALCULATION OF {} OUT OF {} REFERENCES DONE, FILES SAVED. THANK YOU FOR USING THIS SCRIPT \n'.format(len(summary.index),len(Chem.SDMolSupplier(references_file))))

	
if  __name__ == "__main__":
	parser = argparse.ArgumentParser(description='USAGE: python Macrocycle_analizer_V4.py -r references.sdf -c conformers.sdf')
	parser.add_argument('-r','--references',default=None,dest='references',help='sdf file of the references for comparison')
	parser.add_argument('-c','--conformers',default=None,dest='conformers',help='sdf file of the conformers for comparison')
	args = parser.parse_args()
	run_comparison()