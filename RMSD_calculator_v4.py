import argparse
import multiprocessing
from multiprocessing import Pool
import pandas as pd
from rdkit import Chem, RDConfig
from rdkit.Chem import AllChem, rdMolAlign, TorsionFingerprints, Descriptors3D, rdFMCS
from rdkit.Chem.Pharm2D import Gobbi_Pharm2D, Generate
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols


def run_comparison (references=None,conformers=None):

	references 	=args.references
	conformers 	=args.conformers

	templates 	=[]
	lowest_rmsd =[]
	
	for reference in AllChem.SDMolSupplier(references):
		if reference.HasProp('_Name'):
			ref_id=reference.GetProp('_Name').split('_')[0]
			templates.append ([ref_id,reference])

	
	mol_RMSD 		=[]
	mol_references 	=[]
	mol_O3A 		=[]
	mol_minimized 	=[]
	
	for refer in templates:

		try:

			print ('Processing:',refer[0])
			conformer 		=[]
			rmsd 			=[]
			similarity_3D	=[]
			O3A_result 		=[]
			t_angles 		=[]
			r_gyration 		=[]
			i_energy 		=[]
			f_energy 		=[]
			rmsd_minimized 	=[]
		

			for mol in AllChem.SDMolSupplier(conformers):

			

				if refer[0] == mol.GetProp('_Name').split('_')[0]:
				
					mol_copy=mol
					name=str (mol.GetProp('_Name'))
					conformer.append (name)

					#Aligment and RMSD calculation based on Maximum Common Structure SMARTS 
					r=rdFMCS.FindMCS([mol,refer[1]])
					a=refer[1].GetSubstructMatch(Chem.MolFromSmarts(r.smartsString))
					b=mol.GetSubstructMatch(Chem.MolFromSmarts(r.smartsString))
					mapa=list(zip(b,a))

					rms=rdMolAlign.AlignMol (mol,refer[1],atomMap=mapa)
					rmsd.append (rms)
					mol.SetProp('RMSD', str(rms))
					mol_RMSD.append (mol)
					mol_references.append (refer[1])

					# Tortional fingerprint 

					r_list	 	=	Chem.TorsionFingerprints.CalculateTorsionLists (refer[1])
					r_angles	=	Chem.TorsionFingerprints.CalculateTorsionAngles (refer[1],r_list[0],r_list[1])
					c_list		=	Chem.TorsionFingerprints.CalculateTorsionLists (mol)
					c_angles	=	Chem.TorsionFingerprints.CalculateTorsionAngles (mol,c_list[0],c_list[1])
					torsion 	=	Chem.TorsionFingerprints.CalculateTFD (r_angles,c_angles)
					t_angles.append (torsion)


					#Radious of gyration

					radious = Descriptors3D.RadiusOfGyration (mol)
					r_gyration.append (radious)
					mp = AllChem.MMFFGetMoleculeProperties(mol)
					mmff = AllChem.MMFFGetMoleculeForceField(mol, mp)
					energy_value = mmff.CalcEnergy()
					i_energy.append (energy_value)

					# Energy and minimization

					m2=mol
					AllChem.EmbedMolecule(m2)
					AllChem.MMFFOptimizeMolecule(m2,mmffVariant='MMFF94')
					mp=AllChem.MMFFGetMoleculeProperties(m2)
					mmff=AllChem.MMFFGetMoleculeForceField(m2, mp)
					energy_value_minimized=mmff.CalcEnergy()
					f_energy.append (energy_value_minimized)

					m3=Chem.RemoveHs (m2)
					r=rdFMCS.FindMCS([m3,refer[1]])
					a=refer[1].GetSubstructMatch(Chem.MolFromSmarts(r.smartsString))
					b=m3.GetSubstructMatch(Chem.MolFromSmarts(r.smartsString))
					mapa=list(zip(b,a))

					rms_2=rdMolAlign.AlignMol (m3,refer[1],atomMap=mapa)
					rmsd_minimized.append (rms_2)
					m3.SetProp('RMSD', str(rms_2))
					mol_minimized.append (m3)


					O3A=rdMolAlign.GetO3A (mol_copy,refer[1])
					align=O3A.Align()
					O3A_result.append (align)
					mol_copy.SetProp('O3A',str(align))
					mol_O3A.append (mol_copy)


			d={ 'conformer':pd.Series(conformer),
				'RMSD':pd.Series(rmsd),
				'O3A_value':pd.Series (O3A_result),
				'Torsional_Fingerprint':pd.Series (t_angles),
				'Radius_of_Gyration':pd.Series (r_gyration),
				'Initial_Energy': pd.Series (i_energy),
				'Minimization_Energy':pd.Series (f_energy),
				'RMSD_after_minimization':pd.Series(rmsd_minimized)}


			table=pd.DataFrame (d)
			sort =table.sort_values ('RMSD',ascending=True)
			sort = sort.reset_index(drop=True)
			sort.to_csv (refer[0]+'.csv')
			print ('data in file:',refer[0]+'.csv')
			print ('-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- ')

			rog_diff= (float(max(sort['Radius_of_Gyration'])-sort['Radius_of_Gyration'][0]))
			
			lowest_rmsd.append ((sort['conformer'][0],sort['RMSD'][0],sort['O3A_value'][0],
				sort['Torsional_Fingerprint'][0],sort['Radius_of_Gyration'][0],sort['Initial_Energy'][0],
				sort['Minimization_Energy'][0],sort['RMSD_after_minimization'][0],rog_diff))

		except Exception:
			print ('Something wrong with this reference or conformer')
			print ('Omitting')
			pass
		

	print ('SAVING DATA OF LOWEST RMSD OF CONFORMERS ... ... ... ... ... ... ... ...')
	summary=pd.DataFrame(data=lowest_rmsd,columns=['Conformer', 'RMSD', 'O3A_value','Torsional_Fingerprint','Radius_of_Gyration','Initial_Energy','Minimization Energy','RMSD_after_minimization','Dif_Radious_of_Gyration'])
	summary.to_csv ('Lowest_RMSD_Data.csv')
	print ('Lowest RMSD Data in file: Lowest_RMSD_Data.csv')
	print ('***************************************************')

	print ('SAVING STRUCTURES (RMSD, O3A, and MINIMIZATION) ... ... ... ... ... ... ... ... ...')
	output_Ref	= Chem.SDWriter('Aligned_Refrences.sdf')
	output_RMSD	= Chem.SDWriter('RMSD_alignment.sdf')
	output_O3A	= Chem.SDWriter('O3A_alignment.sdf')
	output_Min	= Chem.SDWriter('Minimization.sdf')
	
	mol_references = list(set(mol_references))
	  
	[output_Ref.write (element) for element in mol_references]
	output_Ref.close ()
	[output_RMSD.write (element) for element in mol_RMSD]
	output_RMSD.close()
	[output_O3A.write (element) for element in mol_O3A]
	output_O3A.close ()
	[output_Min.write (element) for element in mol_minimized]
	output_Min.close ()


	print ('Structures in files: Aligned_Refrences.sdf, RMSD_alignment.sdf, O3A_alignment.sdf, and Minimization.sdf ')

	print ('ALL THE CALCULATIONS DONE, FILES SAVED. THANK YOU FOR USING THIS SCRIPT')

if  __name__ == "__main__":
	parser = argparse.ArgumentParser(description='example usage:  python RMSD_calculator.py -r references.sdf -c conformers.sdf')
	parser.add_argument('-r','--references',default=None,dest='references',help='sdf file of the references for comparison')
	parser.add_argument('-c','--conformers',default=None,dest='conformers',help='sdf file of the conformers for comparison')
	args = parser.parse_args()
	run_comparison()