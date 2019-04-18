from rdkit import Chem
from rdkit.Chem import AllChem
import sys

if __name__ == "__main__":

	try:
		print ('''
			*** USAGE: python Macrocycle_extractor.py input_file.sdf output_file.sdf ***
			''')
		file=sys.argv[1]
		out_file= sys.argv[2]

		chem=Chem.SDMolSupplier(file)
		to_file=[]
		for molecule in chem:
			
			found_macrocycles=[]
			if molecule !=None:
				molecule=Chem.AddHs(molecule)
				Chem.Kekulize (molecule, clearAromaticFlags=True)
				all_cycles = Chem.GetSymmSSSR(molecule)

				for cycle in all_cycles:
					if len (cycle) >= 8:
						found_macrocycles.append (cycle)

				longest_macro=[]
				for macro in found_macrocycles:
					longest_macro.append (len(macro))
				if len(macro) == max (longest_macro):
					macrocycle = macro
					macrocycle_final_structure = Chem.RWMol(molecule)


				all_atoms=[]        
				atoms_to_remove=[]
				
				print (molecule.GetProp ('_Name'))
				macrocycle_atoms=list(macrocycle)
				all_atoms+=[atom.GetIdx() for atom in molecule.GetAtoms()]
				atoms_to_remove=(list(set(all_atoms) - set(macrocycle_atoms)))
				print ('Initial Num Atoms:',macrocycle_final_structure.GetNumAtoms())
				print ('Macrocycle length:',len (macrocycle_atoms))
				
				for i in sorted(atoms_to_remove,reverse=True):
					if i < macrocycle_final_structure.GetNumAtoms ():
						macrocycle_final_structure.RemoveAtom (i)
				print ('Final Num Atoms:',macrocycle_final_structure.GetNumAtoms())
				to_file.append (macrocycle_final_structure)
				
				print ('-- -- -- -- -- -- -- -- -- -- --')

		output = Chem.SDWriter(out_file)     
		for element in to_file:
			Chem.RemoveStereochemistry (element)
			#Chem.SanitizeMol(element,sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL^Chem.SANITIZE_ADJUSTHS^Chem.SANITIZE_SETAROMATICITY^Chem.SANITIZE_KEKULIZE)
			output.write (element)
		output.close ()
	except Exception:
		pass