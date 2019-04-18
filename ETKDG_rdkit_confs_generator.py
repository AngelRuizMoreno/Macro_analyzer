from rdkit import Chem
from rdkit.Chem import AllChem
import sys


def generate_conformers (input=None, output=None,numconfs=None):
	molecules=Chem.SDMolSupplier (sys.argv[1])
	w=Chem.SDWriter (sys.argv [2])
	for mol in molecules:
		try:
			if mol is not None:
				print ('Generating conformers for:', mol.GetProp('_Name'))
				confs = AllChem.EmbedMultipleConfs(mol, 
												  clearConfs=True, 
												  numConfs=100, 
												  pruneRmsThresh=1,
												  enforceChirality=True,
												  useExpTorsionAnglePrefs=True,  #This options disables (or enables) the ETKDG conformational algorithm (both most be True or False)
												  useBasicKnowledge=True)     #This options disables (or enables) the ETKDG conformational algorithm (both most be True or False))
				for element in confs:
					mol.SetProp ('_Name',mol.GetProp ('_Name').split('_')[0]+'_'+str(mol.GetConformers ()[element].GetId()))
					w.write (mol,confId=element)
			if mol is None:
				pass
		except Exception:
			pass
	w.close ()
############ MAIN #############

if __name__ == '__main__':
	print (''''USAGE: generate_conformers.py input.sdf output.sdf''')
	generate_conformers (input=sys.argv[1], output=sys.argv[2])
	