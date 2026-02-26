import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import AllChem, EnumerateStereoisomers
from stmol import showmol
import py3Dmol

st.title("Chemical Isomer Analysis (Allenes Support)")
name = st.text_input("Enter Structure Name:", "2,3-pentadiene")

if st.button("Analyze Isomers"):
    try:
        results = pcp.get_compounds(name, 'name')
        if results:
            smiles = results[0].isomeric_smiles
            mol = Chem.MolFromSmiles(smiles)
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol)
            
            opts = EnumerateStereoisomers.StereoEnumerationOptions(tryEmbedding=True)
            isomers = list(EnumerateStereoisomers.EnumerateStereoisomers(mol, options=opts))
            
            st.success(f"Found {len(isomers)} potential stereoisomers.")
            
            # Display 3D
            view = py3Dmol.view(width=400, height=400)
            view.addModel(Chem.MolToMolBlock(mol), 'mol')
            view.setStyle({'stick': {}})
            view.zoomTo()
            showmol(view, height=400, width=400)
        else:
            st.error("Molecule not found.")
    except Exception as e:
        st.error(f"Error: {e}")
