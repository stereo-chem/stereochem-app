import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import AllChem, EnumerateStereoisomers
from stmol import showmol
import py3Dmol

st.set_page_config(page_title="Allene Expert", layout="wide")
st.title("🧪 Full Isomer Analysis (Allenes & More)")

name = st.text_input("Enter Molecule Name:", "2,3-pentadiene")

if st.button("Generate All Stereoisomers"):
    try:
        results = pcp.get_compounds(name, 'name')
        if not results:
            st.error("Molecule not found!")
        else:
            # هناخد الـ SMILES الخام ونشيل منها أي معلومات كيمياء فراغية عشان نجبر الكود يتوقعها هو
            raw_smiles = results[0].canonical_smiles
            base_mol = Chem.MolFromSmiles(raw_smiles)
            base_mol = Chem.AddHs(base_mol)
            
            # --- السحر هنا ---
            # السطر ده بيخلي RDKit تدور على أي مركز "محتمل" يكون كيرالي (ألين أو غيره)
            Chem.FindPotentialStereoBonds(base_mol)
            
            opts = EnumerateStereoisomers.StereoEnumerationOptions(
                tryEmbedding=True, 
                onlyUnassigned=False  # توليد كل الاحتمالات حتى لو مش متحددة في الاسم
            )
            isomers = list(EnumerateStereoisomers.EnumerateStereoisomers(base_mol, options=opts))
            
            st.success(f"Found {len(isomers)} Stereoisomers!")
            
            for i, iso in enumerate(isomers):
                st.subheader(f"Isomer #{i+1}")
                
                # بناء المجسم 3D
                iso = Chem.AddHs(iso)
                AllChem.EmbedMolecule(iso, AllChem.ETKDG())
                AllChem.MMFFOptimizeMolecule(iso)
                
                col1, col2 = st.columns([1, 2])
                with col1:
                    st.write("### 🧬 Analysis:")
                    Chem.AssignStereochemistry(iso, force=True, cleanIt=True)
                    
                    # استخراج الـ R/S و Axial Chirality
                    centers = Chem.FindMolChiralCenters(iso, includeUnassigned=True)
                    for idx, label in centers:
                        st.info(f"**Chiral Center/Axis:** {iso.GetAtomWithIdx(idx).GetSymbol()}{idx} is **{label}**")
                    
                    # استخراج الـ E/Z
                    for bond in iso.GetBonds():
                        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                            stereo = bond.GetStereo()
                            if stereo in [Chem.rdchem.BondStereo.STEREOE, Chem.rdchem.BondStereo.STEREOZ]:
                                label = "E (Trans)" if stereo == Chem.rdchem.BondStereo.STEREOE else "Z (Cis)"
                                st.warning(f"**
