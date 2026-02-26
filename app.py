import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import AllChem, EnumerateStereoisomers
from stmol import showmol
import py3Dmol

st.set_page_config(page_title="Allene & Isomer Pro", layout="wide")
st.title("🧪 Full Isomer Analysis (Allenes & More)")

name = st.text_input("Enter Molecule Name (e.g., 2,3-pentadiene):", "2,3-pentadiene")

if st.button("Generate All Stereoisomers"):
    try:
        # 1. جلب البيانات وتحويلها لجزيء
        results = pcp.get_compounds(name, 'name')
        if not results:
            st.error("Molecule not found!")
        else:
            smiles = results[0].isomeric_smiles
            base_mol = Chem.MolFromSmiles(smiles)
            
            # إضافة الهيدروجينات ضروري جداً للألينات والـ Cis/Trans
            base_mol = Chem.AddHs(base_mol) 
            
            # 2. توليد كل الأيزومرات الممكنة (R/S, E/Z, Axial Chirality)
            opts = EnumerateStereoisomers.StereoEnumerationOptions(tryEmbedding=True, onlyUnassigned=False)
            isomers = list(EnumerateStereoisomers.EnumerateStereoisomers(base_mol, options=opts))
            
            st.success(f"Found {len(isomers)} Stereoisomers!")
            
            # 3. عرض كل أيزومر وتحليله
            for i, iso in enumerate(isomers):
                st.subheader(f"Isomer #{i+1}")
                
                # حساب الإحداثيات الـ 3D مع تحسين الهيكل (Optimization)
                AllChem.EmbedMolecule(iso, AllChem.ETKDG())
                AllChem.MMFFOptimizeMolecule(iso)
                
                col1, col2 = st.columns([1, 2])
                
                with col1:
                    st.write("### 🧬 Stereo Properties:")
                    # تحديث وتحديد الكيمياء الفراغية
                    Chem.AssignStereochemistry(iso, force=True, cleanIt=True)
                    
                    # أ- استخراج المراكز الكيرالية (R/S)
                    centers = Chem.FindMolChiralCenters(iso, includeUnassigned=True)
                    if centers:
                        for idx, label in centers:
                            st.info(f"**Chiral Center:** {iso.GetAtomWithIdx(idx).GetSymbol()}{idx} is **{label}**")
                    
                    # ب- استخراج الروابط المزدوجة (E/Z أو Cis/Trans)
                    has_ez = False
                    for bond in iso.GetBonds():
                        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                            stereo = bond.GetStereo()
                            if stereo in [Chem.rdchem.BondStereo.STEREOE, Chem.rdchem.BondStereo.STEREOZ]:
                                label = "E (Trans)" if stereo == Chem.rdchem.BondStereo.STEREOE else "Z (Cis)"
                                st.warning(f"**Double Bond:** {bond.GetBeginAtomIdx()}-{bond.GetEndAtomIdx()} is **{label}**")
                                has_ez = True
                    
                    if not centers and not has_ez:
                        st.write("No specific stereo-labels detected (Common in symmetric allenes).")

                with col2:
                    # عرض الـ 3D (Balls and Sticks)
                    view = py3Dmol.view(width=600, height=400)
                    mol_block = Chem.MolToMolBlock(iso)
                    view.addModel(mol_block, 'mol')
                    view.setStyle({'stick': {'radius': 0.15}, 'sphere': {'scale': 0.25}})
                    view.zoomTo()
                    showmol(view, height=400, width=600)
                st.divider()

    except Exception as e:
        st.error(f"Error analyzing: {e}")
