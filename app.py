import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import AllChem, EnumerateStereoisomers
from stmol import showmol
import py3Dmol

st.set_page_config(page_title="Chemical Isomer Explorer", layout="wide")
st.title("🧪 Comprehensive Isomer Analysis")
st.markdown("Analysis of **R/S**, **E/Z**, and **Axial Chirality** (Allenes)")

name = st.text_input("Enter Molecule Name (e.g., 2,3-pentadiene or 2-butene):", "2,3-pentadiene")

if st.button("Generate & Analyze All Isomers"):
    try:
        # 1. جلب البيانات من PubChem
        results = pcp.get_compounds(name, 'name')
        if not results:
            st.error("Molecule not found!")
        else:
            smiles = results[0].isomeric_smiles
            # تحويل SMILES إلى جزيء مع إضافة الهيدروجينات للتحليل الهندسي
            base_mol = Chem.MolFromSmiles(smiles)
            base_mol = Chem.AddHs(base_mol)
            
            # 2. إعدادات التوليد الشامل (تشمل الروابط المزدوجة والمراكز الكيرالية)
            opts = EnumerateStereoisomers.StereoEnumerationOptions(
                tryEmbedding=True, 
                onlyUnassigned=False  # توليد كل الاحتمالات حتى لو كانت محددة
            )
            isomers = list(EnumerateStereoisomers.EnumerateStereoisomers(base_mol, options=opts))
            
            st.success(f"Found {len(isomers)} possible stereoisomers (including E/Z and R/S)")

            for i, iso in enumerate(isomers):
                # حساب الإحداثيات الثلاثية لكل أيزومر بشكل منفصل
                iso = Chem.AddHs(iso)
                AllChem.EmbedMolecule(iso, AllChem.ETKDG())
                AllChem.MMFFOptimizeMolecule(iso) # تحسين الشكل الطاقي للمجسم
                
                st.subheader(f"Isomer #{i+1}")
                
                col1, col2 = st.columns([1, 2])
                
                with col1:
                    st.write("#### 📋 Stereochemical Labels:")
                    # تحديد الكيمياء الفراغية (R/S و E/Z)
                    Chem.AssignStereochemistry(iso, force=True, cleanIt=True)
                    
                    # استخراج الـ R/S (Chiral Centers)
                    chiral_centers = Chem.FindMolChiralCenters(iso, includeUnassigned=True)
                    if chiral_centers:
                        for atom_idx, label in chiral_centers:
                            st.info(f"**Center (R/S):** Atom {iso.GetAtomWithIdx(atom_idx).GetSymbol()}{atom_idx} is **{label}**")
                    
                    # استخراج الـ E/Z (Double Bonds)
                    found_ez = False
                    for bond in iso.GetBonds():
                        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                            stereo = bond.GetStereo()
                            if stereo in [Chem.rdchem.BondStereo.STEREOE, Chem.rdchem.BondStereo.STEREOZ]:
                                found_ez = True
                                label = "E (Trans-like)" if stereo == Chem.rdchem.BondStereo.STEREOE else "Z (Cis-like)"
                                st.warning(f"**Bond (E/Z):** Bond between {bond.GetBeginAtomIdx()}-{bond.GetEndAtomIdx()} is **{label}**")
                    
                    if not chiral_centers and not found_ez:
                        st.write("No specific R/S or E/Z labels detected for this structure.")

                with col2:
                    # العرض الثلاثي الأبعاد
                    view = py3Dmol.view(width=600, height=400)
                    mol_block = Chem.MolToMolBlock(iso)
                    view.addModel(mol_block, 'mol')
                    view.setStyle({'stick': {'radius': 0.15}, 'sphere': {'scale': 0.25}})
                    view.setBackgroundColor('#f0f2f6')
                    view.zoomTo()
                    showmol(view, height=400, width=600)
                
                st.divider()

    except Exception as e:
        st.error(f"Error during analysis: {e}")
