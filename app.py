import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers

# 1. إعدادات الصفحة
st.set_page_config(page_title="Chemical Isomer Analysis", layout="wide")

# 2. تصميم الواجهة
st.markdown("""
<h2 style='color: #800000; font-family: serif; border-bottom: 2px solid #dcdde1;'>Chemical Isomer Analysis System</h2>
<div style="background-color: #ffffff; padding: 15px; border: 1px solid #e1e1e1; border-left: 4px solid #800000; margin-bottom: 20px; font-family: sans-serif;">
    <strong style="color: #800000;">Stereoisomerism Reference Guide:</strong><br>
    1. <b style="color: #b22222;">Cis / Trans (Relative):</b> Identical groups on same/opposite sides.<br>
    2. <b style="color: #b22222;">E / Z (Absolute - CIP System):</b> <b>Z (Zusammen)</b> together, <b>E (Entgegen)</b> opposite.<br>
    3. <b style="color: #b22222;">R / S (Optical):</b> Absolute configuration of chiral centers.
</div>
""", unsafe_allow_html=True)

# 3. مدخلات المستخدم
compound_name = st.text_input("Enter Structure Name (e.g., 1,2-dichloroethene or Thalidomide):", "")

if st.button("Analyze Isomers"):
    if not compound_name:
        st.warning("Please enter a compound name first.")
    else:
        try:
            results = pcp.get_compounds(compound_name, 'name')
            
            if not results:
                st.error(f"❌ No compound found for: {compound_name}")
            else:
                base_smiles = results[0].smiles
                mol = Chem.MolFromSmiles(base_smiles)
                
                mol_no_stereo = Chem.Mol(mol)
                for bond in mol_no_stereo.GetBonds():
                    bond.SetStereo(Chem.BondStereo.STEREONONE)
                
                isomers = list(EnumerateStereoisomers(mol_no_stereo))
                
                # --- التعديل هنا: التحقق من الـ Achiral ---
                if len(isomers) <= 1:
                    st.info(f"✨ The compound {compound_name} is Achiral. It does not have geometric (E/Z) or optical (R/S) isomers.")
                # ------------------------------------------

                labels = []
                for i, iso in enumerate(isomers):
                    Chem.AssignStereochemistry(iso, force=True, cleanIt=True)
                    stereo_info = []
                    
                    for bond in iso.GetBonds():
                        stereo = bond.GetStereo()
                        if stereo == Chem.BondStereo.STEREOE: 
                            stereo_info.append("E")
                        elif stereo == Chem.BondStereo.STEREOZ: 
                            stereo_info.append("Z")
                    
                    chiral_centers = Chem.FindMolChiralCenters(iso)
                    for center in chiral_centers:
                        stereo_info.append(f"({center[1]})")
                    
                    # لو مفيش معلومات فراغية بنكتب Achiral تحت الصورة
                    label = f"Isomer {i+1}: " + (", ".join(stereo_info) if stereo_info else "Achiral Structure")
                    labels.append(label)

                st.success(f"Analyzed Structure: {compound_name} | Total Forms: **{len(isomers)}**")
                
                img = Draw.MolsToGridImage(
                    isomers,
                    molsPerRow=3,
                    subImgSize=(350, 350),
                    legends=labels
                )
                
                st.image(img, use_container_width=True)

        except Exception as e:
            st.error(f"An unexpected error occurred: {e}")

st.markdown("---")
st.caption("Powered by RDKit, PubChemPy, and Streamlit.")
