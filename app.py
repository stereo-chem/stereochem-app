import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw, AllChem, EnumerateStereoisomers
from rdkit.Chem.EnumerateStereoisomers import StereoEnumerationOptions
from stmol import showmol
import py3Dmol

# إعدادات واجهة الموبايل
st.set_page_config(page_title="Chemical Isomer Pro", layout="wide")

st.markdown("""
<style>
    .main { background-color: #ffffff; }
    h2 { color: #800000; text-align: center; font-family: 'Serif'; }
    .stButton>button { width: 100%; border-radius: 20px; background-color: #800000; color: white; }
</style>
<h2>Chemical Isomer Analysis (Allenes Support)</h2>
""", unsafe_allow_html=True)

def render_3d(mol, title):
    try:
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        mblock = Chem.MolToMolBlock(mol)
        view = py3Dmol.view(width=400, height=300)
        view.addModel(mblock, 'mol')
        view.setStyle({'stick': {}, 'sphere': {'scale': 0.3}})
        view.zoomTo()
        st.write(f"**{title}**")
        showmol(view, height=300, width=400)
    except:
        st.warning(f"Could not generate 3D for {title}")

# مدخلات المستخدم
compound_name = st.text_input("Enter Structure Name (e.g., 2,3-pentadiene):", "")

if st.button("Analyze Isomers"):
    if compound_name:
        try:
            results = pcp.get_compounds(compound_name, 'name')
            if not results:
                st.error("❌ Compound not found.")
            else:
                base_smiles = results[0].smiles
                mol = Chem.MolFromSmiles(base_smiles)
                
                # إعدادات خاصة للألينات لضمان اكتشاف الكيرالية المحورية
                opts = StereoEnumerationOptions(tryEmbedding=True, onlyUnassigned=False)
                isomers = list(EnumerateStereoisomers(mol, options=opts))
                
                st.info(f"Detected {len(isomers)} potential stereoisomers.")
                
                # عرض النتائج في أعمدة تناسب الشاشة
                for i, iso in enumerate(isomers):
                    with st.container():
                        render_3d(iso, f"Isomer #{i+1}")
                        st.divider()
        except Exception as e:
            st.error(f"Error: {e}")
    else:
        st.warning("Please enter a name.")
