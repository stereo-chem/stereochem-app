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
            
            # --- الجزئية السحرية للألينات ---
            base_mol = Chem.AddHs(base_mol) # إضافة الهيدروجين لرؤية التعامد
            
            # 2. توليد كل الأيزومرات الممكنة (R/S) و (Axial Chirality)
            opts = EnumerateStereoisomers.StereoEnumerationOptions(tryEmbedding=True)
            isomers = list(EnumerateStereoisomers.EnumerateStereoisomers(base_mol, options=opts))
            
            st.success(f"Found {len(isomers)} Stereoisomers!")
            
            # 3. عرض كل أيزومر في عمود منفصل أو تحت بعض
            for i, iso in enumerate(isomers):
                st.subheader(f"Isomer #{i+1}")
                
                # حساب الإحداثيات الـ 3D للأيزومر ده تحديداً
                AllChem.EmbedMolecule(iso, AllChem.ETKDG())
                
                col1, col2 = st.columns([1, 2])
                
                with col1:
                    # استخراج بيانات الكيرالية (R/S)
                    Chem.AssignStereochemistry(iso, force=True, cleanIt=True)
                    st.write("**Chiral Centers/Axes:**")
                    for atom in iso.GetAtoms():
                        if atom.HasProp('_CIPCode'):
                            st.info(f"Atom {atom.GetSymbol()}{atom.GetIdx()}: Label {atom.GetProp('_CIPCode')}")

                with col2:
                    # عرض الـ 3D بنفس ستايل الصورة اللي بعتيها (Balls and Sticks)
                    view = py3Dmol.view(width=600, height=400)
                    mol_block = Chem.MolToMolBlock(iso)
                    view.addModel(mol_block, 'mol')
                    # الستايل ده بيخلي الذرات كور (Balls) والروابط عصيان (Sticks)
                    view.setStyle({'stick': {'radius': 0.15}, 'sphere': {'scale': 0.25}})
                    view.zoomTo()
                    showmol(view, height=400, width=600)
                st.divider()

    except Exception as e:
        st.error(f"Error analyzing: {e}")
