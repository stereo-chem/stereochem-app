import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import AllChem, EnumerateStereoisomers
from stmol import showmol
import py3Dmol
import random

# إعداد واجهة التطبيق
st.set_page_config(page_title="Universal Isomer Analyzer", layout="wide")

# دالة لتوليد لون خلفية عشوائي مريح للعين (فاتح أو داكن حسب الرغبة)
def get_random_color():
    colors = ['#f0f2f6', '#e1e5eb', '#d1d8e0', '#f8f9fa', '#2c3e50', '#34495e', '#1a1a1a']
    return random.choice(colors)

st.title("🧪 Universal Isomer & Stereo Analysis")
st.markdown("تحليل شامل لجميع أنواع الأيزومرات: **E/Z, Cis/Trans, R/S, Axial Chirality**")

# إدخال اسم المركب
name = st.text_input("Enter Molecule Name:", "2,3-pentadiene")

if st.button("Generate All Possible Isomers"):
    # اختيار لون عشوائي للخلفية عند كل ضغطة زر
    bg_color = get_random_color()
    
    try:
        # 1. جلب المركب من قاعدة البيانات
        results = pcp.get_compounds(name, 'name')
        if not results:
            st.error("Molecule not found in PubChem!")
        else:
            smiles = results[0].isomeric_smiles
            # إنشاء الجزيء الأساسي مع إضافة الهيدروجين (ضروري للـ E/Z والـ Allenes)
            base_mol = Chem.MolFromSmiles(smiles)
            base_mol = Chem.AddHs(base_mol)
            
            # 2. توليد جميع الاحتمالات الفراغية الممكنة
            # هذه الدالة تغطي (R/S) و (E/Z) و (Cis/Trans)
            opts = EnumerateStereoisomers.StereoEnumerationOptions(
                tryEmbedding=True, 
                onlyUnassigned=False # توليد الكل حتى لو تم تحديد البعض في الاسم
            )
            isomers = list(EnumerateStereoisomers.EnumerateStereoisomers(base_mol, options=opts))
            
            st.success(f"✅ Found {len(isomers)} potential stereoisomers for '{name}'")
            
            # 3. عرض الأيزومرات
            for i, iso in enumerate(isomers):
                st.markdown(f"### 📍 Isomer #{i+1}")
                
                # حساب الإحداثيات الـ 3D بدقة عالية
                AllChem.EmbedMolecule(iso, AllChem.ETKDG())
                Chem.AssignStereochemistry(iso, force=True, cleanIt=True)
                
                col1, col2 = st.columns([1, 2])
                
                with col1:
                    st.info("**Stereo Labels Detected:**")
                    found_labels = False
                    # البحث عن تسميات R/S على الذرات
                    for atom in iso.GetAtoms():
                        if atom.HasProp('_CIPCode'):
                            st.write(f"🔹 Atom {atom.GetSymbol()}{atom.GetIdx()}: **{atom.GetProp('_CIPCode')}**")
                            found_labels = True
                    
                    # البحث عن تسميات E/Z على الروابط المضاعفة
                    for bond in iso.GetBonds():
                        if bond.HasProp('_CIPCode'):
                            st.write(f"🔸 Bond {bond.GetBeginAtom().GetSymbol()}{bond.GetBeginAtomIdx()}={bond.GetEndAtom().GetSymbol()}{bond.GetEndAtomIdx()}: **{bond.GetProp('_CIPCode')}**")
                            found_labels = True
                            
                    if not found_labels:
                        st.write("No specific R/S or E/Z labels (Check 3D structure for Cis/Trans)")

                with col2:
                    # بناء عرض الـ 3D مع لون الخلفية العشوائي
                    view = py3Dmol.view(width=700, height=450)
                    mol_block = Chem.MolToMolBlock(iso)
                    view.addModel(mol_block, 'mol')
                    
                    # ستايل الكرات والعصي (Ball and Stick)
                    view.setStyle({'stick': {'radius': 0.15}, 'sphere': {'scale': 0.25}})
                    
                    # تطبيق لون الخلفية العشوائي
                    view.setBackgroundColor(bg_color)
                    
                    view.zoomTo()
                    showmol(view, height=450, width=700)
                
                st.divider()

    except Exception as e:
        st.error(f"An error occurred: {e}")
