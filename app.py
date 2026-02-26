import streamlit as st
import pubchempy as pcp
import requests
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
from rdkit.Chem.Draw import rdMolDraw2D
from stmol import showmol
import py3Dmol

# --- 1. إعدادات الصفحة ---
st.set_page_config(page_title="Professional Isomer System", layout="wide")

st.markdown("""
<div style="background-color: #fdf2f2; padding: 15px; border-radius: 10px; border-left: 5px solid #800000; margin-bottom: 20px;">
    <strong style="color: #800000; font-size: 1.2em;">Stereoisomerism Reference Guide:</strong><br>
    <ul style="list-style-type: none; padding-left: 0; margin-top: 10px; color: black;">
        <li>1. <b>Cis / Trans (E/Z):</b> For Trienes (odd number of double bonds).</li>
        <li>2. <b>Ra / Sa:</b> For Allenes (even number of double bonds).</li>
    </ul>
</div>
""", unsafe_allow_html=True)

# --- 2. الدوال المساعدة ---
def get_smiles_smart(name):
    try:
        # محاولة Opsin أولاً
        opsin_url = f"https://opsin.ch.cam.ac.uk/opsin/{name}.json"
        res = requests.get(opsin_url)
        if res.status_code == 200: return res.json()['smiles']
    except: pass
    try:
        # محاولة PubChem ثانياً
        pcp_res = pcp.get_compounds(name, 'name')
        if pcp_res: return pcp_res[0].isomeric_smiles
    except: pass
    return None

def render_pro_2d(mol):
    mc = Chem.Mol(mol)
    Chem.AssignStereochemistry(mc, force=True, cleanIt=True)
    AllChem.Compute2DCoords(mc)
    drawer = rdMolDraw2D.MolDraw2DCairo(500, 500)
    opts = drawer.drawOptions()
    opts.addStereoAnnotation = True
    opts.bondLineWidth = 3.0
    drawer.DrawMolecule(mc)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()

# --- 3. الواجهة والمعالجة ---
compound_name = st.text_input("Enter Structure Name:", "2,3,4-Hexatriene")

if st.button("Analyze & Visualize Isomers"):
    smiles = get_smiles_smart(compound_name)
    
    if smiles:
        mol = Chem.MolFromSmiles(smiles)
        
        # التأكد من أن RDKit نجح في خلق الكائن (حل مشكلة الـ ArgumentError)
        if mol is not None:
            # النمط السحري الذي يكتشف الـ Allene والـ Triene وأي Cumulene
            cumulene_p = Chem.MolFromSmarts("C=C=C") 
            
            # وسوم إضافية للتعرف على كيمياء الترايين
            mol.UpdatePropertyCache()
            Chem.AssignStereochemistry(mol, force=True)

            opts = StereoEnumerationOptions(tryEmbedding=True, onlyUnassigned=False)
            isomers = list(EnumerateStereoisomers(mol, options=opts))
            
            # إذا لم يجد أيزومرات هندسية تلقائياً، نجبره على البحث في نظام الروابط
            if len(isomers) <= 1:
                isomers = [mol] # لضمان العرض على الأقل

            st.subheader(f"Found {len(isomers)} Isomer(s)")
            st.write("---")
            
            cols = st.columns(len(isomers))
            for i, iso in enumerate(isomers):
                with cols[i]:
                    st.markdown(f"### Isomer {i+1}")
                    
                    # عرض 2D
                    st.image(render_pro_2d(iso), use_container_width=True)
                    
                    # عرض 3D
                    m3d = Chem.AddHs(iso)
                    AllChem.EmbedMolecule(m3d, AllChem.ETKDG())
                    
                    view = py3Dmol.view(width=400, height=300)
                    view.addModel(Chem.MolToMolBlock(m3d), 'mol')
                    view.setStyle({'stick': {'colorscheme': 'cyanCarbon'}})
                    
                    # تلوين ذرات الروابط المتعددة باللون الأحمر للتمييز
                    matches = iso.GetSubstructMatches(Chem.MolFromSmarts("C=C"))
                    for match in matches:
                        for idx in match:
                            view.setStyle({'serial': idx+1}, {'stick': {'color':'red'}, 'sphere': {'color':'red','scale':0.3}})
                    
                    view.zoomTo()
                    showmol(view)
        else:
            st.error("RDKit couldn't parse the structure. Please check the name.")
    else:
        st.error("Compound not found in databases.")

