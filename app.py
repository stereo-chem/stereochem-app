import streamlit as st
import pubchempy as pcp
import requests
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
from rdkit.Chem.Draw import rdMolDraw2D
from stmol import showmol
import py3Dmol
import numpy as np

# --- 1. إعدادات الصفحة والستايل ---
st.set_page_config(page_title="Professional Isomer System", layout="wide")

st.markdown("""
<div style="background-color: #fdf2f2; padding: 15px; border-radius: 10px; border-left: 5px solid #800000; margin-bottom: 20px;">
    <strong style="color: #800000; font-size: 1.2em;">Stereoisomerism Reference Guide:</strong><br>
    <ul style="list-style-type: none; padding-left: 0; margin-top: 10px; color: black;">
        <li>1. <b>Cis / Trans & E / Z:</b> Geometric isomers for odd-numbered double bonds (Trienes).</li>
        <li>2. <b>Ra / Sa:</b> Axial stereochemistry for even-numbered double bonds (Allenes).</li>
        <li>3. <b>R / S:</b> Absolute configuration for chiral centers.</li>
    </ul>
</div>
""", unsafe_allow_html=True)

st.markdown("<h2 style='color: #800000; font-family: serif; border-bottom: 2px solid #dcdde1;'>Chemical Isomer Analysis System 2.0 (Allenes & Trienes)</h2>", unsafe_allow_html=True)

# --- 2. الدوال المساعدة ---
def get_smiles_smart(name):
    try:
        opsin_url = f"https://opsin.ch.cam.ac.uk/opsin/{name}.json"
        res = requests.get(opsin_url)
        if res.status_code == 200: return res.json()['smiles']
    except: pass
    try:
        pcp_res = pcp.get_compounds(name, 'name')
        if pcp_res: return pcp_res[0].isomeric_smiles
    except: pass
    return None

def render_pro_2d(mol):
    mc = Chem.Mol(mol)
    Chem.AssignStereochemistry(mc, force=True, cleanIt=True)
    AllChem.Compute2DCoords(mc)
    Chem.WedgeMolBonds(mc, mc.GetConformer())
    drawer = rdMolDraw2D.MolDraw2DCairo(500, 500)
    opts = drawer.drawOptions()
    opts.bondLineWidth = 3.5
    opts.addStereoAnnotation = True
    opts.useMolBlockWedging = False
    opts.fixedBondLength = 35
    opts.explicitMethyl = True
    drawer.DrawMolecule(mc)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()

# --- 3. المعالجة الرئيسية ---
# مثال افتراضي يدعم الترايين: Hexa-2,3,4-triene أو 2,3,4-Hexatriene
compound_name = st.text_input("Enter Structure Name:", "2,3,4-Hexatriene")

if st.button("Analyze & Visualize Isomers"):
    smiles = get_smiles_smart(compound_name)
    if smiles:
        mol = Chem.MolFromSmiles(smiles)
        
        # تحسين الـ SMARTS ليشمل Allene (C=C=C) و Triene (C=C=C=C) وما فوق
        cumulene_p = Chem.MolFromSmarts("C=C=[C,C=C]") 
        
        if mol.HasSubstructMatch(cumulene_p):
            for match in mol.GetSubstructMatches(cumulene_p):
                # تحديد الذرات الطرفية في نظام الروابط المتتالية كـ مراكز فراغية
                mol.GetAtomWithIdx(match[0]).SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
        
        # استخراج جميع الأيزومرات الممكنة (تلقائياً سيعالج E/Z للترايين و Ra/Sa للألين)
        opts = StereoEnumerationOptions(tryEmbedding=True, tryStructs=True, onlyUnassigned=False)
        isomers = list(EnumerateStereoisomers(mol, options=opts))
        
        # في حال كان أيزومر واحد فقط، نقوم بتوليد الصورة المرآتية يدوياً للأمان
        if len(isomers) == 1 and mol.HasSubstructMatch(cumulene_p):
            iso2 = Chem.Mol(isomers[0])
            for a in iso2.GetAtoms():
                if a.GetChiralTag() == Chem.ChiralType.CHI_TETRAHEDRAL_CW:
                    a.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
                elif a.GetChiralTag() == Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
                    a.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
            isomers.append(iso2)

        st.subheader(f"Found {len(isomers)} Isomer(s)")
        
        st.write("---")
        cols = st.columns(len(isomers))
        
        for i, iso in enumerate(isomers):
            with cols[i]:
                # تحديد التسمية بناءً على نوع النظام (زوجي الروابط = Ra/Sa ، فردي الروابط = E/Z)
                # هنا نستخدم تسمية عامة Isomer 1, 2 مع ترك الـ RDKit يكتب الرمز فوق الرسمة
                st.markdown(f"### Isomer {i+1}", unsafe_allow_html=True)
                
                # عرض الـ 2D (سيظهر E/Z أو R/S تلقائياً بفضل addStereoAnnotation)
                st.image(render_pro_2d(iso), use_container_width=True)
                
                # عرض الـ 3D التفاعلي
                m3d = Chem.AddHs(iso)
                AllChem.EmbedMolecule(m3d, AllChem.ETKDG())
                
                view = py3Dmol.view(width=400, height=300)
                view.addModel(Chem.MolToMolBlock(m3d), 'mol')
                
                # تلوين ذرات الكربون في نظام الـ Cumulene (Allene or Triene)
                cumulene_matches = iso.GetSubstructMatches(cumulene_p)
                cumulene_atoms = set()
                for match in cumulene_matches:
                    cumulene_atoms.update(match)
                
                for idx in cumulene_atoms:
                    view.setStyle({'serial': idx+1}, {'stick': {'color':'red'}, 'sphere': {'color':'red','scale':0.3}})
                
                view.setStyle({'stick': {}, 'sphere': {'scale':0.25}})
                view.zoomTo()
                showmol(view)
    else:
        st.error("Compound not found.")
