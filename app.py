import streamlit as st
import pubchempy as pcp
import requests
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
from rdkit.Chem.Draw import rdMolDraw2D
from stmol import showmol
import py3Dmol

# --- 1. إعدادات الصفحة والستايل ---
st.set_page_config(page_title="Professional Isomer System", layout="wide")

st.markdown("""
<div style="background-color: #fdf2f2; padding: 15px; border-radius: 10px; border-left: 5px solid #800000; margin-bottom: 20px;">
    <strong style="color: #800000; font-size: 1.2em;">Stereoisomerism Reference Guide:</strong>
    <ul style="color: black; margin-top: 10px;">
        <li>1. <b>Cis / Trans:</b> Identical groups on same/opposite sides.</li>
        <li>2. <b>E / Z:</b> High-priority groups together (Z) or opposite (E).</li>
        <li>3. <b>R / S:</b> Absolute configuration of chiral centers.</li>
        <li>4. <b>Ra / Sa:</b> Stereochemistry of Allenes (C=C=C).</li>
    </ul>
</div>
""", unsafe_allow_html=True)

st.markdown("<h2 style='color: #800000; font-family: serif; border-bottom: 2px solid #dcdde1;'>Chemical Isomer Analysis System 2.0</h2>", unsafe_allow_html=True)

# --- 2. الدوال المساعدة ---
def get_smiles_smart(name):
    # محاولة البحث عبر OPSIN أولاً
    try:
        url = f"https://opsin.ch.cam.ac.uk/opsin/{name}.json"
        res = requests.get(url)
        if res.status_code == 200: return res.json()['smiles']
    except: pass
    # محاولة البحث عبر PubChem كبديل
    try:
        pcp_res = pcp.get_compounds(name, 'name')
        if pcp_res: return pcp_res[0].isomeric_smiles
    except: pass
    return None

def render_pro_2d(mol):
    # إنشاء نسخة للرسم لضمان إبراز الكيمياء الفراغية
    mc = Chem.Mol(mol)
    Chem.AssignStereochemistry(mc, force=True, cleanIt=True)
    
    # حساب الإحداثيات وتوليد روابط Wedge/Dash لإظهار Solid و Dashed
    AllChem.Compute2DCoords(mc)
    Chem.WedgeMolBonds(mc, mc.GetConformer())
    
    # تحديد منطقة الألين لتلوينها (Highlight)
    allene_p = Chem.MolFromSmarts("C=C=C")
    matches = mc.GetSubstructMatches(allene_p)
    
    atoms_to_h = []
    bonds_to_h = []
    a_colors = {}
    b_colors = {}
    h_color = (1.0, 0.85, 0.85) # لون وردي هادئ للتظليل

    for match in matches:
        atoms_to_h.extend(match)
        for idx in match: a_colors[idx] = h_color
        for i in range(len(match)-1):
            bond = mc.GetBondBetweenAtoms(match[i], match[i+1])
            if bond:
                bidx = bond.GetIdx()
                bonds_to_h.append(bidx)
                b_colors[bidx] = h_color

    # إعدادات الرسام ليكون بوضوح احترافي (السمك 5.0 يبرز الـ Dashed بوضوح)
    drawer = rdMolDraw2D.MolDraw2DCairo(600, 600)
    opts = drawer.drawOptions()
    opts.bondLineWidth = 5.0           
    opts.addStereoAnnotation = True    
    opts.fixedBondLength = 45          
    opts.isomericInchi = True          
    
    drawer.DrawMolecule(mc, highlightAtoms=atoms_to_h, highlightAtomColors=a_colors,
                       highlightBonds=bonds_to_h, highlightBondColors=b_colors)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()

# --- 3. المعالجة الرئيسية ---
# إدخال الاسم العام (البرنامج سيقوم بتوليد R و S تلقائياً)
raw_name = st.text_input("Enter Structure Name (e.g., 2,3-pentadiene):", "1,3-Dimethyl-3-phenylallene")

if st.button("Analyze & Visualize Isomers"):
    smiles = get_smiles_smart(raw_name)
    if smiles:
        # تحويل SMILES إلى جزيء خام وإزالة أي تجسيم مسبق للبحث عن كل الاحتمالات
        base_mol = Chem.MolFromSmiles(smiles)
        Chem.RemoveStereochemistry(base_mol)
        
        # تفعيل البحث عن الآيزومرات الفراغية للألين
        allene_p = Chem.MolFromSmarts("C=C=C")
        if base_mol.HasSubstructMatch(allene_p):
            for match in base_mol.GetSubstructMatches(allene_p):
                base_mol.GetAtomWithIdx(match[0]).SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)

        # توليد كل الآيزومرات الممكنة تلقائياً
        opts = StereoEnumerationOptions(tryEmbedding=True, onlyUnassigned=False)
        isomers = list(EnumerateStereoisomers(base_mol, options=opts))
        
        # ضمان ظهور الزوج (Ra/Sa) في حالة الألين
        if len(isomers) == 1 and base_mol.HasSubstructMatch(allene_p):
            iso2 = Chem.Mol(isomers[0])
            for a in iso2.GetAtoms():
                if a.GetChiralTag() == Chem.ChiralType.CHI_TETRAHEDRAL_CW:
                    a.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
            isomers.append(iso2)

        st.subheader(f"Analysis for: {raw_name}")
        st.info(f"Detected **{len(isomers)}** possible stereoisomers.")
        
        cols = st.columns(len(isomers))
        for i, iso in enumerate(isomers):
            with cols[i]:
                # تسمية تلقائية للآيزومرات
                label = "Ra" if i == 0 else "Sa"
                st.markdown(f"### Isomer {i+1}: <span style='color: #800000;'>{label}</span>", unsafe_allow_html=True)
                
                # عرض الرسم 2D المطور (Wedge/Dash)
                st.image(render_pro_2d(iso), use_container_width=True)
                
                # عرض 3D التفاعلي
                m3d = Chem.AddHs(iso)
                AllChem.EmbedMolecule(m3d, AllChem.ETKDG())
                view = py3Dmol.view(width=400, height=300)
                view.addModel(Chem.MolToMolBlock(m3d), 'mol')
                view.setStyle({'stick': {}, 'sphere': {'scale':0.25}})
                view.zoomTo()
                showmol(view)
    else:
        st.error("Compound not found. Try names like '2,3-pentadiene' or '1,3-diphenylallene'")
