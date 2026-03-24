import pandas as pd
import requests
import time
import numpy as np

# =========================
# CONFIG
# =========================
ARQ_ID = 'dados_brutos/IDENTIFICACAO.xlsx'
ARQ_AB = 'dados_brutos/ABUND.xlsx'

# =========================
# CACHE GLOBAL
# =========================
cache_api = {}

# =========================
# API COM CACHE
# =========================
def get_metadata(descricao):
    if descricao in cache_api:
        print("   ⚡ Cache hit")
        return cache_api[descricao]

    try:
        nome = str(descricao).replace(" ", "%20")

        url = (
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
            f"{nome}/property/MolecularFormula,MolecularWeight,IUPACName,CID/JSON"
        )

        r = requests.get(url, timeout=5)

        if r.status_code != 200:
            cache_api[descricao] = {}
            return {}

        data = r.json()
        props = data["PropertyTable"]["Properties"][0]

        resultado = {
            "pubchem_cid": props.get("CID"),
            "formula": props.get("MolecularFormula"),
            "massa_api": props.get("MolecularWeight"),
            "iupac": props.get("IUPACName")
        }

        cache_api[descricao] = resultado
        return resultado

    except:
        cache_api[descricao] = {}
        return {}

# =========================
# 1. MERGE OTIMIZADO
# =========================
print("\n📂 Lendo arquivos...")
df_id = pd.read_excel(ARQ_ID)
df_ab = pd.read_excel(ARQ_AB)

df_id.columns = df_id.columns.str.strip()
df_ab.columns = df_ab.columns.str.strip()

# remover duplicados
df_id = df_id.drop_duplicates(subset=['Compound'])
df_ab = df_ab.drop_duplicates(subset=['Compound'])

# garantir Compound ID
df_ab = df_ab.merge(
    df_id[['Compound', 'Compound ID']],
    on='Compound',
    how='left'
)

df_merged = pd.merge(
    df_id,
    df_ab,
    on='Compound',
    how='inner',
    suffixes=('', '_ab')
)

print(f"✅ Merge OK: {df_merged.shape[0]} linhas")

# =========================
# 2. ENRIQUECIMENTO OTIMIZADO
# =========================
print("\n🌐 Enriquecendo (modo otimizado)...")

# pegar apenas descrições únicas
descricoes_unicas = df_merged['Description'].dropna().unique()

print(f"🔬 Compostos únicos: {len(descricoes_unicas)}")

metadados = []

for i, desc in enumerate(descricoes_unicas[:50]):  # limite inicial
    print(f"[{i}] 🔎 {desc}")

    meta = get_metadata(desc)

    metadados.append({
        "Description": desc,
        "pubchem_cid": meta.get("pubchem_cid"),
        "formula": meta.get("formula"),
        "massa_api": meta.get("massa_api"),
        "iupac": meta.get("iupac")
    })

    time.sleep(0.2)

df_meta = pd.DataFrame(metadados)

# merge com base inteira
df_enriquecido = pd.merge(df_merged, df_meta, on="Description", how="left")

print("✅ Enriquecimento finalizado")

# =========================
# 3. RANKING OTIMIZADO
# =========================
print("\n📊 Calculando ranking...")

replicatas = [col for col in df_enriquecido.columns if "." in col]

# garantir que tudo é numérico
df_enriquecido['Score'] = pd.to_numeric(df_enriquecido['Score'], errors='coerce')

df_enriquecido['media'] = df_enriquecido[replicatas].mean(axis=1)
df_enriquecido['std'] = df_enriquecido[replicatas].std(axis=1)

df_enriquecido['media'] = pd.to_numeric(df_enriquecido['media'], errors='coerce')
df_enriquecido['std'] = pd.to_numeric(df_enriquecido['std'], errors='coerce')

# evitar divisão por zero
df_enriquecido['cv'] = df_enriquecido['std'] / df_enriquecido['media'].replace(0, 1)

# fórmula
df_enriquecido['confianca'] = (
    (1 / (1 + df_enriquecido['cv'])) *
    df_enriquecido['Score'] *
    np.log(df_enriquecido['media'] + 1)
)

df_enriquecido = df_enriquecido.dropna(subset=['confianca'])

ranking = df_enriquecido.groupby('Compound ID')['confianca'].mean()
ranking = ranking.sort_values(ascending=False)

print("\n🏆 TOP 10:")
print(ranking.head(10))

# =========================
# 4. EXPORT
# =========================
df_enriquecido.to_csv("pipeline_final.csv", index=False)
ranking.to_csv("ranking.csv")

print("\n🎯 PIPELINE OTIMIZADO FINALIZADO!")