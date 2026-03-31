import pandas as pd
import requests
import time
import numpy as np

# =========================
# CONFIG
# =========================
ARQ_ID = 'dados_brutos/IDENTIFICACAO.xlsx'
ARQ_AB = 'dados_brutos/ABUND.xlsx'
ARQ_FINAL = 'dados_brutos/compostos_final.xlsx'

# =========================
# CACHE GLOBAL
# =========================
cache_api = {}

# =========================
# NORMALIZAÇÃO
# =========================
def normalizar_nome(nome):
    if pd.isna(nome):
        return None

    nome = str(nome)
    nome = nome.split("[")[0]
    nome = nome.strip()

    return nome

# =========================
# PUBCHEM
# =========================
def buscar_pubchem(nome):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{nome}/property/MolecularFormula,MolecularWeight,IUPACName,CID/JSON"

    try:
        r = requests.get(url, timeout=10)
        if r.status_code != 200:
            return {}

        props = r.json()["PropertyTable"]["Properties"][0]

        return {
            "pubchem_cid": props.get("CID"),
            "formula": props.get("MolecularFormula"),
            "massa_api": props.get("MolecularWeight"),
            "iupac": props.get("IUPACName")
        }

    except:
        return {}

# =========================
# ChEBI (REST real simplificado)
# =========================
def buscar_chebi(nome):
    try:
        url = f"https://www.ebi.ac.uk/chebi/ws/rest/search?searchTerm={nome}&stars=3"
        r = requests.get(url, timeout=10)

        if r.status_code != 200:
            return {}

        data = r.json()

        if "listElement" not in data:
            return {}

        primeiro = data["listElement"][0]

        return {
            "chebi_id": primeiro.get("chebiId"),
            "chebi_name": primeiro.get("chebiAsciiName"),
            "score_chebi": primeiro.get("searchScore")
        }

    except:
        return {}

# =========================
# METADATA COMPLETO
# =========================
def get_metadata(descricao):
    if descricao in cache_api:
        return cache_api[descricao]

    nome = normalizar_nome(descricao)

    if not nome:
        return {}

    pubchem = buscar_pubchem(nome)
    chebi = buscar_chebi(nome)

    resultado = {
        **pubchem,
        **chebi
    }

    cache_api[descricao] = resultado
    return resultado


# =========================
# 1. MERGE (CORRIGIDO)
# =========================
print("\n📂 Lendo arquivos...")
df_id = pd.read_excel(ARQ_ID)
df_ab = pd.read_excel(ARQ_AB)

df_id.columns = df_id.columns.str.strip()
df_ab.columns = df_ab.columns.str.strip()

df_id = df_id.drop_duplicates(subset=['Compound'])
df_ab = df_ab.drop_duplicates(subset=['Compound'])

# Foi removido o merge do df_ab consigo mesmo que causava o erro.
# Agora fazemos direto o inner join unindo a identificação com a abundância.
df_merged = pd.merge(
    df_id,
    df_ab,
    on='Compound',
    how='inner'
)

print(f"✅ Merge OK: {df_merged.shape[0]} linhas")

# =========================
# 2. ENRIQUECIMENTO
# =========================
print("\n🌐 Enriquecendo dados...")

descricoes = df_merged['Description'].dropna().unique()

metadados = []

for i, desc in enumerate(descricoes):
    print(f"[{i}] 🔎 {desc}")

    meta = get_metadata(desc)

    metadados.append({
        "Description": desc,
        **meta
    })

    time.sleep(0.2)

df_meta = pd.DataFrame(metadados)

df_enriquecido = pd.merge(df_merged, df_meta, on="Description", how="left")

print("✅ Enriquecimento concluído")

# =========================
# 3. RANKING
# =========================
print("\n📊 Calculando ranking...")

replicatas = [col for col in df_enriquecido.columns if "." in col]

df_enriquecido[replicatas] = df_enriquecido[replicatas].apply(pd.to_numeric, errors='coerce')
df_enriquecido['Score'] = pd.to_numeric(df_enriquecido['Score'], errors='coerce')

df_enriquecido['media'] = df_enriquecido[replicatas].mean(axis=1)
df_enriquecido['std'] = df_enriquecido[replicatas].std(axis=1)

df_enriquecido['cv'] = df_enriquecido['std'] / df_enriquecido['media'].replace(0, 1)

df_enriquecido['confianca'] = (
    (1 / (1 + df_enriquecido['cv'])) *
    df_enriquecido['Score'] *
    np.log(df_enriquecido['media'] + 1)
)

ranking = (
    df_enriquecido
    .dropna(subset=['confianca'])
    .groupby('Compound ID')['confianca']
    .mean()
    .sort_values(ascending=False)
)

print("\n🏆 TOP 10:")
print(ranking.head(10))

# =========================
# 4. VALIDAÇÃO
# =========================
print("\n🔍 Validando com compostos_final...")

try:
    df_ref = pd.read_excel(ARQ_FINAL)

    df_ref.columns = df_ref.columns.str.strip()

    df_validacao = pd.merge(
        df_enriquecido,
        df_ref[['Compound ID']],
        on='Compound ID',
        how='inner'
    )

    taxa = len(df_validacao) / len(df_enriquecido) * 100

    print(f"✅ Taxa de correspondência: {taxa:.2f}%")

except Exception as e:
    print("⚠️ Não foi possível validar:", e)

# =========================
# 5. EXPORT
# =========================
df_enriquecido.to_csv("pipeline_final.csv", index=False)
ranking.to_csv("ranking.csv")

print("\n🎯 PIPELINE FINALIZADO COM SUCESSO!")