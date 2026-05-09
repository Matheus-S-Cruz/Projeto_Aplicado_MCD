import pandas as pd
import requests
import time
import numpy as np
import json
import re

from urllib.parse import quote
from concurrent.futures import ThreadPoolExecutor, as_completed
from sklearn.preprocessing import MinMaxScaler

# =========================
# CONFIG
# =========================
ARQ_ID = 'dados_brutos/IDENTIFICACAO.xlsx'
ARQ_AB = 'dados_brutos/ABUND.xlsx'
ARQ_FINAL = 'dados_brutos/compostos_final.xlsx'

ARQ_CACHE = 'cache_pubchem.json'

MAX_WORKERS = 3
DELAY = 0.05
TIMEOUT = 5
SALVAR_CACHE_CADA = 50

# =========================
# SESSION HTTP
# =========================
session = requests.Session()

# =========================
# CACHE
# =========================
try:
    with open(ARQ_CACHE, 'r', encoding='utf-8') as f:
        cache_api = json.load(f)

    print(f"✅ Cache carregado: {len(cache_api)} itens")

except:
    cache_api = {}
    print("⚠️ Novo cache criado")

# =========================
# NORMALIZAÇÃO
# =========================
def normalizar_nome(nome):

    if pd.isna(nome):
        return None

    nome = str(nome)

    nome = nome.split("[")[0]
    nome = nome.strip()

    # remove caracteres problemáticos
    nome = re.sub(r'[^a-zA-Z0-9\s\-\(\),]', '', nome)

    # remove espaços duplos
    nome = re.sub(r'\s+', ' ', nome)

    # remove nomes absurdos
    if len(nome) > 80:
        return None

    # blacklist
    blacklist = [
        "unknown",
        "unidentified",
        "untitled",
        "na",
        "null"
    ]

    if nome.lower() in blacklist:
        return None

    return nome

# =========================
# CACHE SAVE
# =========================
def salvar_cache():

    with open(ARQ_CACHE, 'w', encoding='utf-8') as f:
        json.dump(cache_api, f, ensure_ascii=False, indent=4)

# =========================
# PUBCHEM
# =========================
def buscar_pubchem(nome):

    if not nome:
        return {}

    estrategias = [
        nome,
        nome.split(",")[0],
        " ".join(nome.split()[:4]),
        quote(nome.replace("-", " "))
    ]

    for tentativa_nome in estrategias:

        url = (
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
            f"{tentativa_nome}/property/"
            f"MolecularFormula,MolecularWeight,IUPACName,CID/JSON"
        )

        for tentativa in range(2):

            try:

                r = session.get(
                    url,
                    timeout=TIMEOUT
                )

                if r.status_code != 200:
                    continue

                data = r.json()

                props = (
                    data
                    .get("PropertyTable", {})
                    .get("Properties", [])
                )

                if not props:
                    continue

                props = props[0]

                return {
                    "pubchem_cid": props.get("CID"),
                    "formula": props.get("MolecularFormula"),
                    "massa_api": props.get("MolecularWeight"),
                    "iupac": props.get("IUPACName")
                }

            except Exception:
                time.sleep(0.5)

    return {}

# =========================
# METADATA
# =========================
def get_metadata(descricao):

    if descricao in cache_api:
        return cache_api[descricao]

    nome = normalizar_nome(descricao)

    if not nome:
        cache_api[descricao] = {}
        return {}

    resultado = buscar_pubchem(nome)

    cache_api[descricao] = resultado

    return resultado

# =========================
# THREAD
# =========================
def enriquecer_descricao(desc):

    print(f"🔎 {desc[:80]}")

    meta = get_metadata(desc)

    time.sleep(DELAY)

    return {
        "Description": desc,
        **meta
    }

# =========================
# LEITURA
# =========================
print("\n📂 Lendo arquivos...")

df_id = pd.read_excel(ARQ_ID)
df_ab = pd.read_excel(ARQ_AB)

df_id.columns = df_id.columns.str.strip()
df_ab.columns = df_ab.columns.str.strip()

# =========================
# GARANTE COMPOUND ID
# =========================
if 'Compound ID' not in df_id.columns:

    if 'Compound' in df_id.columns:
        df_id['Compound ID'] = df_id['Compound']

# =========================
# LIMPEZA
# =========================
print("\n🧹 Limpando dados...")

df_id = df_id.drop_duplicates(subset=['Compound'])
df_ab = df_ab.drop_duplicates(subset=['Compound'])

df_id = df_id.dropna(subset=['Description'])

# =========================
# MERGE
# =========================
print("\n🔗 Realizando merge...")

df_merged = pd.merge(
    df_id,
    df_ab,
    on='Compound',
    how='inner'
)

print(f"✅ Merge OK: {df_merged.shape[0]} linhas")

# =========================
# DESCRIPTIONS
# =========================
descricoes = (
    df_merged['Description']
    .dropna()
    .astype(str)
    .unique()
)

print(f"🧪 Descriptions únicas: {len(descricoes)}")

# =========================
# ENRIQUECIMENTO
# =========================
print("\n🌐 Enriquecendo dados...")

metadados = []

contador = 0

with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:

    futures = {
        executor.submit(
            enriquecer_descricao,
            desc
        ): desc

        for desc in descricoes
    }

    for future in as_completed(futures):

        try:

            resultado = future.result()

            metadados.append(resultado)

            contador += 1

            if contador % SALVAR_CACHE_CADA == 0:

                salvar_cache()

                print(f"💾 Cache parcial salvo ({contador})")

        except Exception as e:
            print(f"⚠️ Erro paralelo: {e}")

# salva final
salvar_cache()

print("💾 Cache final salvo")


# =========================
# METADATA DF
# =========================
df_meta = pd.DataFrame(metadados)

# garante colunas mesmo vazias
colunas_metadata = [
    'Description',
    'pubchem_cid',
    'formula',
    'massa_api',
    'iupac'
]

for col in colunas_metadata:

    if col not in df_meta.columns:
        df_meta[col] = np.nan

df_enriquecido = pd.merge(
    df_merged,
    df_meta,
    on='Description',
    how='left'
)

print("✅ Enriquecimento concluído")

# =========================
# REPLICATAS
# =========================
print("\n📊 Detectando replicatas...")

replicatas = []

for col in df_enriquecido.columns:

    nome = str(col).lower()

    if (
        "sample" in nome or
        "rep" in nome or
        "abund" in nome
    ):
        replicatas.append(col)

print(f"✅ Replicatas encontradas: {len(replicatas)}")

# fallback
if len(replicatas) == 0:

    numericas = df_enriquecido.select_dtypes(
        include=[np.number]
    ).columns.tolist()

    remover = [
        'Score',
        'pubchem_cid',
        'massa_api'
    ]

    replicatas = [
        c for c in numericas
        if c not in remover
    ]

    print(f"✅ Fallback replicatas: {len(replicatas)}")

# =========================
# NUMÉRICO
# =========================
df_enriquecido[replicatas] = (
    df_enriquecido[replicatas]
    .apply(pd.to_numeric, errors='coerce')
)

df_enriquecido['Score'] = pd.to_numeric(
    df_enriquecido['Score'],
    errors='coerce'
)

# =========================
# ESTATÍSTICAS
# =========================
print("\n📈 Calculando estatísticas...")

df_enriquecido['media'] = (
    df_enriquecido[replicatas]
    .mean(axis=1)
)

df_enriquecido['std'] = (
    df_enriquecido[replicatas]
    .std(axis=1)
)

df_enriquecido['cv'] = np.where(
    df_enriquecido['media'] > 0,
    df_enriquecido['std'] / df_enriquecido['media'],
    np.nan
)

df_enriquecido['cv'] = (
    df_enriquecido['cv']
    .replace([np.inf, -np.inf], np.nan)
)

# =========================
# NORMALIZAÇÃO SEGURA
# =========================
def normalizar_coluna(df, coluna):

    valores = df[[coluna]].fillna(0)

    if valores.nunique().iloc[0] <= 1:
        return np.zeros(len(df))

    scaler = MinMaxScaler()

    return scaler.fit_transform(valores).flatten()

# log para reduzir impacto de outliers
df_enriquecido['media_log'] = np.log1p(
    df_enriquecido['media']
)

df_enriquecido['media_norm'] = normalizar_coluna(
    df_enriquecido,
    'media_log'
)

df_enriquecido['score_norm'] = normalizar_coluna(
    df_enriquecido,
    'Score'
)

# =========================
# SCORE CIENTÍFICO
# =========================
print("\n🧠 Calculando score científico...")

scaler = MinMaxScaler()

# metadata válida
df_enriquecido['metadata_ok'] = np.where(
    df_enriquecido['pubchem_cid'].notna(),
    1,
    0
)

# remove outliers extremos
limite = df_enriquecido['media'].quantile(0.99)

df_enriquecido = df_enriquecido[
    df_enriquecido['media'] <= limite
]

df_enriquecido['score_norm'] = scaler.fit_transform(
    df_enriquecido[['Score']].fillna(0)
)

cv_temp = df_enriquecido['cv'].fillna(
    df_enriquecido['cv'].median()
)

df_enriquecido['cv_norm'] = scaler.fit_transform(
    cv_temp.values.reshape(-1, 1)
)

# estabilidade
df_enriquecido['estabilidade'] = (
    1 - df_enriquecido['cv_norm']
)

# score final ponderado
df_enriquecido['confianca'] = (
    0.40 * df_enriquecido['media_norm'] +
    0.30 * df_enriquecido['score_norm'] +
    0.20 * df_enriquecido['estabilidade'] +
    0.10 * df_enriquecido['metadata_ok']
)

df_enriquecido['confianca'] = (
    df_enriquecido['confianca']
    .clip(lower=0)
)

# =========================
# RANKING
# =========================
print("\n🏆 Gerando ranking...")

ranking = (
    df_enriquecido
    .groupby('Description')['confianca']
    .mean()
    .dropna()
    .sort_values(ascending=False)
)

print("\n🏆 TOP 10:")
print(ranking.head(10))

# =========================
# RANKING DETALHADO
# =========================
ranking_detalhado = (
    df_enriquecido[
        [
            'Description',
            'confianca',
            'media',
            'cv',
            'Score',
            'pubchem_cid',
            'formula'
        ]
    ]
    .sort_values(
        by='confianca',
        ascending=False
    )
)

# =========================
# VALIDAÇÃO
# =========================
print("\n🔍 Validando...")

try:

    df_ref = pd.read_excel(ARQ_FINAL)

    df_ref.columns = (
        df_ref.columns.str.strip()
    )

    if 'Compound ID' in df_ref.columns:

        df_validacao = pd.merge(
            df_enriquecido,
            df_ref[['Compound ID']],
            on='Compound ID',
            how='inner'
        )

        taxa = (
            len(df_validacao)
            / max(len(df_enriquecido), 1)
        ) * 100

        print(
            f"✅ Taxa de correspondência: "
            f"{taxa:.2f}%"
        )

    else:
        print(
            "⚠️ 'Compound ID' não encontrado"
        )

except Exception as e:
    print(f"⚠️ Erro validação: {e}")

# =========================
# EXPORT
# =========================
print("\n💾 Exportando arquivos...")

df_enriquecido.to_csv(
    "pipeline_final.csv",
    index=False
)

ranking.to_csv(
    "ranking.csv",
    header=True
)

ranking_detalhado.to_csv(
    "ranking_detalhado.csv",
    index=False
)

print("\n🎯 PIPELINE FINALIZADO!")