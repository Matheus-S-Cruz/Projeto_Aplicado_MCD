import pandas as pd
import requests
import time
import numpy as np
import json
import re
import xml.etree.ElementTree as ET

from urllib.parse import quote
from concurrent.futures import ThreadPoolExecutor, as_completed
from sklearn.preprocessing import MinMaxScaler

# =========================
# CONFIG
# =========================
ARQ_ID = 'dados_brutos/IDENTIFICACAO.xlsx'
ARQ_AB = 'dados_brutos/ABUND.xlsx'
ARQ_FINAL = 'dados_brutos/Compostos_final.xlsx'

ARQ_CACHE = 'cache_pubchem.json'

MAX_WORKERS = 3
DELAY = 0.05
TIMEOUT = 10
SALVAR_CACHE_CADA = 50
LIMITE_COMPOSTOS = None  # None = sem limite

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
# PUBCHEM - PROPRIEDADES
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
            f"MolecularFormula,MolecularWeight,IUPACName/JSON"
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
# PUBCHEM - DESCRIÇÃO/USOS
# =========================
def buscar_pubchem_descricao(cid):
    """Busca descrição do composto no PubChem (usos e aplicações)."""

    if not cid:
        return {}

    try:
        url = (
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
            f"{int(cid)}/description/JSON"
        )

        r = session.get(url, timeout=TIMEOUT)

        if r.status_code != 200:
            return {}

        data = r.json()

        informacoes = (
            data
            .get("InformationList", {})
            .get("Information", [])
        )

        for info in informacoes:
            desc = info.get("Description", "")
            if desc and len(desc) > 20:
                return {"uso_descricao": desc[:500]}

        return {}

    except Exception:
        return {}

# =========================
# CHEBI (via PubChem proxy)
# =========================
def buscar_chebi_via_pubchem(cid):
    """Busca ChEBI ID e ontologia usando PubChem como proxy."""

    if not cid:
        return {}

    result = {}

    try:
        # 1. ChEBI ID via sinônimos
        url_syn = (
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/"
            f"compound/cid/{int(cid)}/synonyms/JSON"
        )

        r = session.get(url_syn, timeout=TIMEOUT)

        if r.status_code == 200:
            syns = (
                r.json()
                .get("InformationList", {})
                .get("Information", [{}])[0]
                .get("Synonym", [])
            )

            for s in syns:
                if 'CHEBI:' in s.upper():
                    result['chebi_id'] = s
                    break

        time.sleep(DELAY)

        # 2. Ontologia via classificação ChEBI no PubChem
        url_class = (
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/"
            f"compound/cid/{int(cid)}/classification/JSON"
        )

        r2 = session.get(url_class, timeout=15)

        if r2.status_code == 200:
            data = r2.json()

            hierarchies = (
                data
                .get("Hierarchies", {})
                .get("Hierarchy", [])
            )

            for h in hierarchies:
                src = h.get("SourceName", "")

                if src != "ChEBI":
                    continue

                nodes = h.get("Node", [])

                # extrai nome de cada nó
                def get_node_name(n):
                    info = n.get("Information", {})
                    name_obj = info.get("Name", {})
                    if isinstance(name_obj, dict):
                        swm = name_obj.get(
                            "StringWithMarkup", {}
                        )
                        return (
                            swm.get("String", "")
                            if isinstance(swm, dict)
                            else ""
                        )
                    return str(name_obj)

                # nome ChEBI (nó com Match=true)
                match_node_id = None
                for n in nodes:
                    info = n.get("Information", {})
                    if info.get("Match"):
                        result['chebi_nome'] = get_node_name(n)
                        match_node_id = n.get("NodeID")
                        break

                # pais diretos do nó match
                # (os mais específicos primeiro)
                if match_node_id:
                    # monta mapa de nós
                    node_map = {
                        n.get("NodeID"): n
                        for n in nodes
                    }

                    parent_names = []
                    current = node_map.get(match_node_id)

                    # sobe a hierarquia
                    for _ in range(5):
                        if not current:
                            break
                        parent_ids = current.get(
                            "ParentID", []
                        )
                        if not parent_ids:
                            break
                        pid = (
                            parent_ids[0]
                            if isinstance(parent_ids, list)
                            else parent_ids
                        )
                        parent_node = node_map.get(pid)
                        if not parent_node:
                            break
                        pname = get_node_name(parent_node)
                        if pname and pname not in [
                            'chemical entity',
                            'molecular entity'
                        ]:
                            parent_names.append(pname)
                        current = parent_node

                    if parent_names:
                        result['ontologia'] = (
                            ' | '.join(parent_names)
                        )

                break  # só precisa do ChEBI

        return result

    except Exception:
        return {}

# =========================
# CLASSIFICAÇÃO
# =========================
AMINOACIDOS = [
    'ala', 'arg', 'asn', 'asp', 'cys', 'gln', 'glu',
    'gly', 'his', 'ile', 'leu', 'lys', 'met', 'phe',
    'pro', 'ser', 'thr', 'trp', 'tyr', 'val'
]

PEPTIDE_NOMES = {
    2: 'dipeptídeo', 3: 'tripeptídeo',
    4: 'tetrapeptídeo', 5: 'pentapeptídeo',
    6: 'hexapeptídeo'
}

REGRAS_ONTOLOGIA = [
    (['fatty acid', 'lipid', 'acyl'],
     'Ácido graxo / Lipídio', 'Lipídios', 'Natural'),
    (['amino acid', 'amino-acid'],
     'Aminoácido', 'Aminoácidos', 'Natural'),
    (['carbohydrate', 'sugar', 'glycoside', 'monosaccharide'],
     'Carboidrato / Glicosídeo', 'Carboidratos', 'Natural'),
    (['terpene', 'terpenoid', 'isoprenoid', 'monoterpene', 'sesquiterpene', 'diterpene'],
     'Terpenoide', 'Terpenos', 'Natural'),
    (['alkaloid'],
     'Alcaloide', 'Alcaloides', 'Natural'),
    (['flavonoid', 'phenol', 'polyphenol', 'phenylpropanoid'],
     'Flavonoide / Polifenol', 'Fenólicos', 'Natural'),
    (['steroid', 'sterol'],
     'Esteroide', 'Lipídios', 'Natural'),
    (['nucleotide', 'nucleoside', 'purine', 'pyrimidine'],
     'Nucleotídeo / Nucleosídeo', 'Nucleotídeos', 'Natural'),
    (['vitamin'],
     'Vitamina', 'Cofatores', 'Natural'),
    (['organic acid', 'carboxylic acid'],
     'Ácido orgânico', 'Metabolismo primário', 'Natural'),
    (['drug', 'pharmaceutical', 'xenobiotic'],
     'Fármaco / Xenobiótico', 'Xenobiótico', 'Sintético'),
]


def classificar_composto(nome, pubchem_data, chebi_data):
    """Classifica: natural/sintético, categoria química, metabolismo."""

    nome_lower = (nome or "").lower()
    ontologia = chebi_data.get('ontologia', '').lower()

    # detecta peptídeos pelo nome
    partes = nome_lower.replace("-", " ").split()
    n_amino = sum(1 for p in partes if p in AMINOACIDOS)

    if n_amino >= 2:
        label = PEPTIDE_NOMES.get(
            n_amino, f'peptídeo ({n_amino} aa)'
        )
        return {
            'categoria_quimica': f'Peptídeo ({label})',
            'metabolismo': 'Aminoácidos',
            'tipo_composto': 'Natural'
        }

    # busca nas regras de ontologia
    for palavras, cat, met, tipo in REGRAS_ONTOLOGIA:
        if any(w in ontologia for w in palavras):
            return {
                'categoria_quimica': cat,
                'metabolismo': met,
                'tipo_composto': tipo
            }

    # fallback: nome contém "acid"
    if any(w in nome_lower for w in ['acid', 'ácido']):
        return {
            'categoria_quimica': 'Ácido orgânico',
            'metabolismo': 'Metabolismo primário',
            'tipo_composto': 'Natural'
        }

    # fallback: usa primeiro parent da ontologia
    if ontologia:
        first_parent = chebi_data.get(
            'ontologia', ''
        ).split(' | ')[0]
        return {
            'categoria_quimica': first_parent,
            'metabolismo': 'Metabolismo secundário',
            'tipo_composto': 'Natural'
        }

    return {
        'categoria_quimica': 'Não classificado',
        'metabolismo': 'Indeterminado',
        'tipo_composto': 'Indeterminado'
    }

# =========================
# MODO DE IONIZAÇÃO
# =========================
def inferir_ionizacao(adducts):
    """Infere modo de ionização a partir dos adutos."""

    if pd.isna(adducts):
        return "Indeterminado"

    adducts_str = str(adducts).lower()

    positivos = ['+h', '+na', '+k', '+nh4']
    negativos = ['-h', '+cl', '+fa', '+hac', '+cho2']

    if any(a in adducts_str for a in positivos):
        return "Positivo"

    if any(a in adducts_str for a in negativos):
        return "Negativo"

    return "Indeterminado"

# =========================
# METADATA COMPLETA
# =========================
def get_metadata(descricao):

    if descricao in cache_api:
        return cache_api[descricao]

    nome = normalizar_nome(descricao)

    if not nome:
        cache_api[descricao] = {}
        return {}

    # PubChem propriedades
    resultado_pubchem = buscar_pubchem(nome)

    # PubChem descrição/usos
    resultado_descricao = buscar_pubchem_descricao(
        resultado_pubchem.get('pubchem_cid')
    )

    # ChEBI (via PubChem synonyms + classification)
    resultado_chebi = buscar_chebi_via_pubchem(
        resultado_pubchem.get('pubchem_cid')
    )

    # classificação
    classificacao = classificar_composto(
        nome, resultado_pubchem, resultado_chebi
    )

    resultado = {
        **resultado_pubchem,
        **resultado_descricao,
        **resultado_chebi,
        **classificacao
    }

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

if LIMITE_COMPOSTOS:
    descricoes = descricoes[:LIMITE_COMPOSTOS]

print(f"🧪 Descriptions únicas: {len(descricoes)}")

# =========================
# ENRIQUECIMENTO
# =========================
print("\n🌐 Enriquecendo dados (PubChem + ChEBI)...")

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
    'iupac',
    'uso_descricao',
    'chebi_id',
    'chebi_nome',
    'ontologia',
    'chebi_definicao',
    'categoria_quimica',
    'metabolismo',
    'tipo_composto'
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
# MODO DE IONIZAÇÃO
# =========================
print("\n🔬 Inferindo modo de ionização...")

df_enriquecido['modo_ionizacao'] = (
    df_enriquecido['Adducts'].apply(inferir_ionizacao)
)

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
        'massa_api',
        'Fragmentation Score',
        'Isotope Similarity',
        'Mass Error (ppm)',
        'Neutral mass (Da)',
        'Retention time (min)',
        'Chromatographic peak width (min)',
        'Identifications'
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

for col in ['Score', 'Fragmentation Score', 'Isotope Similarity']:
    if col in df_enriquecido.columns:
        df_enriquecido[col] = pd.to_numeric(
            df_enriquecido[col],
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

# metadata válida (PubChem + ChEBI)
df_enriquecido['metadata_ok'] = np.where(
    df_enriquecido['pubchem_cid'].notna()
    | df_enriquecido['chebi_id'].notna(),
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

# Fragmentation Score normalizado
if 'Fragmentation Score' in df_enriquecido.columns:
    df_enriquecido['frag_norm'] = normalizar_coluna(
        df_enriquecido,
        'Fragmentation Score'
    )
else:
    df_enriquecido['frag_norm'] = 0

# Isotope Similarity normalizado
if 'Isotope Similarity' in df_enriquecido.columns:
    df_enriquecido['isotope_norm'] = normalizar_coluna(
        df_enriquecido,
        'Isotope Similarity'
    )
else:
    df_enriquecido['isotope_norm'] = 0

# score final ponderado (inclui Fragmentation e Isotope)
df_enriquecido['confianca'] = (
    0.25 * df_enriquecido['media_norm'] +
    0.20 * df_enriquecido['score_norm'] +
    0.15 * df_enriquecido['frag_norm'] +
    0.15 * df_enriquecido['isotope_norm'] +
    0.15 * df_enriquecido['estabilidade'] +
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
            'Fragmentation Score',
            'Isotope Similarity',
            'pubchem_cid',
            'formula',
            'chebi_id',
            'ontologia',
            'categoria_quimica',
            'metabolismo',
            'tipo_composto',
            'uso_descricao'
        ]
    ]
    .sort_values(
        by='confianca',
        ascending=False
    )
)

# =========================
# COMPOSTOS_FINAL (FORMATO MODELO)
# =========================
print("\n📋 Gerando planilha no formato Compostos_final...")

# agrupa por Description (um registro por composto)
df_compostos = (
    df_enriquecido
    .sort_values('confianca', ascending=False)
    .drop_duplicates(subset=['Description'])
)

# monta ID sequencial
df_compostos = df_compostos.reset_index(drop=True)
df_compostos['ID_seq'] = [
    f"COMP_{str(i+1).zfill(4)}"
    for i in range(len(df_compostos))
]

compostos_final = pd.DataFrame({
    'ID': df_compostos['ID_seq'],
    'Metabólito/Composto': df_compostos['Description'],
    'Fórmula': df_compostos['formula'],
    'Massa Molecular': df_compostos['massa_api'],
    'IUPAC': df_compostos['iupac'],
    'Modo de Ionização': df_compostos['modo_ionizacao'],
    'Categoria química': df_compostos['categoria_quimica'],
    'Tipo (Natural/Sintético)': df_compostos['tipo_composto'],
    'Metabolismo': df_compostos['metabolismo'],
    'Ontologia (ChEBI)': df_compostos['ontologia'],
    'PubChem CID': df_compostos['pubchem_cid'],
    'ChEBI ID': df_compostos['chebi_id'],
    'Usos/Aplicações': df_compostos['uso_descricao'],
    'Score Confiança': df_compostos['confianca'].round(4),
    'Score Identificação': df_compostos['Score'],
    'Fragmentation Score': df_compostos['Fragmentation Score'],
    'Isotope Similarity': df_compostos['Isotope Similarity'],
    'Média Abundância': df_compostos['media'].round(2),
    'CV (%)': (df_compostos['cv'] * 100).round(2)
})

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
# ESTATÍSTICAS FINAIS
# =========================
total = len(compostos_final)
com_pubchem = compostos_final['PubChem CID'].notna().sum()
com_chebi = compostos_final['ChEBI ID'].notna().sum()
com_classif = (
    compostos_final['Categoria química'] != 'Não classificado'
).sum()
naturais = (
    compostos_final['Tipo (Natural/Sintético)'] == 'Natural'
).sum()
sinteticos = (
    compostos_final['Tipo (Natural/Sintético)'] == 'Sintético'
).sum()

print(f"\n📊 RESUMO:")
print(f"   Total compostos únicos: {total}")
print(f"   Com PubChem CID: {com_pubchem} ({100*com_pubchem/max(total,1):.1f}%)")
print(f"   Com ChEBI ID: {com_chebi} ({100*com_chebi/max(total,1):.1f}%)")
print(f"   Classificados: {com_classif} ({100*com_classif/max(total,1):.1f}%)")
print(f"   Naturais: {naturais} | Sintéticos: {sinteticos}")

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

compostos_final.to_excel(
    "compostos_final_resultado.xlsx",
    index=False
)

print("\n🎯 PIPELINE FINALIZADO!")
print("📄 Arquivos gerados:")
print("   - pipeline_final.csv (dados completos)")
print("   - ranking.csv (ranking resumido)")
print("   - ranking_detalhado.csv (ranking com detalhes)")
print("   - compostos_final_resultado.xlsx (formato modelo)")
