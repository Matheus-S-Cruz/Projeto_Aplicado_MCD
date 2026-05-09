import pandas as pd
import requests
import time

ARQUIVO_ENTRADA = "dados_brutos/IDENTIFICACAO.xlsx"
ARQUIVO_SAIDA = "compostos_enriquecidos.xlsx"

# Quantidade máxima de linhas para leitura
LIMITE_LINHAS = 500

# Delay entre requisições
DELAY = 0.05

# =========================
# NORMALIZAÇÃO
# =========================
def normalizar_nome(nome):
    if pd.isna(nome):
        return None

    nome = str(nome)
    nome = nome.split("[")[0]
    return nome.strip()

# =========================
# PUBCHEM
# =========================
def buscar_pubchem(nome):
    try:
        url = (
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
            f"{nome}/property/MolecularFormula,MolecularWeight,IUPACName,CID/JSON"
        )

        r = requests.get(url, timeout=10)

        if r.status_code != 200:
            return {}

        props = r.json()["PropertyTable"]["Properties"][0]

        return {
            "pubchem_cid": props.get("CID"),
            "formula": props.get("MolecularFormula"),
            "massa_pubchem": props.get("MolecularWeight"),
            "iupac_name": props.get("IUPACName"),
        }

    except Exception as e:
        print(f"Erro PubChem ({nome}): {e}")
        return {}

# =========================
# ChEBI
# =========================
def buscar_chebi(nome):
    try:
        url = (
            f"https://www.ebi.ac.uk/chebi/ws/rest/search"
            f"?searchTerm={nome}&stars=3"
        )

        r = requests.get(url, timeout=10)

        if r.status_code != 200:
            return {}

        # Algumas versões da API retornam XML.
        # Este bloco tenta interpretar JSON apenas se disponível.
        try:
            data = r.json()
        except:
            return {}

        if "listElement" not in data:
            return {}

        primeiro = data["listElement"][0]

        return {
            "chebi_id": primeiro.get("chebiId"),
            "chebi_nome": primeiro.get("chebiAsciiName"),
        }

    except Exception as e:
        print(f"Erro ChEBI ({nome}): {e}")
        return {}

# =========================
# METADATA
# =========================
def get_metadata(nome):
    nome_limpo = normalizar_nome(nome)

    if not nome_limpo:
        return {}

    return {
        **buscar_pubchem(nome_limpo),
        **buscar_chebi(nome_limpo)
    }

# =========================
# MAIN
# =========================
def main():
    print("📖 Lendo planilha...")

    # Lê apenas a coluna necessária
    df = pd.read_excel(
        ARQUIVO_ENTRADA,
        usecols=["Compound"],
        nrows=LIMITE_LINHAS
    )

    # Remove compostos duplicados
    df = df.drop_duplicates(subset=["Compound"])

    print(f"✅ Compostos únicos encontrados: {len(df)}")

    resultados = []

    # Cache para evitar consultas repetidas
    cache = {}

    for i, row in df.iterrows():
        nome = row.get("Compound")

        if pd.isna(nome):
            continue

        nome_limpo = normalizar_nome(nome)

        print(f"[{i}] 🔎 {nome_limpo}")

        # Usa cache se já consultado
        if nome_limpo in cache:
            meta = cache[nome_limpo]
            print("   ↳ usando cache")

        else:
            meta = get_metadata(nome_limpo)
            cache[nome_limpo] = meta

        resultados.append({
            "sinal_bruto": nome,
            "nome_normalizado": nome_limpo,
            **meta
        })

        # Pequena pausa para evitar bloqueio da API
        time.sleep(DELAY)

    print("💾 Salvando resultados...")

    df_final = pd.DataFrame(resultados)
    df_final.to_excel(ARQUIVO_SAIDA, index=False)

    print("✅ Enriquecimento concluído!")
    print(f"📄 Arquivo salvo em: {ARQUIVO_SAIDA}")

if __name__ == "__main__":
    main()