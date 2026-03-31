import pandas as pd
import requests
import time

ARQUIVO_ENTRADA = "dados_brutos/IDENTIFICACAO.xlsx"
ARQUIVO_SAIDA = "compostos_enriquecidos.xlsx"

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
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{nome}/property/MolecularFormula,MolecularWeight,IUPACName,CID/JSON"
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

    except:
        return {}

# =========================
# ChEBI
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
            "chebi_nome": primeiro.get("chebiAsciiName"),
        }

    except:
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
    df = pd.read_excel(ARQUIVO_ENTRADA)

    resultados = []

    for i, row in df.iterrows():
        nome = row.get("Compound")

        if pd.isna(nome):
            continue

        print(f"[{i}] 🔎 {nome}")

        meta = get_metadata(nome)

        resultados.append({
            "sinal_bruto": nome,
            "nome_normalizado": normalizar_nome(nome),
            **meta
        })

        time.sleep(0.2)

    df_final = pd.DataFrame(resultados)
    df_final.to_excel(ARQUIVO_SAIDA, index=False)

    print("✅ Enriquecimento concluído!")

if __name__ == "__main__":
    main()