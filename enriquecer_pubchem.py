import pandas as pd
import requests
import time

ARQUIVO_ENTRADA = "dados_brutos/IDENTIFICACAO.xlsx"
ARQUIVO_SAIDA = "compostos_enriquecidos.xlsx"

# =========================
# PUBCHEM
# =========================
def buscar_pubchem(nome_composto):
    try:
        nome_formatado = str(nome_composto).replace(" ", "%20")

        url = (
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
            f"{nome_formatado}/property/"
            f"MolecularFormula,MolecularWeight,IUPACName,CID/JSON"
        )

        response = requests.get(url, timeout=10)

        if response.status_code != 200:
            return None

        data = response.json()
        props = data["PropertyTable"]["Properties"][0]

        return {
            "cid": props.get("CID"),
            "formula": props.get("MolecularFormula"),
            "massa": props.get("MolecularWeight"),
            "iupac": props.get("IUPACName"),
        }

    except:
        return None


# =========================
# CHEBI (taxonomia simples)
# =========================
def buscar_chebi(nome_composto):
    try:
        url = f"https://www.ebi.ac.uk/chebi/searchId.do?chebiId={nome_composto}"
        response = requests.get(url, timeout=10)

        if response.status_code != 200:
            return None

        # Simples (placeholder — API real é mais complexa SOAP/REST)
        return {
            "chebi_id": None,
            "tipo": None  # Aqui entraria classe química
        }

    except:
        return None


# =========================
# LIMPEZA DE DADOS
# =========================
def normalizar_nome(nome):
    if pd.isna(nome):
        return None
    return str(nome).strip().lower()


# =========================
# MAIN
# =========================
def main():
    print("📂 Lendo arquivo...")
    df = pd.read_excel(ARQUIVO_ENTRADA)

    resultados = []

    for _, row in df.iterrows():
        nome = normalizar_nome(row.get("Compound"))

        if not nome:
            continue

        print(f"🔎 Buscando: {nome}")

        pubchem = buscar_pubchem(nome)
        chebi = buscar_chebi(nome)

        resultados.append({
            "composto": nome,

            # PUBCHEM
            "cid_pubchem": pubchem["cid"] if pubchem else None,
            "formula_pubchem": pubchem["formula"] if pubchem else None,
            "massa_pubchem": pubchem["massa"] if pubchem else None,
            "iupac_name": pubchem["iupac"] if pubchem else None,

            # CHEBI (futuro)
            "chebi_id": chebi["chebi_id"] if chebi else None,
            "classe_quimica": chebi["tipo"] if chebi else None,

            # PLANILHA
            "massa_planilha": row.get("Neutral mass (Da)"),
            "tempo_retencao": row.get("Retention time (min)")
        })

        time.sleep(0.2)

    df_resultado = pd.DataFrame(resultados)

    print("💾 Salvando arquivo...")
    df_resultado.to_excel(ARQUIVO_SAIDA, index=False)

    print("✅ Finalizado!")


if __name__ == "__main__":
    main()