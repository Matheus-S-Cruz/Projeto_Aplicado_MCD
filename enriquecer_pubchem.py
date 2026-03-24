import pandas as pd
import requests
import time

# =========================
# CONFIGURAÇÃO
# =========================
ARQUIVO_ENTRADA = "dados_brutos/IDENTIFICACAO.xlsx"
ARQUIVO_SAIDA = "compostos_enriquecidos.xlsx"

# =========================
# BUSCA POR MASSA (PUBCHEM)
# =========================
def buscar_pubchem_por_massa(massa, tentativas=2):
    url = (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/"
        f"fastformula/{massa}/property/"
        f"MolecularFormula,MolecularWeight,IUPACName,CID/JSON"
    )

    for tentativa in range(tentativas):
        try:
            print(f"   🌐 Tentativa {tentativa+1}...")

            response = requests.get(url, timeout=5)

            if response.status_code != 200:
                print("   ❌ Massa não encontrada")
                return None

            data = response.json()
            props = data["PropertyTable"]["Properties"][0]

            print("   ✅ Encontrado por massa!")

            return {
                "cid": props.get("CID"),
                "formula": props.get("MolecularFormula"),
                "massa": props.get("MolecularWeight"),
                "iupac": props.get("IUPACName"),
            }

        except requests.exceptions.Timeout:
            print("   ⏱️ Timeout! Tentando novamente...")

        except Exception as e:
            print(f"   ⚠️ Erro: {e}")
            return None

    print("   ❌ Falhou após tentativas")
    return None


# =========================
# FUNÇÃO PRINCIPAL
# =========================
def main():
    print("📂 Lendo arquivo...")
    df = pd.read_excel(ARQUIVO_ENTRADA)

    print("📊 Colunas encontradas:")
    print(df.columns)

    # TESTE RÁPIDO (remova depois)
    # df = df.head(10)

    resultados = []

    for index, row in df.iterrows():
        massa = row.get("Neutral mass (Da)")
        tempo_retencao = row.get("Retention time (min)")
        sinal_bruto = row.get("Compound")

        print(f"\n[{index}] 🔎 Sinal: {sinal_bruto}")
        print(f"   ⚖️ Massa: {massa}")

        if pd.isna(massa):
            print("   ⚠️ Massa inválida, pulando...")
            continue

        dados_api = buscar_pubchem_por_massa(massa)

        resultados.append({
            # SINAL BRUTO
            "sinal_bruto": sinal_bruto,

            # DADOS EXPERIMENTAIS
            "massa_planilha": massa,
            "tempo_retencao": tempo_retencao,

            # PUBCHEM
            "cid_pubchem": dados_api["cid"] if dados_api else None,
            "formula_pubchem": dados_api["formula"] if dados_api else None,
            "massa_pubchem": dados_api["massa"] if dados_api else None,
            "iupac_name": dados_api["iupac"] if dados_api else None,
        })

        # evitar bloqueio da API
        time.sleep(0.2)

    print("\n💾 Gerando arquivo final...")
    df_resultado = pd.DataFrame(resultados)
    df_resultado.to_excel(ARQUIVO_SAIDA, index=False)

    print(f"✅ Arquivo gerado: {ARQUIVO_SAIDA}")


# =========================
# EXECUÇÃO
# =========================
if __name__ == "__main__":
    main()