"""
QuimioAnalytics — Motor de ETL (versão refatorada do pipeline.py).

Diferente do pipeline.py (script único com caminhos fixos, mantido como
referência histórica), aqui o processamento é exposto em funções reutilizáveis
que:
  - aceitam arquivos enviados (caminho OU objeto de arquivo, ex.: upload do Streamlit);
  - rodam o merge + enriquecimento (PubChem/ChEBI) + estatística + ranking;
  - persistem o resultado no banco (models.py).

Modos de enriquecimento:
  - usar_cache_apenas=False  -> consulta as APIs (correto, porém lento);
  - usar_cache_apenas=True   -> usa só o cache_pubchem.json (rápido, ideal p/ demo).

Ponto de entrada principal: processar(...).
"""

from __future__ import annotations

import datetime as dt
import json
import os
import re
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd
import requests

import models

# =========================
# CONFIG
# =========================
ARQ_CACHE = "cache_pubchem.json"
HIST_DIR = "historico"
MAX_WORKERS = 3
DELAY = 0.05
TIMEOUT = 10

session_http = requests.Session()


# =========================
# CACHE
# =========================
def carregar_cache(caminho: str = ARQ_CACHE) -> dict:
    try:
        with open(caminho, "r", encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return {}


def salvar_cache(cache: dict, caminho: str = ARQ_CACHE) -> None:
    with open(caminho, "w", encoding="utf-8") as f:
        json.dump(cache, f, ensure_ascii=False, indent=4)


# =========================
# NORMALIZAÇÃO DE NOME
# =========================
def normalizar_nome(nome):
    if pd.isna(nome):
        return None
    nome = str(nome).split("[")[0].strip()
    nome = re.sub(r"[^a-zA-Z0-9\s\-\(\),]", "", nome)
    nome = re.sub(r"\s+", " ", nome)
    if len(nome) > 80:
        return None
    if nome.lower() in {"unknown", "unidentified", "untitled", "na", "null"}:
        return None
    return nome


# =========================
# PUBCHEM / CHEBI
# =========================
def buscar_pubchem(nome):
    if not nome:
        return {}
    estrategias = [nome, nome.split(",")[0], " ".join(nome.split()[:4])]
    for tentativa_nome in estrategias:
        url = (
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
            f"{tentativa_nome}/property/MolecularFormula,MolecularWeight,IUPACName/JSON"
        )
        for _ in range(2):
            try:
                r = session_http.get(url, timeout=TIMEOUT)
                if r.status_code != 200:
                    continue
                props = r.json().get("PropertyTable", {}).get("Properties", [])
                if not props:
                    continue
                p = props[0]
                return {
                    "pubchem_cid": p.get("CID"),
                    "formula": p.get("MolecularFormula"),
                    "massa_api": p.get("MolecularWeight"),
                    "iupac": p.get("IUPACName"),
                }
            except Exception:
                time.sleep(0.5)
    return {}


def buscar_pubchem_descricao(cid):
    if not cid:
        return {}
    try:
        url = (
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
            f"{int(cid)}/description/JSON"
        )
        r = session_http.get(url, timeout=TIMEOUT)
        if r.status_code != 200:
            return {}
        infos = r.json().get("InformationList", {}).get("Information", [])
        for info in infos:
            desc = info.get("Description", "")
            if desc and len(desc) > 20:
                return {"uso_descricao": desc[:500]}
    except Exception:
        pass
    return {}


def buscar_chebi_via_pubchem(cid):
    if not cid:
        return {}
    result = {}
    try:
        url_syn = (
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{int(cid)}/synonyms/JSON"
        )
        r = session_http.get(url_syn, timeout=TIMEOUT)
        if r.status_code == 200:
            syns = (
                r.json().get("InformationList", {}).get("Information", [{}])[0].get("Synonym", [])
            )
            for s in syns:
                if "CHEBI:" in s.upper():
                    result["chebi_id"] = s
                    break
        time.sleep(DELAY)

        url_class = (
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{int(cid)}/classification/JSON"
        )
        r2 = session_http.get(url_class, timeout=15)
        if r2.status_code == 200:
            hierarchies = r2.json().get("Hierarchies", {}).get("Hierarchy", [])
            for h in hierarchies:
                if h.get("SourceName") != "ChEBI":
                    continue
                nodes = h.get("Node", [])

                def get_node_name(n):
                    info = n.get("Information", {})
                    name_obj = info.get("Name", {})
                    if isinstance(name_obj, dict):
                        swm = name_obj.get("StringWithMarkup", {})
                        return swm.get("String", "") if isinstance(swm, dict) else ""
                    return str(name_obj)

                match_node_id = None
                for n in nodes:
                    if n.get("Information", {}).get("Match"):
                        result["chebi_nome"] = get_node_name(n)
                        match_node_id = n.get("NodeID")
                        break

                if match_node_id:
                    node_map = {n.get("NodeID"): n for n in nodes}
                    parent_names = []
                    current = node_map.get(match_node_id)
                    for _ in range(5):
                        if not current:
                            break
                        parent_ids = current.get("ParentID", [])
                        if not parent_ids:
                            break
                        pid = parent_ids[0] if isinstance(parent_ids, list) else parent_ids
                        parent_node = node_map.get(pid)
                        if not parent_node:
                            break
                        pname = get_node_name(parent_node)
                        if pname and pname not in ("chemical entity", "molecular entity"):
                            parent_names.append(pname)
                        current = parent_node
                    if parent_names:
                        result["ontologia"] = " | ".join(parent_names)
                break
    except Exception:
        pass
    return result


# =========================
# CLASSIFICAÇÃO
# =========================
AMINOACIDOS = [
    "ala", "arg", "asn", "asp", "cys", "gln", "glu", "gly", "his", "ile",
    "leu", "lys", "met", "phe", "pro", "ser", "thr", "trp", "tyr", "val",
]
PEPTIDE_NOMES = {2: "dipeptídeo", 3: "tripeptídeo", 4: "tetrapeptídeo", 5: "pentapeptídeo", 6: "hexapeptídeo"}
REGRAS_ONTOLOGIA = [
    (["fatty acid", "lipid", "acyl"], "Ácido graxo / Lipídio", "Lipídios", "Natural"),
    (["amino acid", "amino-acid"], "Aminoácido", "Aminoácidos", "Natural"),
    (["carbohydrate", "sugar", "glycoside", "monosaccharide"], "Carboidrato / Glicosídeo", "Carboidratos", "Natural"),
    (["terpene", "terpenoid", "isoprenoid", "monoterpene", "sesquiterpene", "diterpene"], "Terpenoide", "Terpenos", "Natural"),
    (["alkaloid"], "Alcaloide", "Alcaloides", "Natural"),
    (["flavonoid", "phenol", "polyphenol", "phenylpropanoid"], "Flavonoide / Polifenol", "Fenólicos", "Natural"),
    (["steroid", "sterol"], "Esteroide", "Lipídios", "Natural"),
    (["nucleotide", "nucleoside", "purine", "pyrimidine"], "Nucleotídeo / Nucleosídeo", "Nucleotídeos", "Natural"),
    (["vitamin"], "Vitamina", "Cofatores", "Natural"),
    (["organic acid", "carboxylic acid"], "Ácido orgânico", "Metabolismo primário", "Natural"),
    (["drug", "pharmaceutical", "xenobiotic"], "Fármaco / Xenobiótico", "Xenobiótico", "Sintético"),
]


def classificar_composto(nome, chebi_data):
    nome_lower = (nome or "").lower()
    ontologia = chebi_data.get("ontologia", "").lower()

    partes = nome_lower.replace("-", " ").split()
    n_amino = sum(1 for p in partes if p in AMINOACIDOS)
    if n_amino >= 2:
        label = PEPTIDE_NOMES.get(n_amino, f"peptídeo ({n_amino} aa)")
        return {"categoria_quimica": f"Peptídeo ({label})", "metabolismo": "Aminoácidos", "tipo_composto": "Natural"}

    for palavras, cat, met, tipo in REGRAS_ONTOLOGIA:
        if any(w in ontologia for w in palavras):
            return {"categoria_quimica": cat, "metabolismo": met, "tipo_composto": tipo}

    if any(w in nome_lower for w in ["acid", "ácido"]):
        return {"categoria_quimica": "Ácido orgânico", "metabolismo": "Metabolismo primário", "tipo_composto": "Natural"}

    if ontologia:
        first_parent = chebi_data.get("ontologia", "").split(" | ")[0]
        return {"categoria_quimica": first_parent, "metabolismo": "Metabolismo secundário", "tipo_composto": "Natural"}

    return {"categoria_quimica": "Não classificado", "metabolismo": "Indeterminado", "tipo_composto": "Indeterminado"}


def inferir_ionizacao(adducts):
    if pd.isna(adducts):
        return "Indeterminado"
    s = str(adducts).lower()
    if any(a in s for a in ["+h", "+na", "+k", "+nh4"]):
        return "Positivo"
    if any(a in s for a in ["-h", "+cl", "+fa", "+hac", "+cho2"]):
        return "Negativo"
    return "Indeterminado"


# =========================
# ENRIQUECIMENTO (com cache)
# =========================
def get_metadata(descricao, cache, usar_cache_apenas):
    if descricao in cache:
        return cache[descricao]
    if usar_cache_apenas:
        return {}  # modo demo: não bate na API se não estiver no cache

    nome = normalizar_nome(descricao)
    if not nome:
        cache[descricao] = {}
        return {}

    pub = buscar_pubchem(nome)
    desc = buscar_pubchem_descricao(pub.get("pubchem_cid"))
    chebi = buscar_chebi_via_pubchem(pub.get("pubchem_cid"))
    classif = classificar_composto(nome, chebi)

    resultado = {**pub, **desc, **chebi, **classif}
    cache[descricao] = resultado
    return resultado


# =========================
# UTILITÁRIOS DE DADOS
# =========================
def detectar_colunas_replicata(df: pd.DataFrame) -> list[str]:
    """Colunas no padrão Grupo.Replicata (ex.: '1.1', '2.2')."""
    cols = []
    for c in df.columns:
        if re.fullmatch(r"\d+\.\d+", str(c).strip()):
            cols.append(c)
    return cols


def _norm_minmax(serie: pd.Series) -> np.ndarray:
    s = pd.to_numeric(serie, errors="coerce").fillna(0)
    lo, hi = s.min(), s.max()
    if hi - lo == 0:
        return np.zeros(len(s))
    return ((s - lo) / (hi - lo)).to_numpy()


# =========================
# PROCESSAMENTO PRINCIPAL
# =========================
def processar(
    arquivo_identificacao,
    arquivo_abundancia,
    Session,
    usar_cache_apenas: bool = True,
    limite: int | None = None,
    progresso=None,
) -> dict:
    """Lê os 2 arquivos, enriquece, calcula estatística/ranking e grava no banco.

    Retorna um dicionário-resumo. `progresso` é um callback opcional(str, float)
    para atualizar uma barra de progresso (ex.: st.progress).
    """
    def _log(msg, frac=None):
        if progresso:
            progresso(msg, frac)

    cache = carregar_cache()

    # 1) LEITURA
    _log("Lendo planilhas...", 0.05)
    df_id = pd.read_excel(arquivo_identificacao)
    df_ab = pd.read_excel(arquivo_abundancia)
    df_id.columns = df_id.columns.str.strip()
    df_ab.columns = df_ab.columns.str.strip()

    if "Compound ID" not in df_id.columns and "Compound" in df_id.columns:
        df_id["Compound ID"] = df_id["Compound"]

    # 2) LIMPEZA + MERGE
    df_id = df_id.drop_duplicates(subset=["Compound"]).dropna(subset=["Description"])
    df_ab = df_ab.drop_duplicates(subset=["Compound"])
    df = pd.merge(df_id, df_ab, on="Compound", how="inner", suffixes=("", "_ab"))
    _log(f"Merge concluído: {len(df)} compostos.", 0.15)

    if limite:
        df = df.head(limite)

    # 3) ENRIQUECIMENTO
    descricoes = df["Description"].dropna().astype(str).unique()
    _log(f"Enriquecendo {len(descricoes)} compostos...", 0.2)
    metadados = []
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as ex:
        futures = {
            ex.submit(get_metadata, d, cache, usar_cache_apenas): d for d in descricoes
        }
        for i, fut in enumerate(as_completed(futures)):
            d = futures[fut]
            try:
                metadados.append({"Description": d, **fut.result()})
            except Exception:
                metadados.append({"Description": d})
            if i % 25 == 0:
                _log(f"Enriquecidos {i}/{len(descricoes)}...", 0.2 + 0.4 * i / max(len(descricoes), 1))

    if not usar_cache_apenas:
        salvar_cache(cache)

    df_meta = pd.DataFrame(metadados)
    for col in ["pubchem_cid", "formula", "massa_api", "iupac", "uso_descricao",
                "chebi_id", "ontologia", "categoria_quimica", "metabolismo", "tipo_composto"]:
        if col not in df_meta.columns:
            df_meta[col] = np.nan
    df = pd.merge(df, df_meta, on="Description", how="left")

    # 4) MODO DE IONIZAÇÃO
    df["modo_ionizacao"] = df["Adducts"].apply(inferir_ionizacao) if "Adducts" in df.columns else "Indeterminado"

    # 5) ESTATÍSTICA (média/std/cv sobre as replicatas)
    _log("Calculando estatísticas...", 0.65)
    reps = detectar_colunas_replicata(df)
    for c in reps:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    for c in ["Score", "Fragmentation Score", "Isotope Similarity"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    df["media"] = df[reps].mean(axis=1) if reps else np.nan
    df["std"] = df[reps].std(axis=1) if reps else np.nan
    df["cv"] = np.where(df["media"] > 0, df["std"] / df["media"], np.nan)
    df["cv"] = df["cv"].replace([np.inf, -np.inf], np.nan)

    # 6) SCORE DE CONFIANÇA (fórmula ponderada, normalização min-max manual)
    df["media_norm"] = _norm_minmax(np.log1p(df["media"].fillna(0)))
    df["score_norm"] = _norm_minmax(df.get("Score", pd.Series(0, index=df.index)))
    df["frag_norm"] = _norm_minmax(df.get("Fragmentation Score", pd.Series(0, index=df.index)))
    df["isotope_norm"] = _norm_minmax(df.get("Isotope Similarity", pd.Series(0, index=df.index)))
    cv_fill = df["cv"].fillna(df["cv"].median())
    df["estabilidade"] = 1 - _norm_minmax(cv_fill)
    df["metadata_ok"] = np.where(df["pubchem_cid"].notna() | df["chebi_id"].notna(), 1, 0)
    df["confianca"] = (
        0.25 * df["media_norm"] + 0.20 * df["score_norm"] + 0.15 * df["frag_norm"]
        + 0.15 * df["isotope_norm"] + 0.15 * df["estabilidade"] + 0.10 * df["metadata_ok"]
    ).clip(lower=0)

    # 7) RANKING
    df = df.sort_values("confianca", ascending=False).reset_index(drop=True)
    df["posicao"] = df.index + 1

    # 8) PERSISTÊNCIA
    _log("Gravando no banco...", 0.8)
    n = _persistir(df, reps, Session)

    resumo = {
        "compostos": int(len(df)),
        "com_pubchem": int(df["pubchem_cid"].notna().sum()),
        "com_chebi": int(df["chebi_id"].notna().sum()),
        "medicoes": n["medicoes"],
        "amostras": n["amostras"],
        "modo": "cache" if usar_cache_apenas else "API",
    }

    # 9) HISTÓRICO — salva um retrato (snapshot) desta análise
    _log("Registrando no histórico...", 0.95)
    df_ist = documento_ist(Session)
    ana = registrar_analise(Session, df_ist, resumo, arquivo_identificacao)
    resumo["analise_id"] = ana["id"]
    resumo["analise_nome"] = ana["nome"]

    _log("Concluído.", 1.0)
    return resumo


def _persistir(df: pd.DataFrame, reps: list[str], Session) -> dict:
    """Grava o DataFrame processado nas 8 tabelas do banco (limpa dados antigos)."""
    from models import (
        Abundancia, Amostra, Composto, CompostoTag, MetricaComposto,
        Ontologia, ProcedenciaCampo, RankingComposto, Tag,
    )

    with Session() as s:
        # limpa dados de um processamento anterior (mantém as tags do seed)
        for M in (Abundancia, MetricaComposto, RankingComposto, Ontologia,
                  ProcedenciaCampo, CompostoTag, Composto, Amostra):
            s.query(M).delete()
        s.commit()

        # tags do catálogo (nome -> id)
        tags_por_nome = {t.nome: t for t in s.query(Tag).all()}

        # amostras (uma por grupo das colunas Grupo.Replicata)
        grupos = sorted({int(str(c).split(".")[0]) for c in reps})
        amostras = {}
        for g in grupos:
            am = Amostra(id_amostra=f"AM{g}", nome=f"Amostra {g}", grupo=str(g), tipo="amostra")
            s.add(am)
            amostras[g] = am.id_amostra

        n_medicoes = 0
        for _, row in df.iterrows():
            cid = str(row["Compound"])

            # média por grupo -> amostra mais abundante
            medias_grupo = {}
            for c in reps:
                g = int(str(c).split(".")[0])
                v = row.get(c)
                if pd.notna(v):
                    medias_grupo.setdefault(g, []).append(v)
            medias_grupo = {g: np.mean(vs) for g, vs in medias_grupo.items() if vs}
            grupo_top = max(medias_grupo, key=medias_grupo.get) if medias_grupo else None

            comp = Composto(
                compound_id=cid,
                compound_id_ext=(str(row["Compound ID"]) if pd.notna(row.get("Compound ID")) else None),
                nome=str(row["Compound"]),
                description=(str(row["Description"]) if pd.notna(row.get("Description")) else None),
                pubchem_cid=(int(row["pubchem_cid"]) if pd.notna(row.get("pubchem_cid")) else None),
                chebi_id=(str(row["chebi_id"]) if pd.notna(row.get("chebi_id")) else None),
                formula=(str(row["formula"]) if pd.notna(row.get("formula")) else (str(row["Formula"]) if pd.notna(row.get("Formula")) else None)),
                massa_molecular=(float(row["massa_api"]) if pd.notna(row.get("massa_api")) else (float(row["Neutral mass (Da)"]) if pd.notna(row.get("Neutral mass (Da)")) else None)),
                iupac=(str(row["iupac"]) if pd.notna(row.get("iupac")) else None),
                descricao=(str(row["uso_descricao"]) if pd.notna(row.get("uso_descricao")) else None),
                classe_geral=(str(row["categoria_quimica"]) if pd.notna(row.get("categoria_quimica")) else None),
                subclasse=(row["ontologia"].split(" | ")[0] if isinstance(row.get("ontologia"), str) else None),
                tipo_composto=(str(row["tipo_composto"]) if pd.notna(row.get("tipo_composto")) else None),
                metabolismo=(str(row["metabolismo"]) if pd.notna(row.get("metabolismo")) else None),
            )
            s.add(comp)

            s.add(MetricaComposto(
                compound_id=cid,
                modo_aquisicao=row.get("modo_ionizacao"),
                adducts=(str(row["Adducts"]) if pd.notna(row.get("Adducts")) else None),
                link=(str(row["Link"]) if pd.notna(row.get("Link")) else None),
                score=(float(row["Score"]) if pd.notna(row.get("Score")) else None),
                fragmentation_score=(float(row["Fragmentation Score"]) if pd.notna(row.get("Fragmentation Score")) else None),
                isotope_similarity=(float(row["Isotope Similarity"]) if pd.notna(row.get("Isotope Similarity")) else None),
                mass_error_ppm=(float(row["Mass Error (ppm)"]) if pd.notna(row.get("Mass Error (ppm)")) else None),
                neutral_mass=(float(row["Neutral mass (Da)"]) if pd.notna(row.get("Neutral mass (Da)")) else None),
                mz=(float(row["m/z"]) if pd.notna(row.get("m/z")) else None),
                retention_time_min=(float(row["Retention time (min)"]) if pd.notna(row.get("Retention time (min)")) else None),
                peak_width_min=(float(row["Chromatographic peak width (min)"]) if pd.notna(row.get("Chromatographic peak width (min)")) else None),
                identifications=(int(row["Identifications"]) if pd.notna(row.get("Identifications")) else None),
                media_abundancia=(float(row["media"]) if pd.notna(row.get("media")) else None),
                desvio_padrao=(float(row["std"]) if pd.notna(row.get("std")) else None),
                cv=(float(row["cv"]) if pd.notna(row.get("cv")) else None),
            ))

            s.add(RankingComposto(
                compound_id=cid,
                confianca=(float(row["confianca"]) if pd.notna(row.get("confianca")) else None),
                posicao=int(row["posicao"]),
                original_id=cid,
            ))

            # ontologia (um nível por linha)
            if isinstance(row.get("ontologia"), str) and row["ontologia"]:
                for nivel, termo in enumerate(row["ontologia"].split(" | ")):
                    s.add(Ontologia(compound_id=cid, nivel=nivel, termo=termo.strip()))

            # procedência (origem dos dados enriquecidos)
            if pd.notna(row.get("pubchem_cid")):
                s.add(ProcedenciaCampo(compound_id=cid, campo="formula/massa/descricao", fonte="PubChem", referencia=f"CID {int(row['pubchem_cid'])}"))
            if pd.notna(row.get("chebi_id")):
                s.add(ProcedenciaCampo(compound_id=cid, campo="ontologia/classe", fonte="ChEBI", referencia=str(row["chebi_id"])))

            # tags derivadas da abundância (o formato simples não traz tags)
            media = row.get("media")
            for limiar, nome_tag in [(500, "Abund > 500"), (1000, "Abund > 1000"), (5000, "Abund > 5000"), (10000, "Abund > 10000")]:
                if pd.notna(media) and media > limiar and nome_tag in tags_por_nome:
                    s.add(CompostoTag(compound_id=cid, tag_id=tags_por_nome[nome_tag].id))
            if pd.isna(row.get("Fragmentation Score")) and "Not Fragmented" in tags_por_nome:
                s.add(CompostoTag(compound_id=cid, tag_id=tags_por_nome["Not Fragmented"].id))

            # abundância (uma linha por replicata)
            for c in reps:
                v = row.get(c)
                if pd.isna(v):
                    continue
                g = int(str(c).split(".")[0])
                r = int(str(c).split(".")[1])
                s.add(Abundancia(
                    compound_id=cid,
                    id_amostra=amostras[g],
                    abundancia=float(v),
                    tipo="normalised",
                    replicata=r,
                    is_maxima=(g == grupo_top),
                ))
                n_medicoes += 1

        s.commit()

    return {"medicoes": n_medicoes, "amostras": len(grupos)}


# =========================
# CONSULTA: TABELA FINAL "DOCUMENTO IST"
# =========================
def documento_ist(Session) -> pd.DataFrame:
    """Monta a tabela final (formato Documento IST) via JOIN sobre o banco."""
    from sqlalchemy import select
    from models import Abundancia, Amostra, Composto, MetricaComposto, RankingComposto

    with Session() as s:
        stmt = (
            select(
                Composto.compound_id.label("Composto"),
                Composto.compound_id_ext.label("Composto ID"),
                MetricaComposto.modo_aquisicao.label("Modo de aquisição"),
                MetricaComposto.score.label("Score"),
                MetricaComposto.fragmentation_score.label("Fragmentação"),
                MetricaComposto.media_abundancia.label("Abund. relativa"),
                Composto.description.label("Descrição"),
                Composto.classe_geral.label("Classe geral"),
                Composto.subclasse.label("Subclasse"),
                RankingComposto.confianca.label("Confiança"),
                RankingComposto.posicao.label("Posição"),
            )
            .join(MetricaComposto, MetricaComposto.compound_id == Composto.compound_id, isouter=True)
            .join(RankingComposto, RankingComposto.compound_id == Composto.compound_id, isouter=True)
            .order_by(RankingComposto.posicao)
        )
        df = pd.read_sql(stmt, s.bind)

        # amostra mais abundante (marcada com is_maxima)
        amax = pd.read_sql(
            select(Abundancia.compound_id, Amostra.nome)
            .join(Amostra, Amostra.id_amostra == Abundancia.id_amostra)
            .where(Abundancia.is_maxima.is_(True))
            .distinct(),
            s.bind,
        ).rename(columns={"compound_id": "Composto", "nome": "Amostra mais abundante"})

    df = df.merge(amax, on="Composto", how="left")
    ordem = ["Composto", "Composto ID", "Modo de aquisição", "Score", "Fragmentação",
             "Abund. relativa", "Amostra mais abundante", "Descrição", "Classe geral",
             "Subclasse", "Confiança", "Posição"]
    return df[[c for c in ordem if c in df.columns]]


# =========================
# HISTÓRICO DE ANÁLISES (Opção B — snapshots)
# =========================
def _nome_arquivo(x) -> str:
    """Extrai um nome legível de um caminho ou de um upload do Streamlit."""
    if hasattr(x, "name"):
        return os.path.basename(str(x.name))
    return os.path.basename(str(x))


def registrar_analise(Session, df_ist: pd.DataFrame, resumo: dict, arquivo_identificacao) -> dict:
    """Salva o retrato (snapshot) da análise atual e registra na tabela `analise`."""
    from models import Analise

    os.makedirs(HIST_DIR, exist_ok=True)
    agora = dt.datetime.now()
    caminho = os.path.join(HIST_DIR, f"analise_{agora.strftime('%Y%m%d_%H%M%S')}.pkl")
    df_ist.to_pickle(caminho)

    nome = f"{agora.strftime('%d/%m/%Y %H:%M')} — {_nome_arquivo(arquivo_identificacao)}"
    with Session() as s:
        a = Analise(
            nome=nome,
            modo=resumo.get("modo"),
            n_compostos=resumo.get("compostos"),
            n_medicoes=resumo.get("medicoes"),
            n_amostras=resumo.get("amostras"),
            arquivo=_nome_arquivo(arquivo_identificacao),
            snapshot=caminho,
        )
        s.add(a)
        s.commit()
        return {"id": a.id, "nome": a.nome}


def listar_analises(Session) -> list[dict]:
    """Lista as análises registradas, da mais recente para a mais antiga.

    Auto-limpeza: registros cujo arquivo de snapshot não existe mais (ex.: apagado
    manualmente da pasta `historico/`) são removidos do banco e não aparecem na lista.
    """
    from models import Analise

    with Session() as s:
        registros = s.query(Analise).order_by(Analise.criado_em.desc(), Analise.id.desc()).all()
        validos, orfaos = [], []
        for a in registros:
            if a.snapshot and os.path.exists(a.snapshot):
                validos.append({
                    "id": a.id, "nome": a.nome, "modo": a.modo,
                    "n_compostos": a.n_compostos, "n_medicoes": a.n_medicoes,
                    "n_amostras": a.n_amostras, "snapshot": a.snapshot,
                })
            else:
                orfaos.append(a)
        for a in orfaos:  # remove registros sem arquivo (snapshot apagado na mão)
            s.delete(a)
        if orfaos:
            s.commit()
        return validos


def carregar_analise(snapshot_path: str) -> pd.DataFrame:
    """Carrega o retrato (Documento IST) de uma análise do histórico."""
    return pd.read_pickle(snapshot_path)


def excluir_analise(Session, analise_id: int, snapshot_path: str | None = None) -> None:
    """Remove uma análise do histórico (registro no banco + arquivo do snapshot)."""
    from models import Analise

    with Session() as s:
        a = s.get(Analise, analise_id)
        if a is not None:
            s.delete(a)
            s.commit()
    if snapshot_path and os.path.exists(snapshot_path):
        try:
            os.remove(snapshot_path)
        except OSError:
            pass
