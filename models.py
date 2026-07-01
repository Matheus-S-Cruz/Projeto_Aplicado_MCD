"""
Modelagem do banco local do QuimioAnalytics (SQLAlchemy 2.0).

Esquema estrela:
    composto (dim) ── identidade + enriquecimento químico (Progenesis + APIs)
    amostra  (dim) ── amostras/replicatas
    abundancia (fato) ── composto × amostra × replicata × tipo (normalised/raw)
    metrica_composto (1:1) ── métricas analíticas do Progenesis
    ranking_composto (1:1) ── confiança + campos da ordenação biológica
    ontologia (1:N) ── hierarquia ChEBI
    tag / composto_tag (M:N) ── tags configuradas pela equipe no Progenesis
    procedencia_campo (1:N) ── rastreabilidade da origem de cada informação

A tabela final "Documento IST" (Composto, Composto ID, Modo de aquisição,
Score, Fragmentação, Abund. relativa, Amostra mais abundante, Descrição,
Classe geral, Subclasse) é uma consulta sobre este modelo, não uma tabela.
"""

from __future__ import annotations

import datetime as dt

from sqlalchemy import (
    Boolean,
    Date,
    DateTime,
    Float,
    ForeignKey,
    Integer,
    String,
    Text,
    UniqueConstraint,
    create_engine,
    func,
)
from sqlalchemy.orm import (
    DeclarativeBase,
    Mapped,
    mapped_column,
    relationship,
    sessionmaker,
)


class Base(DeclarativeBase):
    pass


# =========================================================
# DIMENSÃO: COMPOSTO  (identidade + enriquecimento químico)
# =========================================================
class Composto(Base):
    __tablename__ = "composto"

    # --- Progenesis (interno) ---
    compound_id: Mapped[str] = mapped_column(String(80), primary_key=True)  # ex.: 9.44_426.2758n ("Compound")
    compound_id_ext: Mapped[str | None] = mapped_column(String(80))         # "Compound ID" externo (CSID) -> coluna "Composto ID"
    nome: Mapped[str] = mapped_column(String(255), nullable=False)          # "Compound"
    description: Mapped[str | None] = mapped_column(Text)                   # chave de enriquecimento

    # --- APIs externas (PubChem / ChEBI / ClassyFire) ---
    nome_padronizado: Mapped[str | None] = mapped_column(String(255))
    pubchem_cid: Mapped[int | None]
    chebi_id: Mapped[str | None] = mapped_column(String(50))
    formula: Mapped[str | None] = mapped_column(String(120))
    massa_molecular: Mapped[float | None]                                  # Neutral mass (Da)
    iupac: Mapped[str | None] = mapped_column(Text)
    descricao: Mapped[str | None] = mapped_column(Text)                    # usos/descrição (PubChem)

    # classificação química (ClassyFire / heurística)
    classe_geral: Mapped[str | None] = mapped_column(String(120))          # "Classe geral"
    subclasse: Mapped[str | None] = mapped_column(String(120))             # "Subclasse"
    tipo_composto: Mapped[str | None] = mapped_column(String(30))          # Natural/Sintético
    metabolismo: Mapped[str | None] = mapped_column(String(120))

    criado_em: Mapped[dt.datetime] = mapped_column(DateTime, server_default=func.now())

    # relacionamentos
    metrica: Mapped["MetricaComposto | None"] = relationship(
        back_populates="composto", uselist=False, cascade="all, delete-orphan"
    )
    ranking: Mapped["RankingComposto | None"] = relationship(
        back_populates="composto", uselist=False, cascade="all, delete-orphan"
    )
    abundancias: Mapped[list["Abundancia"]] = relationship(
        back_populates="composto", cascade="all, delete-orphan"
    )
    ontologias: Mapped[list["Ontologia"]] = relationship(
        back_populates="composto", cascade="all, delete-orphan"
    )
    procedencias: Mapped[list["ProcedenciaCampo"]] = relationship(
        back_populates="composto", cascade="all, delete-orphan"
    )
    tags: Mapped[list["Tag"]] = relationship(
        secondary="composto_tag", back_populates="compostos"
    )


# =========================================================
# DIMENSÃO: AMOSTRA
# =========================================================
class Amostra(Base):
    __tablename__ = "amostra"

    id_amostra: Mapped[str] = mapped_column(String(80), primary_key=True)
    nome: Mapped[str] = mapped_column(String(120), nullable=False)         # ex.: Amostra_2
    tipo: Mapped[str | None] = mapped_column(String(50))                   # ex.: branco, óleo NEG
    grupo: Mapped[str | None] = mapped_column(String(80))                  # condição experimental (Grupo de Grupo.Replicata)
    projeto: Mapped[str | None] = mapped_column(String(120))
    data_coleta: Mapped[dt.date | None] = mapped_column(Date)

    abundancias: Mapped[list["Abundancia"]] = relationship(back_populates="amostra")


# =========================================================
# FATO: ABUNDÂNCIA  (composto × amostra × replicata × tipo)
# AJUSTE 3: campo `tipo` (normalised/raw) — o export traz os dois.
# =========================================================
class Abundancia(Base):
    __tablename__ = "abundancia"
    __table_args__ = (
        UniqueConstraint(
            "compound_id", "id_amostra", "replicata", "tipo", name="uq_medicao"
        ),
    )

    id_medicao: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    compound_id: Mapped[str] = mapped_column(
        ForeignKey("composto.compound_id", ondelete="CASCADE")
    )
    id_amostra: Mapped[str] = mapped_column(
        ForeignKey("amostra.id_amostra", ondelete="CASCADE")
    )

    abundancia: Mapped[float] = mapped_column(Float, nullable=False)
    tipo: Mapped[str] = mapped_column(String(12), default="normalised")    # normalised | raw
    replicata: Mapped[int | None]                                          # 2º número de Grupo.Replicata
    is_maxima: Mapped[bool] = mapped_column(Boolean, default=False)        # "amostra mais abundante"

    composto: Mapped["Composto"] = relationship(back_populates="abundancias")
    amostra: Mapped["Amostra"] = relationship(back_populates="abundancias")


# =========================================================
# MÉTRICAS ANALÍTICAS (1:1 com composto) — colunas do Progenesis
# AJUSTE 2: anova_p, q_value, max_fold_change, highest/lowest_mean,
#           maximum_abundance, min_cv_pct, charge, accepted,
#           theoretical_isotope_distribution, link.
# =========================================================
class MetricaComposto(Base):
    __tablename__ = "metrica_composto"

    compound_id: Mapped[str] = mapped_column(
        ForeignKey("composto.compound_id", ondelete="CASCADE"), primary_key=True
    )

    # identificação / qualidade
    accepted: Mapped[bool | None]                                          # "Accepted?"
    modo_aquisicao: Mapped[str | None] = mapped_column(String(50))         # Acquisition Mode / Positivo|Negativo
    adducts: Mapped[str | None] = mapped_column(String(50))
    charge: Mapped[float | None]
    link: Mapped[str | None] = mapped_column(Text)                         # URL base externa (ChemSpider etc.)

    # scores de identificação (ordem de prioridade da ordenação biológica)
    score: Mapped[float | None]
    fragmentation_score: Mapped[float | None]
    isotope_similarity: Mapped[float | None]
    mass_error_ppm: Mapped[float | None]

    # massa / cromatografia
    neutral_mass: Mapped[float | None]
    mz: Mapped[float | None]
    retention_time_min: Mapped[float | None]
    peak_width_min: Mapped[float | None]
    identifications: Mapped[int | None]

    # estatística (Progenesis)
    anova_p: Mapped[float | None]
    q_value: Mapped[float | None]
    max_fold_change: Mapped[float | None]
    highest_mean: Mapped[str | None] = mapped_column(String(80))           # amostra de maior média
    lowest_mean: Mapped[str | None] = mapped_column(String(80))            # amostra de menor média

    # isotopia / abundância
    theoretical_isotope_distribution: Mapped[str | None] = mapped_column(Text)
    isotope_distribution: Mapped[str | None] = mapped_column(Text)
    maximum_abundance: Mapped[float | None]
    min_cv_pct: Mapped[float | None]                                       # Minimum CV%

    # estatísticas derivadas pelo pipeline
    media_abundancia: Mapped[float | None]
    desvio_padrao: Mapped[float | None]
    cv: Mapped[float | None]

    composto: Mapped["Composto"] = relationship(back_populates="metrica")


# =========================================================
# RANKING / ORDENAÇÃO BIOLÓGICA (1:1 com composto)
# Prioridade validada pelo IST: fragmentation_score > score >
# isotope_similarity (>80) > mass_error > fórmula.
# =========================================================
class RankingComposto(Base):
    __tablename__ = "ranking_composto"

    compound_id: Mapped[str] = mapped_column(
        ForeignKey("composto.compound_id", ondelete="CASCADE"), primary_key=True
    )
    confianca: Mapped[float | None]                                        # score ponderado do pipeline
    posicao: Mapped[int | None]                                            # rank global
    rank_group: Mapped[int | None]                                         # campo exigido no PDF
    is_tied: Mapped[bool] = mapped_column(Boolean, default=False)
    feature_group: Mapped[str | None] = mapped_column(String(80))
    original_id: Mapped[str | None] = mapped_column(String(80))

    composto: Mapped["Composto"] = relationship(back_populates="ranking")


# =========================================================
# ONTOLOGIA ChEBI (1:N) — um registro por nível da hierarquia
# =========================================================
class Ontologia(Base):
    __tablename__ = "ontologia"

    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    compound_id: Mapped[str] = mapped_column(
        ForeignKey("composto.compound_id", ondelete="CASCADE")
    )
    nivel: Mapped[int | None]                                              # 0 = mais específico
    termo: Mapped[str] = mapped_column(Text, nullable=False)

    composto: Mapped["Composto"] = relationship(back_populates="ontologias")


# =========================================================
# TAGS (M:N) — configuradas pela equipe no Progenesis (marcadas com 'x')
# AJUSTE 1: confirmado como colunas reais → associação M:N.
# =========================================================
class Tag(Base):
    __tablename__ = "tag"

    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    nome: Mapped[str] = mapped_column(String(60), unique=True)             # "Abund > 1000", "Not Fragmented"...
    descricao: Mapped[str | None] = mapped_column(Text)

    compostos: Mapped[list["Composto"]] = relationship(
        secondary="composto_tag", back_populates="tags"
    )


class CompostoTag(Base):
    __tablename__ = "composto_tag"

    compound_id: Mapped[str] = mapped_column(
        ForeignKey("composto.compound_id", ondelete="CASCADE"), primary_key=True
    )
    tag_id: Mapped[int] = mapped_column(
        ForeignKey("tag.id", ondelete="CASCADE"), primary_key=True
    )


# =========================================================
# RASTREABILIDADE DA ORIGEM (requisito do dashboard)
# =========================================================
class ProcedenciaCampo(Base):
    __tablename__ = "procedencia_campo"

    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    compound_id: Mapped[str] = mapped_column(
        ForeignKey("composto.compound_id", ondelete="CASCADE")
    )
    campo: Mapped[str] = mapped_column(String(60))                         # ex.: "formula", "classe_geral"
    fonte: Mapped[str] = mapped_column(String(40))                         # Progenesis | PubChem | ChEBI | ClassyFire
    referencia: Mapped[str | None] = mapped_column(String(255))            # CID / URL consultada

    composto: Mapped["Composto"] = relationship(back_populates="procedencias")


# =========================================================
# TAGS FIXAS (seed) — conforme PDF de consolidação da reunião
# =========================================================
TAGS_PADRAO: list[tuple[str, str]] = [
    ("Branco", "Indica presença de espaço em branco na amostra."),
    ("Abund > 500", "Abundância acima de 500."),
    ("Abund > 1000", "Abundância acima de 1000."),
    ("Abund > 5000", "Abundância acima de 5000."),
    ("Abund > 10000", "Abundância acima de 10000."),
    ("Anova p-value <= 0.05", "Indica significância estatística."),
    ("Max Fold Change >= 2", "Compostos com Fold Change maior ou igual a 2."),
    ("Not Fragmented", "Indica compostos que não sofreram fragmentação."),
]


def seed_tags(session) -> None:
    """Insere as tags padrão se ainda não existirem."""
    existentes = {nome for (nome,) in session.query(Tag.nome).all()}
    novas = [
        Tag(nome=nome, descricao=desc)
        for nome, desc in TAGS_PADRAO
        if nome not in existentes
    ]
    if novas:
        session.add_all(novas)
        session.commit()


def init_db(url: str = "sqlite:///quimioanalytics.db", echo: bool = False):
    """Cria o schema e retorna uma sessionmaker pronta para uso."""
    engine = create_engine(url, echo=echo)
    Base.metadata.create_all(engine)
    Session = sessionmaker(bind=engine)
    with Session() as s:
        seed_tags(s)
    return Session


if __name__ == "__main__":
    Session = init_db(echo=True)
    print("Banco criado: quimioanalytics.db")
    with Session() as s:
        print("Tags cadastradas:", [t.nome for t in s.query(Tag).all()])
