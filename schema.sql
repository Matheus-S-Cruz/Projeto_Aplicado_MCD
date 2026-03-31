-- =========================
-- DIMENSÃO: COMPOSTOS
-- =========================
CREATE TABLE dim_composto (
    compound_id VARCHAR(50) PRIMARY KEY,
    nome_composto VARCHAR(255) NOT NULL,
    descricao TEXT,
    pubchem_cid INT,
    chebi_id VARCHAR(50),
    formula VARCHAR(100),
    massa_api DOUBLE PRECISION, -- 🔥 ADICIONADO (obrigatório)
    role_biologico VARCHAR(100),
    classe_taxonomica VARCHAR(100)
);

-- =========================
-- DIMENSÃO: AMOSTRAS
-- =========================
CREATE TABLE dim_amostra (
    id_amostra VARCHAR(50) PRIMARY KEY,
    nome_amostra VARCHAR(100) NOT NULL,
    tipo VARCHAR(50),
    projeto VARCHAR(100),
    data_coleta DATE
);

-- =========================
-- TABELA FATO: ABUNDÂNCIA
-- =========================
CREATE TABLE fact_abundancia (
    id_medicao SERIAL PRIMARY KEY,
    
    compound_id VARCHAR(50) NOT NULL,
    id_amostra VARCHAR(50) NOT NULL,

    abundancia_valor DOUBLE PRECISION NOT NULL,
    sinal_analitico VARCHAR(255),

    replicata INT, -- 🔥 MELHORADO (antes era texto)

    score_identificacao DOUBLE PRECISION,

    -- 🔗 CHAVES ESTRANGEIRAS
    FOREIGN KEY (compound_id)
        REFERENCES dim_composto(compound_id)
        ON DELETE CASCADE
        ON UPDATE CASCADE,

    FOREIGN KEY (id_amostra)
        REFERENCES dim_amostra(id_amostra)
        ON DELETE CASCADE
        ON UPDATE CASCADE,

    -- 🔒 EVITA DUPLICIDADE DE MEDIÇÃO
    CONSTRAINT unique_medicao
    UNIQUE (compound_id, id_amostra, replicata)
);

-- =========================
-- DIMENSÃO: ONTOLOGIA
-- =========================
CREATE TABLE dim_ontologia (
    id SERIAL PRIMARY KEY,
    compound_id VARCHAR(50) NOT NULL,
    ontologia TEXT NOT NULL,

    FOREIGN KEY (compound_id)
        REFERENCES dim_composto(compound_id)
        ON DELETE CASCADE
        ON UPDATE CASCADE
);

-- =========================
-- ÍNDICES (PERFORMANCE)
-- =========================
CREATE INDEX idx_fact_compound 
ON fact_abundancia(compound_id);

CREATE INDEX idx_fact_amostra 
ON fact_abundancia(id_amostra);

CREATE INDEX idx_fact_composto_amostra 
ON fact_abundancia(compound_id, id_amostra);

CREATE INDEX idx_ontologia_compound 
ON dim_ontologia(compound_id);

CREATE INDEX idx_composto_id 
ON dim_composto(compound_id);

-- =========================
-- VIEW (OBRIGATÓRIA NA AULA)
-- =========================
CREATE VIEW v_analito_info AS
SELECT
    nome_composto AS analito,
    massa_api AS massa_molecular,
    role_biologico AS uso_conhecido,
    formula,
    pubchem_cid,
    chebi_id
FROM dim_composto;