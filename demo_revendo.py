import streamlit as st
import pandas as pd
from Bio import Entrez
import itertools
import time
from urllib.error import HTTPError, URLError


# ============================================================
# CONFIGURAÇÃO DO PUBMED / NCBI
# ============================================================

Entrez.email = "gouveia.unb@gmail.com"
Entrez.tool = "RevisaoBibliograficaStreamlit"

# Se houver uma API key cadastrada nos Secrets do Streamlit,
# ela será utilizada automaticamente.
try:
    if "NCBI_API_KEY" in st.secrets:
        Entrez.api_key = st.secrets["NCBI_API_KEY"]
except Exception:
    pass

# Número de tentativas automáticas do Biopython
Entrez.max_tries = 3

# Tempo entre tentativas em caso de erro
Entrez.sleep_between_tries = 5


# ============================================================
# CONTROLE DE REQUISIÇÕES
# ============================================================

def wait_between_requests():
    """
    Mantém as requisições abaixo do limite do NCBI.
    Sem API key, utiliza aproximadamente 2 requisições/segundo.
    """
    time.sleep(0.5)


# ============================================================
# GERAR COMBINAÇÕES
# ============================================================

def generate_combinations(terms, min_terms=1):
    """
    Gera todas as combinações possíveis de termos
    utilizando o operador AND.
    """

    combinations = []

    # Remove espaços e termos vazios
    clean_terms = [
        term.strip()
        for term in terms
        if term.strip()
    ]

    for i in range(min_terms, len(clean_terms) + 1):
        for subset in itertools.combinations(clean_terms, i):
            combinations.append(
                " AND ".join(
                    [f'"{term}"' for term in subset]
                )
            )

    return combinations


# ============================================================
# BUSCAR ARTIGOS NO PUBMED
# ============================================================

@st.cache_data(show_spinner=False, ttl=3600)
def fetch_article_counts(term, year, max_results=3):

    query = f"{term} AND {year}[DP]"

    for attempt in range(3):

        try:

            wait_between_requests()

            handle = Entrez.esearch(
                db="pubmed",
                term=query,
                retmax=max_results,
                usehistory="n"
            )

            record = Entrez.read(handle)
            handle.close()

            count = int(record["Count"])
            ids = [str(x) for x in record["IdList"]]

            return count, ids

        except HTTPError as e:

            if attempt < 2:

                # 429 = excesso de requisições
                # 403 = acesso recusado
                # 500/502/503 etc. = erro temporário do servidor

                if e.code in [403, 429, 500, 502, 503, 504]:

                    time.sleep(5 * (attempt + 1))
                    continue

            st.error(
                f"Erro ao consultar o PubMed "
                f"(HTTP {e.code}). "
                f"Tente novamente em alguns minutos."
            )

            return 0, []

        except URLError:

            if attempt < 2:
                time.sleep(5 * (attempt + 1))
                continue

            st.error(
                "Não foi possível estabelecer conexão "
                "com o PubMed."
            )

            return 0, []

        except Exception as e:

            st.error(
                f"Erro inesperado ao consultar o PubMed: "
                f"{type(e).__name__}: {e}"
            )

            return 0, []

    return 0, []


# ============================================================
# BUSCAR DETALHES DOS ARTIGOS
# ============================================================

@st.cache_data(show_spinner=False, ttl=3600)
def fetch_article_details(ids):

    if not ids:
        return []

    for attempt in range(3):

        try:

            wait_between_requests()

            handle = Entrez.efetch(
                db="pubmed",
                id=",".join(ids),
                rettype="xml",
                retmode="text"
            )

            records = Entrez.read(handle)
            handle.close()

            articles = []

            for article in records["PubmedArticle"]:

                # ------------------------------------------------
                # PMID
                # ------------------------------------------------

                pubmed_id = str(
                    article["MedlineCitation"]["PMID"]
                )

                # ------------------------------------------------
                # Título
                # ------------------------------------------------

                title = str(
                    article["MedlineCitation"]["Article"]["ArticleTitle"]
                )

                # ------------------------------------------------
                # Revista
                # ------------------------------------------------

                journal = str(
                    article["MedlineCitation"]
                    ["Article"]
                    ["Journal"]
                    ["Title"]
                )

                # ------------------------------------------------
                # Ano
                # ------------------------------------------------

                pub_date = (
                    article["MedlineCitation"]
                    ["Article"]
                    ["Journal"]
                    ["JournalIssue"]
                    ["PubDate"]
                )

                year = pub_date.get("Year", "N/A")

                # Caso o PubMed não tenha Year, tenta
                # recuperar MedlineDate.
                if year == "N/A":
                    year = pub_date.get(
                        "MedlineDate",
                        "N/A"
                    )

                # ------------------------------------------------
                # Keywords
                # ------------------------------------------------

                keywords_list = article[
                    "MedlineCitation"
                ].get("KeywordList", [])

                if keywords_list:

                    keywords = ", ".join(
                        str(keyword)
                        for keyword_group in keywords_list
                        for keyword in keyword_group
                    )

                else:

                    keywords = "N/A"

                # ------------------------------------------------
                # DOI
                # ------------------------------------------------

                article_ids = (
                    article["PubmedData"]
                    .get("ArticleIdList", [])
                )

                doi = "N/A"

                for article_id in article_ids:

                    try:

                        if (
                            article_id.attributes
                            .get("IdType") == "doi"
                        ):
                            doi = str(article_id)
                            break

                    except Exception:
                        continue

                if doi != "N/A":

                    doi_link = (
                        f"https://doi.org/{doi}"
                    )

                else:

                    doi_link = "N/A"

                # ------------------------------------------------
                # Adiciona artigo
                # ------------------------------------------------

                articles.append({

                    "PubMed ID": pubmed_id,

                    "Title": title,

                    "Year": year,

                    "Journal": journal,

                    "Keywords": keywords,

                    "DOI": doi_link
                })

            return articles

        except HTTPError as e:

            if attempt < 2:

                if e.code in [403, 429, 500, 502, 503, 504]:

                    time.sleep(5 * (attempt + 1))
                    continue

            st.error(
                f"Erro ao buscar detalhes dos artigos "
                f"(HTTP {e.code})."
            )

            return []

        except URLError:

            if attempt < 2:
                time.sleep(5 * (attempt + 1))
                continue

            st.error(
                "Não foi possível estabelecer conexão "
                "com o PubMed."
            )

            return []

        except Exception as e:

            st.error(
                f"Erro ao processar artigos: "
                f"{type(e).__name__}: {e}"
            )

            return []

    return []


# ============================================================
# DOWNLOAD TXT
# ============================================================

def download_as_txt(df, filename):

    txt_data = df.to_string(index=False)

    st.download_button(
        label=f"Baixar {filename} (TXT)",
        data=txt_data.encode("utf-8"),
        file_name=f"{filename}.txt",
        mime="text/plain",
        key=f"{filename}_txt"
    )


# ============================================================
# DOWNLOAD XLSX
# ============================================================

def download_as_xlsx(df, filename):

    from io import BytesIO

    output = BytesIO()

    with pd.ExcelWriter(
        output,
        engine="openpyxl"
    ) as writer:

        df.to_excel(
            writer,
            index=False,
            sheet_name="Resultados"
        )

    output.seek(0)

    st.download_button(
        label=f"Baixar {filename} (XLSX)",
        data=output,
        file_name=f"{filename}.xlsx",
        mime=(
            "application/vnd.openxmlformats-officedocument."
            "spreadsheetml.sheet"
        ),
        key=f"{filename}_xlsx"
    )


# ============================================================
# APLICAÇÃO PRINCIPAL
# ============================================================

def main():

    st.title(
        "Revisão Bibliográfica do PubMed"
    )

    st.write(
        "Informe os termos de pesquisa separados por vírgula."
    )

    # --------------------------------------------------------
    # Entradas
    # --------------------------------------------------------

    terms_input = st.text_input(
        "Termos de pesquisa:",
        "machine learning, drug discovery, cancer"
    )

    start_year = st.number_input(
        "Ano de início:",
        min_value=1900,
        max_value=2100,
        value=2020
    )

    end_year = st.number_input(
        "Ano de fim:",
        min_value=1900,
        max_value=2100,
        value=2024
    )

    max_articles = st.selectbox(
        "Número máximo de artigos relevantes:",
        [10, 25, 50]
    )

    # --------------------------------------------------------
    # Botão de pesquisa
    # --------------------------------------------------------

    if st.button(
        "Enviar pesquisa",
        type="primary"
    ):

        if start_year > end_year:

            st.error(
                "O ano de início deve ser menor "
                "ou igual ao ano de fim."
            )

            return

        terms = [
            term.strip()
            for term in terms_input.split(",")
            if term.strip()
        ]

        if not terms:

            st.error(
                "Informe pelo menos um termo de pesquisa."
            )

            return

        # Limpa resultados anteriores
        st.session_state.pop(
            "df_single_term_or_more",
            None
        )

        st.session_state.pop(
            "df_relevant",
            None
        )

        st.session_state["terms"] = terms
        st.session_state["start_year"] = start_year
        st.session_state["end_year"] = end_year
        st.session_state["max_articles"] = max_articles

        # ====================================================
        # TABELA 1
        # ====================================================

        st.subheader(
            "1. Quantidade de Artigos por Ano"
        )

        combinations = generate_combinations(
            terms,
            min_terms=1
        )

        results = []

        progress = st.progress(0)

        total_queries = (
            len(combinations)
            * (end_year - start_year + 1)
        )

        current_query = 0

        for combo in combinations:

            row = {
                "combination": combo
            }

            for year in range(
                start_year,
                end_year + 1
            ):

                count, _ = fetch_article_counts(
                    combo,
                    year,
                    max_results=3
                )

                row[year] = count

                current_query += 1

                progress.progress(
                    min(
                        current_query / total_queries,
                        1.0
                    )
                )

            results.append(row)

        progress.empty()

        df = pd.DataFrame(results)

        st.session_state[
            "df_single_term_or_more"
        ] = df

        st.dataframe(
            df,
            use_container_width=True
        )

        download_as_txt(
            df,
            "quantidade_artigos_por_ano"
        )

        download_as_xlsx(
            df,
            "quantidade_artigos_por_ano"
        )

    # ========================================================
    # MOSTRAR TABELA 1 SE JÁ EXISTIR
    # ========================================================

    elif "df_single_term_or_more" in st.session_state:

        st.subheader(
            "1. Quantidade de Artigos por Ano"
        )

        df = st.session_state[
            "df_single_term_or_more"
        ]

        st.dataframe(
            df,
            use_container_width=True
        )

        download_as_txt(
            df,
            "quantidade_artigos_por_ano"
        )

        download_as_xlsx(
            df,
            "quantidade_artigos_por_ano"
        )

    # ========================================================
    # TABELA 2
    # ========================================================

    if (
        "terms" in st.session_state
        and "df_single_term_or_more" in st.session_state
        and "df_relevant" not in st.session_state
    ):

        st.subheader(
            "2. Artigos Mais Relevantes"
        )

        terms = st.session_state["terms"]
        start_year = st.session_state["start_year"]
        end_year = st.session_state["end_year"]
        max_articles = st.session_state["max_articles"]

        detailed_articles = []

        # ----------------------------------------------------
        # Prioriza combinações com maior número de termos
        # ----------------------------------------------------

        stop_search = False

        for num_terms in range(
            len(terms),
            0,
            -1
        ):

            combinations = generate_combinations(
                terms,
                min_terms=num_terms
            )

            for combo in combinations:

                if stop_search:
                    break

                for year in range(
                    start_year,
                    end_year + 1
                ):

                    if len(detailed_articles) >= max_articles:
                        stop_search = True
                        break

                    # Busca IDs
                    _, ids = fetch_article_counts(
                        combo,
                        year,
                        max_results=max_articles
                    )

                    if not ids:
                        continue

                    # Limita a quantidade de IDs
                    remaining = (
                        max_articles
                        - len(detailed_articles)
                    )

                    ids_to_fetch = ids[:remaining]

                    articles = fetch_article_details(
                        tuple(ids_to_fetch)
                    )

                    detailed_articles.extend(
                        articles
                    )

                if stop_search:
                    break

            if stop_search:
                break

        # ----------------------------------------------------
        # Remove duplicados
        # ----------------------------------------------------

        unique_articles = {}

        for article in detailed_articles:

            pmid = article["PubMed ID"]

            unique_articles[pmid] = article

        detailed_articles = list(
            unique_articles.values()
        )

        detailed_articles = detailed_articles[
            :max_articles
        ]

        df_relevant = pd.DataFrame(
            detailed_articles
        )

        st.session_state[
            "df_relevant"
        ] = df_relevant

        if not df_relevant.empty:

            st.dataframe(
                df_relevant,
                use_container_width=True
            )

            download_as_txt(
                df_relevant,
                "artigos_mais_relevantes"
            )

            download_as_xlsx(
                df_relevant,
                "artigos_mais_relevantes"
            )

        else:

            st.warning(
                "Nenhum artigo foi encontrado "
                "para os critérios informados."
            )


# ============================================================
# EXECUÇÃO
# ============================================================

if __name__ == "__main__":
    main()
