import streamlit as st
import pandas as pd
from Bio import Entrez
import itertools

# Configuração do e-mail do Entrez
Entrez.email = "gouveia.unb@gmail.com"  # Substitua pelo seu e-mail

# Função para gerar combinações de termos (com N termos ou menos)
def generate_combinations(terms, min_terms=1):
    """Gera todas as combinações possíveis de termos com operadores booleanos."""
    combinations = []
    for i in range(min_terms, len(terms) + 1):
        for subset in itertools.combinations(terms, i):
            combinations.append(' AND '.join([f'"{term.strip()}"' for term in subset]))
    return combinations

# Função para buscar contagem de artigos no PubMed
def fetch_article_counts(term, year, max_results=3):
    query = f"{term} AND {year}[DP]"
    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results, usehistory="y")
    record = Entrez.read(handle)
    handle.close()
    return int(record["Count"]), record["IdList"]

# Função para buscar detalhes dos artigos no PubMed
def fetch_article_details(ids):
    handle = Entrez.efetch(db="pubmed", id=",".join(ids), rettype="xml", retmode="text")
    records = Entrez.read(handle)
    handle.close()
    
    articles = []
    for article in records["PubmedArticle"]:
        pubmed_id = article["MedlineCitation"]["PMID"]
        title = article["MedlineCitation"]["Article"]["ArticleTitle"]
        journal = article["MedlineCitation"]["Article"]["Journal"]["Title"]
        year = article["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["PubDate"].get("Year", "N/A")
        
        # Extraindo keywords (se disponíveis)
        keywords_list = article["MedlineCitation"].get("KeywordList", [])
        keywords = ", ".join([kw for kw_group in keywords_list for kw in kw_group]) if keywords_list else "N/A"

        # Extraindo DOI
        article_ids = article["PubmedData"].get("ArticleIdList", [])
        doi = next((id_ for id_ in article_ids if id_.attributes["IdType"] == "doi"), "N/A")
        doi_link = f"[{doi}](https://doi.org/{doi})" if doi != "N/A" else "N/A"

        articles.append({
            "PubMed ID": pubmed_id,
            "Title": title,
            "Year": year,
            "Journal": journal,
            "Keywords": keywords,
            "DOI": doi_link
        })
    return articles

# Função principal de Streamlit
def main():
    st.title("Revisão Bibliográfica do PubMed")

    # Entrada do usuário
    terms_input = st.text_input("Insira os termos de pesquisa separados por vírgula:", "machine learning, drug discovery, cancer")
    start_year = st.number_input("Ano de início:", min_value=1900, max_value=2100, value=2020)
    end_year = st.number_input("Ano de fim:", min_value=1900, max_value=2100, value=2024)
    max_articles = st.selectbox("Número máximo de artigos relevantes:", [10, 25, 50])

    if st.button("Enviar pesquisa"):
        terms = terms_input.split(",")

        # Tabela 1: Quantidade de artigos por ano
        combinations_single_term_or_more = generate_combinations(terms, min_terms=1)
        results_single_term_or_more = []

        for combo in combinations_single_term_or_more:
            row = {"combination": combo}
            for year in range(start_year, end_year + 1):
                count, _ = fetch_article_counts(combo, year)
                row[year] = count
            results_single_term_or_more.append(row)

        st.subheader("Quantidade de Artigos por Ano (Combinações de Termos)")
        df_single_term_or_more = pd.DataFrame(results_single_term_or_more)
        st.write(df_single_term_or_more)

        # Tabela 2: Artigos mais relevantes
        detailed_articles = []

        # Priorizar combinações com mais termos
        for num_terms in range(len(terms), 0, -1):  # De N para 1
            combinations_multiple_terms = generate_combinations(terms, min_terms=num_terms)
            for combo in combinations_multiple_terms:
                for year in range(start_year, end_year + 1):
                    _, ids = fetch_article_counts(combo, year, max_results=max_articles)
                    if ids:
                        articles = fetch_article_details(ids[:max_articles])
                        detailed_articles.extend(articles)

                        # Interrompe a busca ao atingir o limite de artigos
                        if len(detailed_articles) >= max_articles:
                            break
                if len(detailed_articles) >= max_articles:
                    break
            if len(detailed_articles) >= max_articles:
                break

        if detailed_articles:
            st.subheader("Artigos Mais Relevantes (Prioridade para Combinações com N Termos)")
            df_relevant = pd.DataFrame(detailed_articles)
            st.write(df_relevant)

            # Baixar tabelas relevantes
            def download_as_txt(df, filename):
                txt_data = df.to_string(index=False)
                st.download_button(
                    label=f"Baixar {filename} (TXT)",
                    data=txt_data.encode('utf-8'),
                    file_name=f"{filename}.txt",
                    mime="text/plain",
                )

            def download_as_xlsx(df, filename):
                xlsx_file = pd.ExcelWriter(f"{filename}.xlsx", engine="openpyxl")
                df.to_excel(xlsx_file, index=False, sheet_name="Sheet1")
                xlsx_file.close()
                with open(f"{filename}.xlsx", "rb") as f:
                    xlsx_data = f.read()
                st.download_button(
                    label=f"Baixar {filename} (XLSX)",
                    data=xlsx_data,
                    file_name=f"{filename}.xlsx",
                    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                )

            download_as_txt(df_relevant, "artigos_mais_relevantes")
            download_as_xlsx(df_relevant, "artigos_mais_relevantes")

        else:
            st.write("Nenhum artigo relevante encontrado para os termos fornecidos.")

if __name__ == "__main__":
    main()
