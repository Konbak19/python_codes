from Bio import Entrez
from bs4 import BeautifulSoup
import xlsxwriter
import xml.etree.ElementTree as ET
import os

# Prompt the user for gene and disease lists (comma-separated)
gene_list_str = input("Enter a comma-separated list of genes: ")
disease_list_str = input("Enter a comma-separated list of diseases: ")

# Convert strings to lists
gene_list = gene_list_str.split(",")
disease_list = disease_list_str.split(",")

# Identifying the connected user
Entrez.email = 'xryspa19@mail.com'  # Replace with your actual email
Entrez.tool = 'Demoscript'

# Create a folder for comparison results (if it doesn't exist)
comparison_folder = 'gene_disease_comparisons'
if not os.path.exists(comparison_folder):
    os.makedirs(comparison_folder)

# Loop through each gene in the list
for gene in gene_list:
    gene = gene.strip()  # Remove any leading/trailing whitespaces
    print(gene)

    # Loop through each disease in the list
    for disease in disease_list:
        disease = disease.strip()

        # Construct search query (gene AND disease)
        query = f'({gene}[Title/Abstract] AND {disease} [Title/Abstract])'

        # Include sort parameter for "relevance"
        info = Entrez.esearch(db="pubmed", retmax=10, term=query, sort="relevance")

        # Parsing the XML data
        record = Entrez.read(info)

        # Retrieving records in XML format
        fetch = Entrez.efetch(db='pubmed', resetmode='xml', id=record['IdList'], rettype='', retmode='xml')

        # Writing records in XML file
        with open('articles.xml', 'wb') as f:
            f.write(fetch.read())

        # Create a filename for the Excel file
        filename = f"{comparison_folder}/{gene} vs {disease}.xlsx"

        # Create an Excel workbook
        workbook = xlsxwriter.Workbook(filename)
        worksheet = workbook.add_worksheet()

        # Define headers for the Excel sheet
        headers = ['PubMed ID', 'Title', 'Abstract']

        # Write headers to the worksheet
        for col, header in enumerate(headers):
            worksheet.write(0, col, header)

        # Check if 'abstract_txt' folder exists, if not, create it
        abstract_folder_path = 'abstract_txt'
        if not os.path.exists(abstract_folder_path):
            os.makedirs(abstract_folder_path)

        try:
            # Parse XML content and write data to Excel/text file
            with open('articles.xml', "r") as xml_file:
                soup = BeautifulSoup(xml_file, 'lxml')
                all_articles = soup.findAll('pubmedarticle')

                row = 1

                for article in all_articles:
                    pmid = article.pmid.text
                    title = article.articletitle.text
                    abstract = article.abstracttext.text if article.abstracttext else "N/A"

                    worksheet.write(row, 0, pmid)
                    worksheet.write(row, 1, title)
                    worksheet.write(row, 2, abstract)

                    # Create a separate folder for abstracts of this query
                    abstract_folder_path = os.path.join('abstracts_txt', f"{gene}_vs_{disease}")
                    if not os.path.exists(abstract_folder_path):
                        os.makedirs(abstract_folder_path)

                    abstract_file_path = os.path.join(abstract_folder_path, f"{pmid}_abstract.txt")
                    with open(abstract_file_path, 'w', encoding='utf-8') as abstract_file:
                        abstract_file.write(f"{abstract}")

                    row += 1

        except Exception as e:
            print(f"An error occurred for {gene} vs {disease}: {str(e)}")

        # Close the Excel workbook
        workbook.close()

print("Excel files with comparison results created in the 'gene_disease_comparisons' folder.")
print("Abstracts are stored in separate folders within the 'abstracts_txt' directory.")
