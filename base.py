from Bio import Entrez
from bs4 import BeautifulSoup
import xlsxwriter
import scispacy
import spacy

# Prompt the user for search terms
#search_terms = input("Enter search terms (e.g., cancer mutation gene): ")

# Identifying the connected user
Entrez.email = 'xryspa19@mail.com'
Entrez.tool = 'Demoscript'

# Construct the search query using user input
query = '("alzheimer"[Title/Abstract] OR "alzheimers"[Title/Abstract] OR "parkinson"[Title/Abstract] OR "schizophrenia"[Title/Abstract]) AND ("variant"[Title/Abstract] OR "mutation"[Title/Abstract] OR "polymorphism"[Title/Abstract] OR "snp"[Title/Abstract] OR "variation"[Title/Abstract] OR "mutant"[Title/Abstract]) AND ("microrna"[Title/Abstract] OR "mirna"[Title/Abstract]) NOT Review[ptyp]'
info = Entrez.esearch(db="pubmed", retmax=1000, term=query)
# parsing the XML data
record = Entrez.read(info)
# retrieving records in XML format
fetch = Entrez.efetch(db='pubmed', resetmode='xml', id=record['IdList'], rettype='', retmode='xml')
# writing records in XML file
with open('articles.xml', 'wb') as f:
    f.write(fetch.read())

# Create an Excel workbook and add a worksheet
workbook = xlsxwriter.Workbook('pubmed_articles.xlsx')
worksheet = workbook.add_worksheet()

# Define headers for the Excel file
headers = ['PubMed ID', 'Title']

# Write the headers to the first row of the worksheet
for col, header in enumerate(headers):
    worksheet.write(0, col, header)

try:
    # Open and parse the XML content
    with open('articles.xml', "r") as xml_file:
        soup = BeautifulSoup(xml_file, 'lxml')
        all_articles = soup.findAll('pubmedarticle')

        # Start row counter
        row = 1

        # Finding the information for every article
        for article in all_articles:
            # Searching for PubMed ID
            pmid = article.pmid.text

            # Searching for the article's title
            title = article.articletitle.text

            # Write data to Excel
            worksheet.write(row, 0, pmid)
            worksheet.write(row, 1, title)

            # Increment row counter
            row += 1

    # Close the Excel workbook
    workbook.close()
    print("Excel file successfully created.")

except Exception as e:
    print(f"An error occurred: {str(e)}")
    
nlp = spacy.load("en_ner_bionlp13cg_md")

text="
