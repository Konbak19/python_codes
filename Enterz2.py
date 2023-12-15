from Bio import Entrez
from bs4 import BeautifulSoup
import xlsxwriter
import xml.etree.ElementTree as ET
import os

# Prompt the user for search terms
search_terms = input("Enter search terms (e.g., cancer mutation gene): ")

# Identifying the connected user
Entrez.email = 'xryspa19@mail.com'
Entrez.tool = 'Demoscript'

# Construct the search query using user input
query = f'({search_terms}) AND ("mutation"[Title/Abstract] OR "gene"[Title/Abstract])'

info = Entrez.esearch(db="pubmed", retmax=10, term=query)

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
headers = ['PubMed ID', 'Title', 'Abstract']

def pmid2pmcid(email, pmid):
    Entrez.email = email

    handle = Entrez.elink(dbfrom="pubmed", db="pmc", linkname="pubmed_pmc", id=pmid, retmode="text")

    handle_read = handle.read()
    handle.close()

    root = ET.fromstring(handle_read)

    pmcid = ""

    for link in root.iter('Link'):
        for id in link.iter('Id'):
            pmcid = id.text
    return pmcid

# Write the headers to the first row of the worksheet
for col, header in enumerate(headers):
    worksheet.write(0, col, header)

# Check if the 'abstract_txt' folder exists, if not, create it
abstract_folder_path = 'abstract_txt'
if not os.path.exists(abstract_folder_path):
    os.makedirs(abstract_folder_path)

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

            # Searching for the article's abstract
            abstract = article.abstracttext.text if article.abstracttext else "N/A"

            # Write data to Excel
            worksheet.write(row, 0, pmid)
            worksheet.write(row, 1, title)
            worksheet.write(row, 2, abstract)

            # Save the abstract to a text file
            abstract_file_path = os.path.join(abstract_folder_path, f"{pmid}_abstract.txt")
            with open(abstract_file_path, 'w', encoding='utf-8') as abstract_file:
                abstract_file.write(f"PubMed ID: {pmid}\nTitle: {title}\nAbstract:\n{abstract}\n")

            # Increment row counter
            row += 1

    # Close the Excel workbook
    workbook.close()
    print("Excel file successfully created.")

except Exception as e:
    print(f"An error occurred: {str(e)}")
