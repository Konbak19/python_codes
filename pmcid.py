import spacy
from bs4 import BeautifulSoup
import pandas as pd
import xlsxwriter

# Load spaCy model for Named Entity Recognition
nlp = spacy.load("en_core_web_sm")

# Function to perform Named Entity Recognition on a given text
def perform_ner(text):
    doc = nlp(text)
    entities = [(ent.text, ent.label_) for ent in doc.ents]
    return entities

# Opening and parsing the XML file
with open('articles.xml', 'r', encoding='utf-8') as xml_file:
    soup = BeautifulSoup(xml_file, 'lxml')

# Create an Excel workbook and add a worksheet
workbook = xlsxwriter.Workbook('ner_results.xlsx')
worksheet = workbook.add_worksheet()

# Define headers for the Excel file
headers = ['PMC ID']

# Initialize lists to store results
pmcids = []
titles = []
ner_results = []

# Write the headers to the first row of the worksheet
for col, header in enumerate(headers):
    worksheet.write(0, col, header)

try:
    # Collecting all the articles
    all_articles = soup.find_all("articleid", attrs={"idtype": "pmc"})

    # Start row counter
    row = 1

    for article in all_articles:
        # Extracting PMC ID
        pmcid = article.text
        pmcids.append(pmcid)

        # Performing Named Entity Recognition on the full text
        # You can add the NER processing here if needed

        # Write data to Excel
        worksheet.write(row, 0, pmcid)

        # Increment row counter
        row += 1

    # Close the Excel workbook
    workbook.close()
    print("Excel file successfully created.")

except Exception as e:
    print(f"An error occurred: {str(e)}")

print("Results saved to ner_results.xlsx.")
