from Bio import Entrez
import xlsxwriter

# Identifying the connected user
Entrez.email = 'janedoe@mail.com'
Entrez.tool = 'Demoscript'

# Defining the search query
query = '("cancer"[Title/Abstract])'

# Retrieve PubMed IDs for articles matching the query
info = Entrez.esearch(db="pubmed", retmax=10, term=query)
record = Entrez.read(info)
pmids = record['IdList']

# Create an Excel workbook and add a worksheet
workbook = xlsxwriter.Workbook('pubmed_articles2.xlsx')
worksheet = workbook.add_worksheet()

# Define headers for the Excel file
headers = ['PubMed ID', 'Title']

# Write the headers to the first row of the worksheet
for col, header in enumerate(headers):
    worksheet.write(0, col, header)

# Fetch and write article data to the Excel file
row = 1
for pmid in pmids:
    article_info = Entrez.efetch(db='pubmed', id=pmid, retmode='xml')
    article_record = Entrez.read(article_info)
    
    # Extract PubMed ID and title
    pmid = article_record['PubmedArticle'][0]['PubmedData']['ArticleIdList'][0]
    title = article_record['PubmedArticle'][0]['Article']['ArticleTitle']
    
    # Write data to Excel
    worksheet.write(row, 0, pmid)
    worksheet.write(row, 1, title)
    
    row += 1

# Close the Excel workbook
workbook.close()
print("Excel file successfully created.")


