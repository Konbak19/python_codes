import os
import urllib.request
import xlsxwriter
from Bio import Entrez

def get_assembly_summary(id):
    """Get esummary for an entrez id"""
    esummary_handle = Entrez.esummary(db="assembly", id=id, report="full")
    esummary_record = Entrez.read(esummary_handle)
    return esummary_record

def get_assemblies(term, download=True, path='assemblies'):
    """Download genbank assemblies for a given search term.
    Args:
        term: search term, usually organism name
        download: whether to download the results
        path: folder to save to
    """

    # Provide your own mail here
    Entrez.email = "A.N.Other@example.com"
    handle = Entrez.esearch(db="assembly", term=term, retmax='10')
    record = Entrez.read(handle)
    ids = record['IdList']
    print(f'found {len(ids)} ids')
    links = []
    for id in ids:
        # Get summary
        summary = get_assembly_summary(id)
        # Get FTP link
        url = summary['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_RefSeq']
        if url == '':
            continue
        label = os.path.basename(url)
        # Get the fasta link - change this to get other formats
        link = os.path.join(url, f'{label}_genomic.fna.gz')
        print(link)
        links.append(link)
        if download:
            # Download link
            urllib.request.urlretrieve(link, f'{label}.fna.gz')

    return ids, links

# Get results
ids, links = get_assemblies("mycobacterium tuberculosis", download=True)

# Save results to Excel
workbook = xlsxwriter.Workbook('assembly_results.xlsx')
worksheet = workbook.add_worksheet()

# Define headers for the Excel file
headers = ['ID', 'Link']

# Write the headers to the first row of the worksheet
for col, header in enumerate(headers):
    worksheet.write(0, col, header)

# Write the results to the Excel file
for row,id in enumerate(ids, 1):
    worksheet.write(row, 0, id)
    #worksheet.write(row, 1, links[0])

# Close the Excel workbook
workbook.close()

print("Excel file successfully created.")



