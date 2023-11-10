from Bio import Entrez, SeqIO
import xlsxwriter

Entrez.email = ""

def search_records():
    search_type = input("Do you want to search by Gene Name (G), Organism (O), or both (B)? ").upper()

    if search_type == "G":
        gene_name = input("Enter the Gene Name: ")
        query = f'{gene_name}[Gene Name]'
    elif search_type == "O":
        organism = input("Enter the Organism: ")
        query = f'"{organism}"[Organism]'
    elif search_type == "B":
        gene_name = input("Enter the Gene Name: ")
        organism = input("Enter the Organism: ")
        query = f'{gene_name}[Gene Name] AND "{organism}"[Organism]'
    else:
        print("Invalid choice. Please enter G, O, or B.")
        return

    handle = Entrez.esearch(db="nucleotide", term=query, retmax="40")
    rec_list = Entrez.read(handle)
    handle.close()

    if rec_list['Count'] == "0":
        print(f"No records found for {gene_name} in {organism}.")
        return

    id_list = rec_list['IdList']
    handle = Entrez.efetch(db='nucleotide', id=id_list, rettype='gb')

    recs = list(SeqIO.parse(handle, 'gb'))
    handle.close()

    workbook = xlsxwriter.Workbook('search_results.xlsx')
    worksheet = workbook.add_worksheet()

    headers = ['Name', 'Description']

    for col, header in enumerate(headers):
        worksheet.write(0, col, header)

    for row, rec in enumerate(recs, 1):
        worksheet.write(row, 0, rec.name)
        worksheet.write(row, 1, rec.description)

    workbook.close()
    print("Search results saved to search_results.xlsx")

search_records()



