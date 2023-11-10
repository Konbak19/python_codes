import spacy
from bs4 import BeautifulSoup
import xlsxwriter

# Load the spaCy model (you may need to install spaCy and download the language model)
nlp = spacy.load("en_ner_bionlp13cg_md")

# Create a dictionary to store categories
categories = {
    "cancer": [],
    "mutation": [],
    "gene": [],
    "other": [],
}

# Open and parse the XML file
with open('articles.xml', "r") as xml_file:
    soup = BeautifulSoup(xml_file, 'xml')
    all_articles = soup.findAll('pubmedarticle')

    # Process and categorize each article
    for article in all_articles:
        # Get the article's title
        title = article.articletitle.text

        # Print the title for debugging
        print(f"Title: {title}")

        # Create a spaCy document from the title
        doc = nlp(title)

        # Check for keywords in the title and categorize
        if any(token.text.lower() == "cancer" for token in doc):
            categories["cancer"].append(title)
        elif any(token.text.lower() == "mutation" for token in doc):
            categories["mutation"].append(title)
        elif any(token.text.lower() == "gene" for token in doc):
            categories["gene"].append(title)
        else:
            categories["other"].append(title)

# Create Excel files for each category
for category, articles in categories.items():
    if articles:
        # Create an Excel workbook
        workbook = xlsxwriter.Workbook(f'{category}_articles.xlsx')
        worksheet = workbook.add_worksheet()

        # Define headers for the Excel file
        headers = ['Title']

        # Write the headers to the first row of the worksheet
        for col, header in enumerate(headers):
            worksheet.write(0, col, header)

        # Write the articles to the Excel file
        for row, article in enumerate(articles, 1):
            worksheet.write(row, 0, article)

        # Close the Excel workbook
        workbook.close()

# Print a message for each category
for category, articles in categories.items():
    if articles:
        print(f"{len(articles)} {category} articles saved to {category}_articles.xlsx")

print("Excel files successfully created.")

