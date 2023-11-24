import spacy
from collections import Counter
import pandas as pd
from spacy import displacy
import os
import fitz  # Use fitz if PyMuPDF is installed
from sklearn.feature_extraction.text import TfidfVectorizer

# Load the spaCy model
nlp = spacy.load("en_core_sci_sm")

# Specify the PDF file path (change this to your PDF file)
pdf_path = '/home/konstantinos/nihms957693.pdf'

# Check if the file exists
if not os.path.isfile(pdf_path):
    print(f"The file '{pdf_path}' does not exist.")
else:
    # Extract text from the PDF using PyMuPDF (fitz)
    pdf_text = ""
    with fitz.open(pdf_path) as pdf_document:
        for page_number in range(pdf_document.page_count):
            page = pdf_document.load_page(page_number)
            pdf_text += page.get_text()

    # Process the extracted text with spaCy
    doc = nlp(pdf_text)

    # Get the frequency of each word in the text
    word_frequencies = Counter(token.text for token in doc if token.is_alpha)

    # Calculate the frequency score
    total_words = len(doc)
    word_frequencies = {word: (count / total_words) for word, count in word_frequencies.items()}

    # Create a DataFrame from the word frequencies
    df = pd.DataFrame(list(word_frequencies.items()), columns=['Word', 'Frequency Score'])

    # Sort the DataFrame by Frequency Score in descending order
    df = df.sort_values(by='Frequency Score', ascending=False)

    # Calculate TF-IDF scores
    corpus = [pdf_text]  # Use the entire document as the corpus
    tfidf_vectorizer = TfidfVectorizer()
    tfidf_matrix = tfidf_vectorizer.fit_transform(corpus)
    tfidf_scores = dict(zip(tfidf_vectorizer.get_feature_names_out(), tfidf_matrix.toarray()[0]))

    # Calculate the TF-IDF relevance score for each word
    df['TF-IDF Score'] = df['Word'].map(tfidf_scores)
    
    # Normalize the TF-IDF scores
    max_tfidf_score = df['TF-IDF Score'].max()
    df['Normalized TF-IDF Score'] = df['TF-IDF Score'] / max_tfidf_score

    # Export the DataFrame to an Excel file
    excel_file_path = 'word_frequencies_tfidf.xlsx'
    df.to_excel(excel_file_path, index=False)

    # Display the DataFrame
    print(df)

    # Display the sentences, entities, and dependency parses
    print(list(doc.sents))
    print(doc.ents)
    displacy.render(next(doc.sents), style='dep', jupyter=True)

    print(f"Word frequencies exported to {excel_file_path}")


    # Sort the DataFrame by Frequency Score in descending order
    df = df.sort_values(by='Frequency Score', ascending=False)

    # Export the DataFrame to an Excel file
    excel_file_path = 'word_frequencies.xlsx'
    df.to_excel(excel_file_path, index=False)

    # Display the frequency of each word
    print(df)

    # Display the sentences, entities, and dependency parses
    print(list(doc.sents))
    print(doc.ents)
    displacy.render(next(doc.sents), style='dep', jupyter=True)

    print(f"Word frequencies exported to {excel_file_path}")
