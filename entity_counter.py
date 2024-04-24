import nltk
import spacy
import os
from collections import Counter

nltk.download('punkt')

def separate_into_sentences(text):
    """
    Separates a text into individual sentences using NLTK's sentence tokenizer.
    """
    sentences = nltk.sent_tokenize(text)
    return sentences

def extract_disease_entities(sentences, nlp_model):
    """
    Extracts disease entities using SpaCy's named entity recognition.
    """
    doc = nlp_model(' '.join(sentences))
    disease_entities = [X.text for X in doc.ents if X.label_ == 'DISEASE']
    return disease_entities

def extract_gene_entities(sentences, nlp_model):
    """
    Extracts gene entities using SpaCy's named entity recognition.
    """
    doc = nlp_model(' '.join(sentences))
    gene_entities = [X.text for X in doc.ents if X.label_ == 'GENE_OR_GENE_PRODUCT']
    return gene_entities

# Load SpaCy models
nlp_bc5cdr = spacy.load('en_ner_bc5cdr_md')
nlp_bionlp = spacy.load('en_ner_bionlp13cg_md')

# Define input and output folder paths
input_folder = 'abstracts_txt'
output_folder = 'entity_results'

# Create output folder if it doesn't exist
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# Loop through the input folder and its subfolders
for root, dirs, files in os.walk(input_folder):
    for filename in files:
        if filename.endswith('.txt'):
            # Get the relative path of the input file
            relative_path = os.path.relpath(root, input_folder)
            
            # Create corresponding output folder
            output_subfolder = os.path.join(output_folder, relative_path)
            if not os.path.exists(output_subfolder):
                os.makedirs(output_subfolder)
            
            file_path = os.path.join(root, filename)
            with open(file_path, 'r', encoding='utf-8') as f:
                abstract_text = f.read()
            sentences = separate_into_sentences(abstract_text)
            gene_entities = extract_gene_entities(sentences, nlp_bionlp)
            disease_entities = extract_disease_entities(sentences, nlp_bc5cdr)
            
    
            
            # Write the gene and disease frequencies to a separate file
            output_file_freq = os.path.join(output_subfolder, f"{os.path.splitext(filename)[0]}_entity_frequencies.txt")
            with open(output_file_freq, 'w') as out_freq_f:
                gene_freq = Counter(gene_entities)
                disease_freq = Counter(disease_entities)
                
                out_freq_f.write("Gene/Gene Product Frequencies:\n")
                for gene, freq in gene_freq.items():
                    out_freq_f.write(f"{gene}: {freq}\n")
                
                out_freq_f.write("\nDisease Frequencies:\n")
                for disease, freq in disease_freq.items():
                    out_freq_f.write(f"{disease}: {freq}\n")

                print(f"Frequencies of genes and diseases saved to {output_file_freq}")
