import scispacy
import spacy
import os

# Folder path containing txt files with biomedical abstracts
input_folder = '/home/konstantinos/abstract_txt'

# Output folder for saving processed abstracts
output_folder = '/home/konstantinos/sci_outputs'

# Load the SciSpacy models
nlp_bionlp13cg = spacy.load("en_ner_bionlp13cg_md")
nlp_bc5cdr = spacy.load("en_ner_bc5cdr_md")

# Iterate through each text file in the input folder
for file_name in os.listdir(input_folder):
    file_path = os.path.join(input_folder, file_name)

    # Read the content of the abstract from the text file
    with open(file_path, 'r', encoding='utf-8') as f:
        abstract_text = f.read()

    # Process the abstract with both models
    doc_bionlp13cg = nlp_bionlp13cg(abstract_text)
    doc_bc5cdr = nlp_bc5cdr(abstract_text)

    # Save the results to a new text file in the output folder
    output_file_path = os.path.join(output_folder, f"{file_name}_entities.txt")
    with open(output_file_path, 'w', encoding='utf-8') as output_file:
        output_file.write("Original Abstract:\n")
        output_file.write(f"{abstract_text}\n\n")

        output_file.write("Entities from en_ner_bionlp13cg_md:\n")
        for ent in doc_bionlp13cg.ents:
            output_file.write(f"{ent.text} - {ent.label_}\n")

        output_file.write("\nEntities from en_ner_bc5cdr_md:\n")
        for ent in doc_bc5cdr.ents:
            output_file.write(f"{ent.text} - {ent.label_}\n")

    print(f"Entities extracted from {file_name} and saved to {output_file_path}")
