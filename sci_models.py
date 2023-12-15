import os
import spacy
from collections import Counter
import math
import pandas as pd

# Folder path containing txt files with biomedical abstracts
folder_path = '/home/konstantinos/abstract_txt'

# Output folder for saving processed abstracts and entities
output_folder = '/home/konstantinos/sci_outputs'
abstract_entities_folder = os.path.join(output_folder, 'abstract_entities_txt')

# Check if the 'abstract_entities_txt' folder exists, if not, create it
if not os.path.exists(abstract_entities_folder):
    os.makedirs(abstract_entities_folder)

# Load the SciSpacy models
nlp_bionlp13cg = spacy.load("en_ner_bionlp13cg_md")
nlp_bc5cdr = spacy.load("en_ner_bc5cdr_md")


def calculate_pmi(co_occurrence, entity_counts, total_abstracts):
    pmi_values = {}
    for entities, count in co_occurrence.items():
        entity1, entity2 = entities
        count_entity1 = entity_counts[entity1]
        count_entity2 = entity_counts[entity2]

        # Check for zero denominator to avoid division by zero
        if count_entity1 == 0 or count_entity2 == 0:
            pmi = 0.0
        else:
            pmi = math.log2((count * total_abstracts) / (count_entity1 * count_entity2))
        pmi_values[(entity1, entity2)] = pmi
    return pmi_values


def sigmoid(x):
    return 2 / (1 + math.exp(-x)) - 1


def calculate_normalized_pmi(co_occurrence, entity_counts, total_abstracts):
    pmi_values = {}

    for (entity1, entity2), count in co_occurrence.items():
        count_entity1 = entity_counts[entity1]
        count_entity2 = entity_counts[entity2]

        # Check for zero denominator to avoid division by zero
        if count_entity1 == 0 or count_entity2 == 0:
            normalized_pmi = 0.0
        else:
            p_entity1 = count_entity1 / total_abstracts
            p_entity2 = count_entity2 / total_abstracts

            # Calculate raw PMI
            raw_pmi = math.log2((count * total_abstracts) / (count_entity1 * count_entity2))

            # Calculate PMI Ratio for normalization
            pmi_ratio = raw_pmi / (-math.log2(p_entity1 * p_entity2))

            # Apply sigmoid function for scaling to the range [-1, 1]
            normalized_pmi = sigmoid(pmi_ratio)

        pmi_values[(entity1, entity2)] = normalized_pmi

    return pmi_values


# Initialize counters
entity_counts_bionlp13cg = Counter()
entity_counts_bc5cdr = Counter()
co_occurrence_bionlp13cg = Counter()
co_occurrence_bc5cdr = Counter()
total_abstracts = 0

# Iterate through each text file in the input folder
for file_name in os.listdir(folder_path):
    file_path = os.path.join(folder_path, file_name)

    # Read the content of the abstract from the text file
    with open(file_path, 'r', encoding='utf-8') as f:
        abstract_text = f.read()

    # Process the abstract with both models
    doc_bionlp13cg = nlp_bionlp13cg(abstract_text)
    doc_bc5cdr = nlp_bc5cdr(abstract_text)

    # Update entity counts
    entity_counts_bionlp13cg.update([(ent.text, ent.label_) for ent in doc_bionlp13cg.ents])
    entity_counts_bc5cdr.update([(ent.text, ent.label_) for ent in doc_bc5cdr.ents])

    # Update co-occurrence counts
    entities_bionlp13cg = [(ent.text, ent.label_) for ent in doc_bionlp13cg.ents]
    entities_bc5cdr = [(ent.text, ent.label_) for ent in doc_bc5cdr.ents]

    for entity1 in entities_bionlp13cg:
        for entity2 in entities_bc5cdr:
            co_occurrence_bionlp13cg[(entity1, entity2)] += 1

    for entity1 in entities_bc5cdr:
        for entity2 in entities_bionlp13cg:
            co_occurrence_bc5cdr[(entity1, entity2)] += 1

    total_abstracts += 1

    # Create a text file for each abstract with entities and their types
    abstract_entities_file_path = os.path.join(abstract_entities_folder, f"{file_name}_entities.txt")
    with open(abstract_entities_file_path, 'w', encoding='utf-8') as abstract_entities_file:
        abstract_entities_file.write(f"Abstract:\n{abstract_text}\n\n")
        abstract_entities_file.write("Entities from en_ner_bionlp13cg_md:\n")
        for entity, entity_type in entities_bionlp13cg:
            abstract_entities_file.write(f"{entity} ({entity_type})\n")
        abstract_entities_file.write("\nEntities from en_ner_bc5cdr_md:\n")
        for entity, entity_type in entities_bc5cdr:
            abstract_entities_file.write(f"{entity} ({entity_type})\n")

# Calculate PMI
pmi_values_bionlp13cg = calculate_pmi(co_occurrence_bionlp13cg, entity_counts_bionlp13cg, total_abstracts)
pmi_values_bc5cdr = calculate_pmi(co_occurrence_bc5cdr, entity_counts_bc5cdr, total_abstracts)

# Calculate normalized PMI
pmi_normalized_values_bionlp13cg = calculate_normalized_pmi(co_occurrence_bionlp13cg, entity_counts_bionlp13cg,
                                                            total_abstracts)
pmi_normalized_values_bc5cdr = calculate_normalized_pmi(co_occurrence_bc5cdr, entity_counts_bc5cdr, total_abstracts)

# Convert the results to DataFrames
df_pmi_bionlp13cg = pd.DataFrame(list(pmi_values_bionlp13cg.items()), columns=['Entities', 'PMI'])
df_pmi_bc5cdr = pd.DataFrame(list(pmi_values_bc5cdr.items()), columns=['Entities', 'PMI'])
df_pmi_normalized_bionlp13cg = pd.DataFrame(list(pmi_normalized_values_bionlp13cg.items()), columns=['Entities', 'Normalized PMI'])
df_pmi_normalized_bc5cdr = pd.DataFrame(list(pmi_normalized_values_bc5cdr.items()), columns=['Entities', 'Normalized PMI'])

# Save the DataFrames to Excel files
excel_file_path_bionlp13cg = os.path.join(output_folder, "pmi_results_bionlp13cg.xlsx")
excel_file_path_bc5cdr = os.path.join(output_folder, "pmi_results_bc5cdr.xlsx")
excel_file_path_bionlp13cg_normalized = os.path.join(output_folder, "pmi_results_bionlp13cg_normalized.xlsx")
excel_file_path_bc5cdr_normalized = os.path.join(output_folder, "pmi_results_bc5cdr_normalized.xlsx")

df_pmi_bionlp13cg.to_excel(excel_file_path_bionlp13cg, index=False)
df_pmi_bc5cdr.to_excel(excel_file_path_bc5cdr, index=False)
df_pmi_normalized_bionlp13cg.to_excel(excel_file_path_bionlp13cg_normalized, index=False)
df_pmi_normalized_bc5cdr.to_excel(excel_file_path_bc5cdr_normalized, index=False)

print(f"PMI results saved to {excel_file_path_bionlp13cg}, {excel_file_path_bc5cdr}, {excel_file_path_bionlp13cg_normalized}, and {excel_file_path_bc5cdr_normalized}")
print(f"Abstracts and entities saved to {abstract_entities_folder}")
