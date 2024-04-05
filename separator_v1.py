import nltk
import scispacy
import spacy

nltk.download('punkt')

def separate_into_sentences(text):
    """
    Separates a text into individual sentences using NLTK's sentence tokenizer.
    """
    sentences = nltk.sent_tokenize(text)
    return sentences

text = "Currently, despite the use of a preventive vaccine for several decades as well as the use of effective and well-tolerated viral suppressive medications since 1998, approximately 250 million people remain infected with the virus that causes hepatitis B worldwide. Hepatitis C virus (HCV) and hepatitis B virus (HBV) are the leading causes of liver cancer and overall mortality globally, surpassing malaria and tuberculosis. Linkage to care is estimated to be very poor both in developing countries and in high-income countries, such as the United States, countries in Western Europe, and Japan. In the United States, by CDC estimates, only one-third of HBV-infected patients or less are aware of their infection.BRCA1. "
sentences = separate_into_sentences(text)


nlp_bc5cdr = spacy.load('en_ner_bc5cdr_md')
nlp_bionlp = spacy.load('en_ner_bionlp13cg_md')

def extract_disease_entities(sentences):

    """

    Extracts disease entities using SpaCy's named entity recognition.


    """

    doc = nlp_bc5cdr(' '.join(sentences))

    disease_entities = [(X.text, X.label) for X in doc.ents if X.label_ == 'DISEASE']

    return disease_entities


def extract_gene_entities(sentences):



    doc_bionlp = nlp_bionlp(' '.join(sentences))

    gene_entities = [(X.text, X.label) for X in doc_bionlp.ents if X.label_ == 'GENE_OR_GENE_PRODUCT']

    return gene_entities

gene_entities = extract_gene_entities(sentences)
disease_entities = extract_disease_entities(sentences)
print('Gene/Gene Product Entities:')

for ent in gene_entities:

    print(ent)

print('\nDisease Entities:')

for ent in disease_entities:

    print(ent)
