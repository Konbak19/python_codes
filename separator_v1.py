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

text = "Germ line mutations of the BRCA1 gene confer a high risk of breast cancer and ovarian cancer to female mutation carriers. The BRCA1 protein is involved in the regulation of DNA repair. How specific tumor-associated mutations affect the molecular function of BRCA1, however, awaits further elucidation. Cell lines that harbor BRCA1 gene mutations are invaluable tools for such functional studies. Up to now, the HCC1937 cell line was the only human breast cancer cell line with an identified BRCA1 mutation. In this study, we identified three other BRCA1 mutants from among 41 human breast cancer cell lines by sequencing of the complete coding sequence of BRCA1. Cell line MDA-MB-436 had the 5396 + 1G>A mutation in the splice donor site of exon 20. Cell line SUM149PT carried the 2288delT mutation and SUM1315MO2 carried the 185delAG mutation. All three mutations were accompanied by loss of the other BRCA1 allele. The 185delAG and 5396 + 1G>A mutations are both classified as pathogenic mutations. In contrast with wild-type cell lines, none of the BRCA1 mutants expressed nuclear BRCA1 proteins as detected with Ab-1 and Ab-2 anti-BRCA1 monoclonal antibodies. These three new human BRCA1 mutant cell lines thus seem to be representative breast cancer models that could aid in further unraveling of the function of BRCA1. "
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
