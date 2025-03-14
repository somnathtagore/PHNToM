# PHaNToM
PHaNToM: Systematic Pan-cancer assessment and validation of Hyper, Hypo and Neo-morph mutations

Cancer subtypes present extreme heterogeneity in terms of their molecular alterations’ repertoire, both across distinct samples and within single cells in the same sample, despite their relatively conserved transcriptional states. This has emerged as a major determinant of resistance to standard of care therapies including targeted therapies. The effect of individual mutations can be characterized as either hyper-morphic (Gain of Functions, GOF) or hypo-morphic (Loss of Functions, LOF), depending on whether they induce an increase or decrease of the encoded protein’s enzymatic or regulatory activity. However, neo-morphic (NEO) mutations have also been reported, which are responsible for producing entirely novel functional effects by either abrogating physiologic protein-protein or transcriptional interactions or by introducing non-physiologic ones. 

We have developed a novel computational method called PHaNToM (Protein-activity based identification of Hyper-, Hypo-, and Neo-morphic effecTors of Mutations) for characterizing individual mutations as either hypo, hyper, or neo-morphic by leveraging the activity of their downstream effectors—including transcription factors (TFs) and co-factors (co-TFs)—as a highly-reproducible gene reporter assay. TF and co-TF activity can be accurately measured by the VIPER algorithm (Virtual Inference of Protein-activity by Enriched Regulon analysis), based on the differential expression of their regulatory targets, including at the single cell level.

# Tutorial

In the tutorial section, there are three notebooks, 1) 1.consensus_signature.ipynb (for generating the consensus GOF or LOF signature), 2) 2.aREA.mutations.ipynb (for performing analytical rank-based enrichment analysis), 3) 3.mutation_classfication.ipynb (for classifying VUFS mutations into GOF, LOF, NEU or NEO). Associated files are are also provided.

# Scripts

In this section, associated scripts are provided, if notebooks are not preferred.
