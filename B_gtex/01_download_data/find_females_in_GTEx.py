# In this script, we want to return the GTEx ID that are females

import pandas as pd

data = pd.read_csv("c://Users/tuyen/Documents/postdoc_asu/projects/GTEx/GTEx_v8/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt", sep="\t")

# print (data.head())

# Subset for Females
females = data[data['SEX'] == 2]
print (females.head())

# Save to a file
females.to_csv("c://Users/tuyen/Documents/postdoc_asu/projects/GTEx/GTEx_v8/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS_females.txt", sep="\t", index=False)