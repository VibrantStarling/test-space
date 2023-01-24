# Tackling the veupath task and explaining my thought process (and things I'd like to do if I had time)

Ortholog mapping finds related genes across different species that have evolved from a single common ancestor. Orthology can be used to analyse differences in metabolism, mutation rates (e.g., gene duplication events, PHI scores), and phylogenetic relations. Orthologs can be identified with alignment algorithms such as DIAMOND or BLAST and are usually applied through programmes like OrthoFinder, Anvi'o or orthoMCL to assess similarities between genes. They can be used to identify single copy core gene clusters for further phylogenetic analysis.

Each ortholog is asigned a GO term based on predicted function, which are then further grouped into 3 main categories: biological, cellular, and molecular function.
Comparing GO term enrichment between genomes can provide insights into how those organisms might be interacting with their environment.

## Comparing organisms on toxoDB

I transformed orthologs for each comparison and had a quick looked look at KEGG metabolic pathways. VEG seems to cause the number of KEGG pathways found to decrease from ~2000 to 150. Without further analysis, I'd guess VEG either less complete or less annotated rather than phylogenetically distant, as VEG and ME49 share more ortholog groups with each other than with GT1 (see figure 1).

*3 way comparison*
https://toxodb.org/toxo/app/workspace/strategies/import/6d6601110c504d7b

*T. gondii ME49-VEG*
https://toxodb.org/toxo/app/workspace/strategies/import/fd466ca72324136c

*T. gondii ME49-GT1*
https://toxodb.org/toxo/app/workspace/strategies/import/27696e440b545c29

*T. gondii VEG-GT1*
https://toxodb.org/toxo/app/workspace/strategies/import/f8a1ce788685c46c

*and just to get Gt1 on it's own*
https://toxodb.org/toxo/app/workspace/strategies/import/64c52cb1d0aa5a0d


I also had quick a look on OrthoMCL just to see what it did (I've never used it before, sorry!)

*all three (includes the ME49-VEG comparison)*
https://orthomcl.org/orthomcl/app/workspace/strategies/import/f592f739fdfd3ef1

*ME49-GT1*
https://orthomcl.org/orthomcl/app/workspace/strategies/import/a625de77430d062a

*and VEG-GT1*
https://orthomcl.org/orthomcl/app/workspace/strategies/import/d1c14da847563e5d


## Finding common orthologs
Next I downloaded the orthogroup data given by toxoDB for each genome and extracted information on which orthologs are shared

```
# load modules

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn3_unweighted, venn3_circles
import seaborn as sns

# load in the gene datasets from toxodb with all ortholog identity data

GT1 = pd.read_csv("Toxoplasma_gondii_GT1.txt",delimiter="\t")
ME49 = pd.read_csv("Toxoplasma_gondii_ME49.txt", delimiter = "\t")
VEG = pd.read_csv("Toxoplasma_gondii_VEGG.txt", delimiter ="\t")

# combine all three data sets into one
df = ME49.merge(GT1,how='outer')
df = df.merge(VEG, how='outer')

# get a list of the ortholog IDs
orthologs = df['Ortholog Group'].unique()

# make a dataframe to store count data
orthopresence = pd.DataFrame({'orthologs':orthologs})
orthopresence['GT1'] = 0
orthopresence['ME49'] = 0
orthopresence['VEG'] = 0
orthopresence.set_index('orthologs',inplace=True )


df2 = df[['Organism', 'Ortholog Group']]
```
Finding the othogroups that are shared between the three genomes by counting the occurance of each ortholog by toxoplasma strain.
Sadly this takes quite a while to run and I know could probably be improved and I can actually do this faster with OrthoMCL, but here we are.

```
for i in orthologs:
    for index, row in df2.iterrows():
        if (row[1] == i) and (row[0] == 'Toxoplasma gondii VEG'):
            orthopresence['VEG'][i] += 1
        elif (row[1] == i) and (row[0] == 'Toxoplasma gondii ME49'):
            orthopresence['ME49'][i] += 1
        elif (row[1] == i) and (row[0] == 'Toxoplasma gondii GT1'):
            orthopresence['GT1'][i] += 1


# save the count data as a new file        

orthopresence.to_csv('orthopresence.csv', header=True)

orthopresence = pd.read_csv('orthopresence.csv')

# sort into categories, add totals
orthopresence['category'] = 'a'
orthopresence['total'] = orthopresence['ME49']+orthopresence['GT1']+ orthopresence['VEG']

ind = pd.np.where((orthopresence['ME49']>0) & (orthopresence['GT1']>0) & (orthopresence['VEG']>0))
ind = ind[0].tolist()
orthopresence.loc[ind,'category'] = 'ME49-VEG-GT1'

ind = pd.np.where((orthopresence['ME49']>0) & (orthopresence['GT1']==0) & (orthopresence['VEG']>0))
ind = ind[0].tolist()
orthopresence.loc[ind,'category'] = 'ME49-VEG'

ind = pd.np.where((orthopresence['ME49']>0) & (orthopresence['GT1']>0) & (orthopresence['VEG']==0))
ind = ind[0].tolist()
orthopresence.loc[ind,'category'] = 'ME49-GT1'

ind = pd.np.where((orthopresence['ME49']==0) & (orthopresence['GT1']>0) & (orthopresence['VEG']>0))
ind = ind[0].tolist()
orthopresence.loc[ind,'category'] = 'GT1-VEG'

ind = pd.np.where((orthopresence['ME49']==0) & (orthopresence['GT1']>0) & (orthopresence['VEG']==0))
ind = ind[0].tolist()
orthopresence.loc[ind,'category'] = 'GT1'

ind = pd.np.where((orthopresence['ME49']>0) & (orthopresence['GT1']==0) & (orthopresence['VEG']==0))
ind = ind[0].tolist()
orthopresence.loc[ind,'category'] = 'ME49'

ind = pd.np.where((orthopresence['ME49']==0) & (orthopresence['GT1']==0) & (orthopresence['VEG']>0))
ind = ind[0].tolist()
orthopresence.loc[ind,'category'] = 'VEG'

ind = pd.np.where((orthopresence['orthologs']=='N/A (orthology not determined for transcripts that are not protein-coding)'))
ind = ind[0].tolist()
orthopresence.loc[ind,'category'] = 'no-ortholog'

nogroup = orthopresence.loc[orthopresence['category']=='no-ortholog']
groups = orthopresence.groupby('category')['orthologs'].count().drop(index='no-ortholog')

```

## Visualisation
I've opted for a venn diagram here for the sake of time.
Previously, I have used complexupset in R for something like this to make an upsetplot.

Ideally, if I had time,  I'd also like to do an enrichment analysis between each genome to see what functions are common/different between each set.
I'd take the GO enrichment analysis from the toxodb, but with visulaisation like g:GOSt displaying fold enrichment and p-values https://biit.cs.ut.ee/gprofiler/page/docs

```
v = venn3_unweighted(subsets=(groups['GT1'], groups['ME49'], groups['ME49-GT1'], groups['VEG'], groups['GT1-VEG'], groups['ME49-VEG'], groups['ME49-VEG-GT1']), 
          set_labels = ('GT1', 'ME49', 'VEG'))
plt.savefig('venn_diagram.png', format='png')
plt.show()

```
![Venn diagram](https://github.com/VibrantStarling/test-space/blob/main/venn_diagram.png)

Figure 1. A venn diagram showing the number of orthologs common between each of the three Toxoplasma gondii strains

Table 1. Non coding genes with no orthology: 

| Genome        | Genes with no orthology
| ------------- |-------------| 
| GT1      | 177 |
| ME49      | 598      |  
| VEG | 153      |   
               


