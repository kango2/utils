- Downloaded gene set files from
`http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/6.1/msigdb_v6.1_files_to_download_locally.zip on 31st May 2018.`

- Downloaded GSEA program from `http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/software/gsea-3.0.jar on 31st May 2018.`

Selected following genesets:
* c5.bp = Biological Processes GO terms
* c5.cc = Cellular Component GO terms
* c5.mf = Molecular Functions GO terms
* c2.cp.biocarta = Biocarta Genesets
* c2.cp.kegg = Kegg Pathway Genesets
* c2.cp.reactome = Reactome Pathway Genesets

B, H and SM stands for different tissue types.

Created Phenotype labels file `phenotypelabels.cls`

Place following files/directories in one directory.

1. msigdb_v6.1_files_to_download_locally
2. gsea-3.0.jar
3. Gene expression table (CPM_noduplicates.txt)
4. runGSEA.sh
5. phenotypelabels.cls

Type following at the command line:

`sh runGSEA.sh`

**NOTES**
- Modify `runGSEA.sh` for adding or removing gene sets.
- Modify `phenotypelabels.cls` for modifying phenotype labels and make subsequent change in the `runGSEA.sh` file as required to accomodate for the change in phenotype labels.

