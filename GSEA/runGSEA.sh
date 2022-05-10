for tissue in B H SM; 
	do
	for geneset in c5.bp c5.cc c5.mf c2.cp.biocarta c2.cp.kegg c2.cp.reactome
	do 
		java -Xmx1024m -cp gsea-3.0.jar xtools.gsea.Gsea \
		-res CPM_noduplicates.txt \
		-cls phenotypelabels.cls#$tissue"H_versus_"$tissue"PH" \
		-gmx msigdb_v6.1_files_to_download_locally/msigdb_v6.1_GMTs/$geneset.v6.1.symbols.gmt \
		-collapse false \
		-mode Max_probe \
		-norm meandiv \
		-nperm 1000 \
		-permute gene_set \
		-rnd_type no_balance \
		-scoring_scheme weighted \
		-rpt_label $tissue.$geneset \
		-metric Signal2Noise \
		-sort real \
		-order descending \
		-create_gcts false \
		-create_svgs false \
		-include_only_symbols true \
		-make_sets true \
		-median false \
		-num 100 \
		-plot_top_x 20 \
		-rnd_seed 1229 \
		-save_rnd_lists false \
		-set_max 500 \
		-set_min 15 \
		-zip_report false \
		-out . \
		-gui false
	done
done
