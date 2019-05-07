debug_noisy:
	g++ DrLink_v1.sp.cpp convert_bam.cpp analyze_molecule.cpp read_chrsize.cpp detect_recombsite.cpp find_intersection.cpp sample_molecule.cpp sample_molecule_v2.cpp separate_kb_molecule.cpp separate_chr_v2.cpp check_allele_order.cpp split_string.cpp insert_raw_bp.cpp globals.cpp postprocess_bpPrediction.cpp output_co.cpp -O3 -o DrLink -lz ./gzlib/libgzstream.a
