import tools

ref_cftr_seq = (tools.fasta_reader("cftr_seq.fa"))["NM_000492.3 Homo sapiens CF transmembrane conductance regulator (CFTR), CDS"] #The coding sequence of the CFTR transcript referenced by LOVD

cftr_regions = {"Domains": [[81, 365], [423, 646], [859, 1155], [1210, 1443]], 
                "Binding sites": [[401, 401], [434, 434], [458, 465], [493, 493], [1219, 1219], [1244, 1251]], 
                "Motifs": [[1478, 1480]], 
                "Regions": [[654, 831], [1386, 1480], [1452, 1480]]
                }

vus_list = ["p.Asp1445Asn", #Because LOVD offers no convenient way of extracting variome data, 60 random sequences with amino acid variants were manually chosen and placed in a list
                "p.Glu1433Lys",
                "p.Met1407Thr",
                "p.Thr1396Pro",
                "p.Asp1377His",
                "p.Ile1366Thr",
                "p.Ala1364Val",
                "p.Arg1358Ser",
                "p.Gln1352His",
                "p.Gln1352His",
                "p.Phe1337Val",
                "p.Ile1328Thr",
                "p.Asp1312Tyr",
                "p.Asn1303Lys",
                "p.Ala1285Val",
                "p.Asp1270Asn",
                "p.Thr1246Ile",
                "p.Ser1235Arg",
                "p.Asn1229=",
                "p.Thr1220Ile",
                "p.Thr1220Ile",
                "p.Ser1188Leu",
                "p.Asp1152His",
                "p.Ile1139Val",
                "p.Tyr1092His",
                "p.Arg1070Gln",
                "p.Arg1070Trp",
                "p.Thr1057Arg",
                "p.Tyr1014Cys",
                "p.Ser912Leu",
                "p.Val920Met",
                "p.Tyr919Cys",
                "p.Ser912Leu",
                "p.Glu826Lys",
                "p.Ile807Met",
                "p.Glu585Lys",
                "p.Ile601Phe",
                "p.Glu585Lys",
                "p.Thr582Ser",
                "p.Val754Met",
                "p.Ala349Val",
                "p.Arg334Gln",
                "p.Arg297Gln",
                "p.Arg297Trp",
                "p.Lys294Ile",
                "p.Met156Val",
                "p.Ala120Thr",
                "p.Arg117His",
                "p.Leu88Phe",
                "p.Arg75Gln",
                "p.Arg74Gln",
                "p.Glu56Lys",
                "p.Ser42Phe",
                "p.Arg31Leu",
                "p.Arg31His",
                "p.Arg31Cys",
                "p.Val11Ile",
                "p.Asn287Tyr",
                "p.Ala280Thr",
                "p.Glu217Gly",
                "p.Gly241Arg",
                "p.Tyr301Cys",
                "p.Leu333Phe",
                "p.Leu320Val",
                "p.Arg334Trp",
                "p.Thr451Ser",
                "p.Asp443Tyr",
                "p.Pro477Thr",
                "p.Ser485Cys",
                "p.Ser489Leu",
                "p.Ile601Phe",
                "p.Asp614Gly"]

mutated_seqs = tools.mutation_generator(ref_cftr_seq, vus_list) #Full CFTR sequences of the patients (each differ from the reference sequence by one nucleotide)
tools.list_to_fasta(mutated_seqs, "mutated_sequences")

#Analysis based on 1000 sequences obtained by BLAST and aligned with MUSCLE
blast_dict = tools.fasta_reader("1000_aligned_ref_sequences_MUSCLE.fa")
blast_seqs = list(blast_dict.values())

healthy_profile_1000 = tools.profiler(blast_seqs)
modified_healthy_profile_1000 = tools.profile_modifier(healthy_profile_1000, ref_cftr_seq, blast_dict)

tools.vus_interpreter(vus_list, modified_healthy_profile_1000, blast_dict, cftr_regions)

#Analysis based on 800 sequences obtained by PSIBLAST and aligned with MAFFT
psiblast_dict = tools.fasta_reader("800_aligned_ref_sequences_MAFFT.fa")
psiblast_seqs = list(psiblast_dict.values())

healthy_profile_800 = tools.profiler(psiblast_seqs)
modified_healthy_profile_800 = tools.profile_modifier(healthy_profile_800, ref_cftr_seq, psiblast_dict, "NP_000483.3 cystic fibrosis transmembrane conductance regulator [Homo sapiens]")

tools.vus_interpreter(vus_list, modified_healthy_profile_800, psiblast_dict, cftr_regions, "NP_000483.3 cystic fibrosis transmembrane conductance regulator [Homo sapiens]")