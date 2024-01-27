import pyhgvs as hgvs
import math
import copy
from operator import itemgetter
import blosum as bl

bl_mat = bl.BLOSUM(62)
bl_dict = dict(bl_mat)

"""
genetic_code = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'
}
"""

amino_acids = {'Cys': 'C', 'Asp': 'D', 'Ser': 'S', 'Gln': 'Q', 'Lys': 'K',
     'Ile': 'I', 'Pro': 'P', 'Thr': 'T', 'Phe': 'F', 'Asn': 'N', 
     'Gly': 'G', 'His': 'H', 'Leu': 'L', 'Arg': 'R', 'Trp': 'W', 
     'Ala': 'A', 'Val':'V', 'Glu': 'E', 'Tyr': 'Y', 'Met': 'M',
     'Asx': 'B', 'Glx': 'Z', 'Xaa': "X", "Xle": "J"}

chemical_properties = {"Aliphatic": ["G", "A", "V", "L", "I", "P", "M"],
                        "Aromatic": ["F", "Y", "W"],
                        "Polar uncharged": ["S", "T", "C", "N", "Q"],
                        "Acidic": ["D", "E"],
                        "Basic": ["R", "H", "K"],
                        "Unknown": ["B", "Z", "X", "J"]}

def fasta_reader(filename): #Taken from the labs
    seqDict = {}

    fileIn = open(filename, "r")

    sequence = ""
    header = ""

    for line in fileIn:
        if line[0] == ">":
            if sequence:
                seqDict[header] = sequence
            header = line[1:].strip()
            sequence = ""
        else:
            sequence += line.strip()
    seqDict[header] = sequence

    fileIn.close()

    return seqDict

def shannon_entropy(msa_dict): #Adapted from the function with the same name from Lab 4 of BIO 310 (Spring 2022)
    msa_length = len(list(msa_dict.values())[0])
    number_of_proteins = len(list(msa_dict.keys()))

    entropy_list = [] 
    aa_count_dict = {}
    for i in range(msa_length):
        aa_count_dict[i] = {}
        for seq in msa_dict.values():
            aa = seq[i]
            if not aa in aa_count_dict[i].keys():
                aa_count_dict[i][aa] = 1
            else:
                aa_count_dict[i][aa] += 1
    for pos in aa_count_dict.keys():
        for aa,count in aa_count_dict[pos].items():
            P_i = count/float(number_of_proteins)                                   
            entropy_i = P_i*(math.log(P_i,2))
            entropy_list.append(entropy_i)
    sh_entropy = -(sum(entropy_list))
    
    return sh_entropy

"""
def SPcalc(msa_dict, substitution_matrix): #Adapted from the function with the same name from Lab 4 of BIO 310 (Spring 2022)
    SP = 0

    msa_length = len(list(msa_dict.values())[0])
    number_of_proteins = len(list(msa_dict.keys()))
    protein_sequences = list(msa_dict.values())

    for i in range(msa_length):
        for j in range(number_of_proteins):
            for k in range(j+1,number_of_proteins):
                if protein_sequences[j][i] == "-" or protein_sequences[k][i] == "-":
                    SP += 0
                else:
                    SP += substitution_matrix[protein_sequences[j][i] + protein_sequences[k][i]]
    return SP
"""

def list_to_fasta(seq_list, file_name): #Writes a list of sequences into a FASTA file
    fasta_file = open(file_name + ".fa", "a")
    
    for i in range(0, len(seq_list)):
        fasta_file.write(">" + str(i) + "\n")
        fasta_file.write(seq_list[i] + "\n")

    fasta_file.close()

def mutation_generator(ref_seq, mut_list, type = "aa", ref_transcript_name = "NM_000492.3"): #Generates and returns the mutated sequences described in the database; assumes point substitutions
    #Note that nucleotide numbering in HGVS format starts from 1.

    mutated_seqs = []

    for mutation in mut_list:
        ref_copy = ref_seq
        temp_copy = [*ref_copy]

        hgvs_info = hgvs.HGVSName(ref_transcript_name + ":" + mutation)

        if type == "aa":
            pos = hgvs_info.start
            new_aa = amino_acids[hgvs_info.alt_allele]

            temp_copy[pos - 1] = new_aa
        elif type == "nuc":
            pos = hgvs_info.cdna_start.coord
            new_nuc = hgvs_info.alt_allele

            temp_copy[pos - 1] = new_nuc
        else:
            pass
        
        new_seq = "".join(temp_copy)
        mutated_seqs.append(new_seq)

    return mutated_seqs

def variant_interpreter(mutation, type = "aa", ref_transcript_name = "NM_000492.3"): #Interprets a variant in HGVS notation
    pos = 0
    subst = ""
    ref = ""

    hgvs_info = hgvs.HGVSName(ref_transcript_name + ":" + mutation)
    
    if type == "aa":
        pos = hgvs_info.start - 1
        ref = amino_acids[hgvs_info.ref_allele]
        subst = amino_acids[hgvs_info.alt_allele]
    elif type == "nuc":
        pos = hgvs_info.cdna_start.coord - 1
        ref = hgvs_info.ref_allele
        subst = hgvs_info.alt_allele
    else:
        pass

    score = bl_dict[amino_acids[hgvs_info.ref_allele] + subst]
    
    return pos, ref, subst, score

"""
def translator(DNA): #Translates a DNA sequence to a protein sequence; actually redundant
    protein_seq = ""

    i = 0
    while (i <= len(DNA) - 3):
        codon = DNA[i:i+3]
        protein_seq += genetic_code[codon]

        i += 3
    
    return protein_seq
"""

def gap_pass(sequence, loc): #Returns the gapped position of a location in a gapless sequence
    pos = 0
    counter = 0

    for x in sequence:
        if pos < loc:
            if x == "-":
                counter += 1
                continue
            else:
                counter += 1
                pos += 1
    
    return counter

def profiler(seq_list, type = "aa", result = "percentage"): #Profiles a list of sequences; it assumes that each sequence has the same length and dismisses gaps
    profile = []
    chrList = {}

    if type == "aa":
        chrList = {
            "A": 0, "C": 0, "D": 0,
            "E": 0, "F": 0, "G": 0,
            "H": 0, "I": 0, "K": 0,
            "L": 0, "M": 0, "N": 0,
            "P": 0, "Q": 0, "R": 0,
            "S": 0, "T": 0, "V": 0,
            "W": 0, "Y": 0, "B": 0,
            "Z": 0, "X":0, "J": 0
        }
    elif type == "nuc":
        chrList = {
            "A": 0,
            "T": 0,
            "G": 0,
            "C": 0
        }
    else:
        pass

    for i in range(0, len(seq_list[0])):
        profile.append(chrList.copy())

    for j in range(0, len(seq_list)):
        for k in range(0, len(seq_list[j])):
            if seq_list[j][k] in profile[k]: #If it is not a gap
                profile[k][seq_list[j][k]] += 1
            else: #If it is a gap
                pass

    if result == "percentage":
        for m in range(0, len(profile)):
            sum = 0

            for key in profile[m]:
                sum += profile[m][key]
            
            for key in profile[m]:
                profile[m][key] = round(profile[m][key] / sum, 2)
    else:
        pass

    return profile

def profile_modifier(profile, ref_seq, seq_dict, ref_seq_key = "NP 000483.3 cystic fibrosis transmembrane conductance regulator Homo sapiens"): #Returns a profile that retains only the positions that are not gapped in the sequence of interest
    modified_profile = []

    for i in range(0, len(ref_seq)):
        loc = gap_pass(seq_dict[ref_seq_key], i)
        modified_profile.append(profile[loc])

    return modified_profile

def vus_interpreter(vus_list, profile, msa_dict, regions, variant_key = "NP 000483.3 cystic fibrosis transmembrane conductance regulator Homo sapiens"): #Provides guiding information on whether a VUS is pathogenic; assumes that the VUSes describe amino acid sequences
    for vus in vus_list: #Inspects each VUS and prints the most common character vs. the substitution
        prediction = ""

        pos, ref, subst, bl_score = variant_interpreter(vus)

        max_freq_dict = dict(sorted(profile[pos].items(), key = itemgetter(1), reverse = True)[:4])
        max_freq_keys = sorted(max_freq_dict, key = max_freq_dict.get, reverse = True)
        max_freq_values = sorted(max_freq_dict.values(), reverse = True)
        max_freq = max(max_freq_values)
        
        max_chrs = []

        for key in profile[pos]:
            if profile[pos][key] == max_freq:
                max_chrs.append(key)
            else:
                max_chrs.append(key)

        entropy_before = shannon_entropy(msa_dict)

        modified_msa_dict = copy.deepcopy(msa_dict)

        replaced_seq = msa_dict[variant_key]
        replaced_seq_list = [*replaced_seq]
        actual_pos = gap_pass(replaced_seq, pos)
        replaced_seq_list[actual_pos] = subst
        replaced_seq = "".join(replaced_seq_list)

        modified_msa_dict[variant_key] = replaced_seq

        entropy_after = shannon_entropy(modified_msa_dict)
        entropy_change = entropy_after - entropy_before

        max_freq_chem_groups = []
        chemical_group_after = ""
        ref_chemical_group = ""

        for aa in max_freq_keys:
            chemical_group_before = ""
            for key in chemical_properties:
                if aa in chemical_properties[key]:
                    chemical_group_before = key
                else:
                    pass

            max_freq_chem_groups.append(chemical_group_before)
        
        for key in chemical_properties:
            if (ref in chemical_properties[key]) and (subst in chemical_properties[key]):
                chemical_group_after = key
                ref_chemical_group = key
            elif ref in chemical_properties[key]:
                ref_chemical_group = key
            elif subst in chemical_properties[key]:
                chemical_group_after = key
            else:
                pass
        
        inCriticalPart = False
        for key in regions:
            for region in regions[key]:
                start_pos = gap_pass(msa_dict[variant_key], region[0])
                end_pos = gap_pass(msa_dict[variant_key], region[1])

                if start_pos == end_pos:
                    if pos == start_pos:
                        inCriticalPart = True
                        break
                    else:
                        pass
                elif start_pos <= pos <= end_pos:
                    inCriticalPart = True
                    break
                else:
                    pass

        if inCriticalPart == True:
            if ref == subst:
                prediction = "Non-pathogenic"
            elif bl_score >= 0 and ref_chemical_group == chemical_group_after:
                prediction = "Likely non-pathogenic"
            else:
                prediction = "Likely pathogenic"
        else:
            if ref == subst:
                prediction = "Non-pathogenic"
            elif bl_score >= 0 and ref_chemical_group == chemical_group_after:
                prediction = "Likely non-pathogenic"
            elif bl_score >= 0 and subst in max_freq_keys[0:2]:
                prediction = "Likely non-pathogenic"
            elif bl_score >= 0 and chemical_group_after in max_freq_chem_groups[0:2]:
                prediction = "Likely non-pathogenic"
            elif bl_score >= 0:
                prediction = "Possibly (<likely) non-pathogenic"
            else:
                prediction = "Likely pathogenic"

        print("The most common characters in the profile: ")
        for i in range(0, len(max_freq_keys)):
            print("-" + max_freq_keys[i] + " with " + str(max_freq_values[i] * 100) + "%")
        print("Chemical properties of the most common characters: ")
        for j in range(0, len(max_freq_keys)):
            print("-" + max_freq_chem_groups[j] + " with " + str(max_freq_values[j] * 100) + "%")
        print("Substitution: " + ref + " -> " + subst)
        print("Position in the aligned sequence: " + str(actual_pos))
        print("Chemical property of the reference character: " + ref_chemical_group)
        print("Chemical property of the substitution character: " + chemical_group_after)
        print("Change in entropy: " + str(entropy_change))
        print("Substitution score (BLOSUM62): " + str(bl_score))
        print("In a critical part of the protein? " + str(inCriticalPart))
        print("Automated prediction: " + prediction + "\n")