grep "^>" 1000_aligned_ref_sequences.fa | sort | uniq -c | sort -nr | head #Checks whether there are duplicate sequence names in 1000_aligned_ref_sequences.fa
