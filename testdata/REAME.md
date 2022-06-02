## Kmer subtraction

Jellyfish count dump files can be used to identify (subtract) kmers that exist in one file but not the other. Jellyfish can accomplish this task but doesn't take counts into account. Therefore, Arthur Georges has created a script that can first filter kmer count dump files for specified threshold and identify unique kmers in each file. 
