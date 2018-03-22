
FindingInfo - Python package for computing linguistic complexity and Shannon's entropy in a genome sequence
-----------------------------------------------------------------------------------------------------------

**Author: Aditi Gupta; Rutgers University, Newark NJ.**

![FindingInfo Logo](/findingInfo_logo.png)


PC: Todd Richman


**Usage:**
```
$ python findinginfo.py -f <genome sequence in fasta format> -w <window size>
```

**Questions/How to cite:**

Please email me at ag1349@njms.rutgers.edu or aditi9783@gmail.com if you have any questions.
If you use FindingInfo, please contact me for information on how to cite this work (publication coming soon).

**Computing linguistic complexity**

Linguistic complexity measures the extent to which a sequence contains the non-repetitive combinations of letters from the alphabet [1]. For a sequence of length n, its complexity score is defined as follows:

![LC equation](/LinguisticComplex_eqn.png)<!-- .element height="10%" width="10%" -->

where Ui is the ratio of the actual number to the maximum possible number of all combinations of letters in a subsequence of length i. The complexity score is between 0 and 1 with low-scores indicating presence of repetitive combinations of letters in the sequence. For DNA sequence, the alphabet is the set of nucleotides. For computing complexity scores for the M. tuberculosis H37Rv reference genome, we split the genome in overlapping windows of length 21. Thus complexity score of a given site considers 10 positions upstream and downstream of the site in addition to the site itself.

**Computing Shannon’s entropy**

Shannon’s entropy quantifies the nucleotide diversity of a sequence from the frequencies of letters in the alphabet [2,3]. For a DNA sequence, Shannon’s entropy is defined as:

![H equation](/ShannonEntropy_eqn.png)<!-- .element height="10%" width="10%" -->

where pj is the frequency of nucleotide j in the sequence. Thus, pj = nj/n, where nj is the number of times nucleotide j appears in sequence of length n. We split the M. tuberculosis H37Rv reference genome in overlapping windows of size 21, same as for computing linguistic complexity. The Shannon’s entropy for a given genomic site this considered nucleotide frequencies from 10 bases upstream to 10 bases downstream of the site. By setting the logarithm base to 4, we obtained the entropy values between 0 and 1 with a homo-polymer sequence (no nucleotide diversity) having a score of 0 and a sequence with equal occurrence of all nucleotides getting a score of 1.

**References**
1.	Popov, O., Segal, D.M. & Trifonov, E.N. Linguistic complexity of protein sequences as compared to texts of human languages. Biosystems 38, 65-74 (1996).
2.	Adami, C. Information theory in molecular biology. Physics of Life Reviews 1, 3-22 (2004).
3.	Shannon, C.E. A Mathematical Theory of Communication. Vol. 27 379–423 (Bell System Technical
Journal, 1948).

