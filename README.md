## Python code for the bioinformatics analyses from Alderson, Adriaenssens, et al.

### Searches a proteome for specific motifs and performs statistical tests.<br /> 

Given a supplied proteome and queried peptide motif(s), the script will search for the number of instances of the motif(s) and calculate the expected number of instances based on amino acid frequencies in the proteome. If a list of disordered regions is supplied, then a comparison between structured and disordered regions of the proteome can be perofrmed. Flexible input within the amino-acid motif is permitted; for example, "X" allows any residue at that position and multiple residues indicate "or". Residue positions are separated by the _ character, and so the IxI/V motif would be inputted as, "IV_X_IV". <br />

To execute this script, the user must have **biopython**, **itertools**, **numpy**, **pandas**, and **scipy** installed in their Python distribution. The user must provide a proteome file downloaded from [UniProt](https://www.uniprot.org/) or in FASTA format. For example outputs, see the images below: <br />
<br />
* Statistical tests to check for deviations from expectations based on amino-acid frequencies in the proteome
<p align="left">
  <img src="figures/structured_chisq.png" width="400px" height="auto"/>
  <img src="figures/disordered_chisq.png" width="400px" height="auto"/>
</p>
<br />

* Comparison of IxI/V motif frequency in structured and disordered regions of the proteome <br />
<p align="left">
  <img src="output/compare_fractions_IV_X_IV_motifs.png" width="908px" height="auto"/>
</p>

See these pages for more information or installation of the various Python packages:<br />
[biopython](https://biopython.org/wiki/Download)<br />
[itertools](https://docs.python.org/3/library/itertools.html)<br />
[numpy](https://docs.scipy.org/doc/numpy-1.10.1/user/install.html) <br />
[numpy via a pre-built package](https://scipy.org/install.html) <br />
[matplotlib](https://matplotlib.org/faq/installing_faq.html)<br />
[pandas](https://pypi.org/project/pandas/)<br />
[scipy](https://www.scipy.org/install.html)<br />

