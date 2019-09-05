# KmersTYNSA
## Description
A bunch of methods to work in Python with nucleotides sequences using kmers hash.
### Methods/Functions
```python
def load_seq(sequence_file_name):
```
> Loads the desired sequence in memory.
---
```python
def change_char(string, position, replacement_char):
```
> Allows to change a letter for another in a string.
---

```python
def kmer_index(sequence, kmer_size, alphabet):
```
> Get the kmer index from a string of nucleotides of kmer size (there's 4 possibilities for each word in the nucleotides alphabet).
---
```python
def kmer_index2word(kmer_hash, kmer_size, alphabet):
```
> Translates a kmer index to the respective nucleotides.
---
```python
def compute_frequency(sequence, kmer_size, alphabet):
```
> Index a whole sequence and calculates the frecuency of each index.
---
```python
def generate_dictionary(sequence, kmer_size, alphabet, o_file_name):
```
> Generates a dictionary, sort it by kmer index, group them for easier calculations and dump it to a file in pickle format.
---
```python
def group_dictionary(dic):
```
> Given a sorted dictionary makes a new one grouping positions by the kmer hash (eg. ACCT, [12, 231, 235]).
---
```python
def generate_hits(i_file_name_1, i_file_name_2, o_file_name):
```
> Compares two dictionaries files to get the hits between them and dump it to a file in pickle format.
---
```python
def fragment_extension_from_hits(hits_file_name, output_file_name, sequence_1, sequence_2, kmer_size):
```
> Find the seeds and extend them to get and simplify new fragments and dumps it to a file in pickle format to later use. This approach uses the sequence length as a score to select the fragments.
---
```python
def fragment_extension_from_hits2(hits_file_name, output_file_name, sequence_1, sequence_2, kmer_size):
```
> Find the seeds and extend them to get and simplify new fragments and dumps it to a file in pickle format to later use. This approach uses the best number of hits as a score to select the fragments.
---
```python
def fragment_extension_from_hits_to_csv(hits_file_name, output_file_name, sequence_1, sequence_2, kmer_size):
```
> Find the seeds and extend them to get and simplify new fragments and dumps it to a csv file. This approach uses the best number of hits as a score to select the fragments.
