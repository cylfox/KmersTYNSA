# KmersTYNSA
## Description
Find seeds for alignments in long nucleotides sequences
### Methods/Functions
```python
def load_seq(sequence_file_name):
```
> Loads the desired sequence in memory
---

```python
def load_seq(sequence_file_name):
```
> 
---

```python
def change_char(string, position, replacement_char):
```
> 
---

```python
def kmer_index(sequence, kmer_size, alphabet):
```
> 
---
```python
def kmer_index2word(kmer_hash, kmer_size, alphabet):
```
> 
---
```python
def compute_frequency(sequence, kmer_size, alphabet):
```
> 
---
```python
def generate_dictionary(sequence, kmer_size, alphabet, o_file_name):
```
> 
---
```python
def group_dictionary(dic):
```
> 
---
```python
def generate_hits(i_file_name_1, i_file_name_2, o_file_name):
```
> 
---
```python
def fragment_extension_from_hits(hits_file_name, output_file_name, sequence_1, sequence_2, kmer_size):
```
> 
---
```python
def fragment_extension_from_hits2(hits_file_name, output_file_name, sequence_1, sequence_2, kmer_size):
```
> 
---
```python
def fragment_extension_from_hits_to_csv(hits_file_name, output_file_name, sequence_1, sequence_2, kmer_size):
```
> 
