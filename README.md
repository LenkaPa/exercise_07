# Motif Search

## The Brute-Force Motif Search
### Task 1
* In R, create a function `Score()`, that calculates the score for a consensus string.

* Input:
    * An array of starting indexes.
    * `DNAStringSet` object of sequences (for example file `seq_score.fasta`).
    * Motif length.
   
* Output:
    * The score for the consensus string.

#import
library(Biostrings)

#function
Score <- function(start_positions, dna_set, motif_length) {
  n_seq <- length(dna_set)

  # Validate indexing
  for (i in seq_len(n_seq)) {
    max_start <- width(dna_set[[i]]) - motif_length + 1
    if (start_positions[i] > max_start || start_positions[i] < 1) {
      stop(paste("Invalid start position for sequence", i, 
                 "- valid range is 1 to", max_start))
    }
  }

  motifs <- character(n_seq)
  for (i in seq_len(n_seq)) {
    motifs[i] <- as.character(subseq(dna_set[[i]],
                                     start = start_positions[i],
                                     width = motif_length))
  }

  bases <- c("A", "C", "G", "T")
  pfm <- matrix(0, nrow=4, ncol=motif_length,
                dimnames=list(bases, 1:motif_length))

  for (i in seq_len(n_seq)) {
    for (j in seq_len(motif_length)) {
      base <- substr(motifs[i], j, j)
      pfm[base, j] <- pfm[base, j] + 1
    }
  }

  score <- sum(apply(pfm, 2, max))
  return(score)
}


#example
motif_length <- 6
dna <- readDNAStringSet("/Users/lenapavlikova/Library/Mobile Documents/com~apple~CloudDocs/ŠKOLA/7. semestr/MPA-PRG/exercise_08/exercise_07/seq_score.fasta")

# Make sure these are within allowed range for your sequences!
start_positions <- rep(1, length(dna)) 
Score(start_positions, dna, motif_length)





 
### Task 2
* In R, create function `NextLeaf()` according to the following pseudocode.

* Input:
    * `s` An array of starting indexes `s = (s1 s2 … st)`, where *t* is the number of sequences.
    * `t` Number of sequences.
    * `k` `k = n - l + 1`, where `n` is length of sequences and `l` is motif length.

* Output:
    * `s` An array of starting indexes that corresponds to the next leaf in the tree.

```
NextLeaf(s, t, k)
1   for i ← t to 1
2     if s[i] < k
3       s[i] ← s[i] + 1
4       return s
5     s[i] ← 1
6   return s
```

### Task 3
* In R, create a function `BFMotifSearch()` according to the following pseudocode.

* Input:
    * `DNA` `DNAStringSet` object of sequences (for example file `seq_motif.fasta`).
    * `t` Number of sequences.
    * `n` Length of each sequence.
    * `l` Motif length.

* Output:
    * `bestMotif` An array of starting positions for each sequence with the best score for the consensus string.

```
BFMotifSearch(DNA, t, n, l)
1   s ← (1, 1, ... , 1)
2   bestScore ← Score(s, DNA, l)
3   while forever
4     s ← NextLeaf(s, t, n − l + 1)
5       if Score(s, DNA, l) > bestScore
6         bestScore ← Score(s, DNA, l)
7         bestMotif ← (s1, s2, . . . , st)
8       if s = (1, 1, . . . , 1)
9         return bestMotif
```

## The Branch-and-Bound Motif Search
### Task 4
* In R, create a function `NextVertex()` according to the following pseudocode.

* Input:
    * `s` An array of starting indexes `s = (s1 s2 … st)`, where *t* is the number of sequences.
    * `i` Level of vertex.
    * `t` Number of sequences.
    * `k` `k = n - l + 1`, where `n` is length of sequences and `l` is motif length.

* Output:
    * `s` The next vertex in the tree.
    * Current level of vertex.

```
NextVertex(s, i, t, k)
1   if i < t
2     s[i + 1] ← 1
3     return (s, i + 1)
4   else
5     for j ← t to 1
6       if s[j] < k
7         s[j] ← s[j] + 1
8         return (s, j)
9   return (s, 0)
```

### Task 5
* In R, create a function `ByPass()` according to the following pseudocode.

* Input:
    * `s = (s1 s2 … st)`; an array of starting indexes, where *t* is the number of sequences
    * `i`; level of vertex
    * `t`; number of DNA sequences
    * `k = n - l + 1`, where `n` is length of DNA sequences and `l` is motif length

* Output:
    * the next leaf after a skip of a subtree
    * current level of vertex

```
ByPass(s, i, t, k)
1   for j ← i to 1
2     if s[j] < k
3       s[j] ← s[j] + 1
4       return (s, j)
5   return (s, 0)
```

### Task 6
* In R, create a function `BBMotifSearch()` according to the following pseudocode.

* Input:
    * `DNA` `DNAStringSet` object of sequences (for example file `seq_motif.fasta`).
    * `t` Number of sequences.
    * `n` Length of each sequence.
    * `l` Motif length.
* Output:
    * `bestMotif` An array of starting positions for each sequence with the best score for the consensus string.

* Modify function `Score()` to calculate score for the consensus string of the first `i` sequences of `DNA`.

```
BBMotifSearch(DNA, t, n, l)
1   s ← (1, ... , 1)
2   bestScore ← 0
3   i ← 1
4   while i > 0
5     if i < t
6       optimisticScore ← Score(s, i, DNA, l) + (t - i) * l
7       if optimisticScore < bestScore
8         (s, i) ← ByPass(s, i, t, n - l + 1)
9       else
10        (s, i) ← NextVertex(s, i, t, n − l + 1)
11    else
12      if Score(s, t, DNA, l) > bestScore
13        bestScore ← Score(s, t, DNA, l)
14        bestMotif ← (s1, s2, ... , st)
15      (s, i) ← NextVertex(s, i, t, n − l + 1)
16  return bestMotif
```


<details>
<summary>Download files from GitHub</summary>
<details>
<summary>Basic Git settings</summary>

> * Configure the Git editor
>     ```bash
>     git config --global core.editor notepad
>     ```
> * Configure your name and email address
>     ```bash
>     git config --global user.name "Zuzana Nova"
>     git config --global user.email z.nova@vut.cz
>     ```
> * Check current settings
>     ```bash
>     git config --global --list
>     ```
>
</details>

* Create a fork on your GitHub account. 
  On the GitHub page of this repository find a <kbd>Fork</kbd> button in the upper right corner.
  
* Clone forked repository from your GitHub page to your computer:
    ```bash
    git clone <fork repository address>
    ```
* In a local repository, set new remote for a project repository:
```bash
git remote add upstream https://github.com/mpa-prg/exercise_07.git
```

#### Send files to GitHub
Create a new commit and send new changes to your remote repository.
* Add file to a new commit.
    ```bash
    git add <file_name>
    ```
* Create a new commit, enter commit message, save the file and close it.
    ```bash
    git commit
    ```
* Send a new commit to your GitHub repository.
    ```bash
    git push origin main
    ```

</details>
