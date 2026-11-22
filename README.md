# Sequence Motifs: PWM and Information Content

R script that builds several position-specific matrices from a set of aligned DNA sequences (the POU2F1 binding site example).

## What it does

From 8 aligned sequences of 8 nucleotides, the script computes and prints:

1. **Counts matrix (PPSM):** how many times each nucleotide appears at each position.
2. **Position weight matrix (PWM):** counts divided by the column total.
3. **PWM with pseudocounts (p = 1):** avoids zero probabilities.
4. **Log-likelihood matrix:** natural log of the pseudocount PWM divided by the background probability (0.25).
5. **Information content matrix:** bits of information per nucleotide and position.

## Requirements

- R (tested with R 4.3.3)

No extra packages or input files are needed.

## Usage

Open `MOTIFS_SCRIPT.R` in RStudio and run the whole script (`Source`), or from the terminal, in the repository folder:

```bash
Rscript MOTIFS_SCRIPT.R
```

To use other sequences, edit the `motifs` vector at the top of the script. All sequences must have the same length. The pseudocount value can be changed in the `pseudocounts` variable.

## Files

- `MOTIFS_SCRIPT.R`: main script

## Author

Èric Vargas Riera
