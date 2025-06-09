README.txt  
===========

Project Title
-------------
CpG Density and Genomic Feature Visualization in Plasmid DNA


Overview:
---------
This Python project visualizes CpG dinucleotide patterns, annotated genomic features, and nucleosome positioning signals in a plasmid DNA sequence provided in GenBank format (`.gbk`). It consists of:

1. Parsing a DNA sequence and extracting features from a GenBank file.
2. Identifying all CpG dinucleotide positions in the sequence.
3. Calculating CpG density using a sliding window approach.
4. Scanning for nucleosome-like DNA periodicity (AT-rich repeats every ~10 bp).
5. Visualizing the data in a 4-panel figure using `matplotlib`:
   - **Panel A**: CpG density (per 200 bp window).
   - **Panel B**: Individual CpG site positions.
   - **Panel C**: Annotated gene and feature locations.
   - **Panel D**: Predicted nucleosome affinity based on AT dinucleotide periodicity.


How It Works:
-------------
- **CpG Detection**: Finds all instances of the dinucleotide "CG" in the DNA.
- **Sliding Window Density**: Calculates how many "CG" pairs exist in every 200 bp window (stepping by 50 bp).
- **Feature Parsing**: Extracts labeled genomic features (e.g., `CDS`, `gene`, `misc_feature`) with names and coordinates.
- **Nucleosome Affinity**: Uses a basic heuristic to estimate nucleosome formation potential by scanning for "AA", "TT", "TA", or "AT" motifs every 10 bp in a 147 bp sliding window.
- **Visualization**: Displays four aligned plots showing patterns across the plasmid for analysis.


Requirements:
-------------
Install all required packages using the provided `requirements.txt`. Main dependencies:

- biopython
- matplotlib
- numpy


Usage:
------
1. Place your `.gbk` GenBank file and script in the same directory or update the file path.
2. Run the Python script.
3. A matplotlib figure will be shown with all four panels for analysis.

