## 🧬 insilico peptide design

Many therapeutic peptides function by binding to protein targets, inducing conformational changes that result in activation or inhibition of those proteins.

**Examples include:**
- **GLP-1 analogues** (used in diabetes treatment)
- **Peptide hormones** (involved in signal transduction)
- **Protein–protein interaction mimetics**

Such interactions are inherently structural. As a result, *in silico* design methods provide an efficient way to explore potential peptide candidates before committing to costly laboratory experiments.

### 📦 Pipeline Structure

This repository contains modular scripts implementing each major step of the *in silico* peptide design process:

<pre>

src/
│
├── fetch_data.py                # Download protein sequence & structure
├── prepare_structures.py        # Clean & prepare protein models
├── design_library.py            # Generate peptide sequence libraries
├── docking_prep.py              # Prepare files for docking tools
├── run_docking.sh               # Example docking script (e.g., Vina/HADDOCK)
├── parse_docking.py             # Extract poses & scores
├── scoring.py                   # Rank peptides by binding metrics
├── refine.py                    # MD-based refinement of top peptides
├── select_for_synthesis.py      # Automatic hit selection
├── utils.py                     # Shared helper functions
│
├── data/
│   └── GLP1_template.fasta      # Example peptide template
│
└── README.md                    # This file

</pre>


Each step is implemented as a separate script or notebook.
