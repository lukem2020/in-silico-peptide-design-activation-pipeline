## 🧬 insilico peptide design

**🔗 [View Docking Results on Tamarind Bio](https://app.tamarind.bio/jobs/d724cdc7-c7c0-4405-9884-a33c9d79115d)**

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
├── run_docking.sh               # Example docking script
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

### End-to-End Pipeline

#### 1. Data Acquisition

**Objective**: Retrieve reference sequences and structures for both the target protein and the parent peptide.

- **Tasks**
  - **Get target protein sequence** (e.g., GLP-1R or your chosen receptor) from UniProt as FASTA.
  - **Optionally get the parent peptide/ligand** (e.g., the natural hormone) sequence.
  - **Download 3D structure** of the target from RCSB PDB or AlphaFold DB.
  - **Store everything** in a standardized `data/` directory layout.

#### 2. Structure Preparation

**Objective**: Convert raw protein/peptide structures into models ready for modeling and docking.

- **Tasks**
  - **Clean structures** (remove water, ions, alternate locations).
  - **Add hydrogens and protonate** at physiological pH.
  - **Minimize the 3D structure** of the parent peptide.
  - **Export to required formats** (PDB, PDBQT, MOL2).

#### 3. Peptide Library Design

**Objective**: Systematically generate sequence variants derived from a parent peptide.

- **Tasks**
  - **Define which positions to mutate** and which amino acids are allowed at each position.
  - **Enumerate all combinatorial variants** consistent with these rules.
  - **Compute basic biophysical descriptors**:
    - Net charge  
    - Average hydrophobicity (Kyte–Doolittle)  
    - Sequence length
  - **Filter variants to remove**:
    - Highly charged sequences  
    - Extremely hydrophobic or extremely hydrophilic sequences  
    - Variants with unexpected length changes
  - **Export the selected variants** to FASTA for downstream 3D modeling.

#### 4. Docking Preparation

**Objective**: Prepare receptor and peptide variants for automated docking.

- **Tasks**
  - **Define a docking grid** centered on the known/putative binding site.
  - **Convert structures** (receptor and each variant peptide) to AutoDock Vina PDBQT format.
  - **Organize variant structures** into per-sequence directories.
  - **Create a batch docking script** to run Vina across all variants.

#### 5. Automated Docking

**Objective**: Evaluate each peptide variant’s binding to the target receptor.

- **Tasks**
  - **Run AutoDock Vina** for each receptor–peptide pair.
  - **Record outputs**, including:
    - Best docking score  
    - Binding poses  
    - Log files
  - **Store all results** in a structured directory tree.

#### 6. Post-Docking Analysis

**Objective**: Extract and summarize useful information from the docking results.

- **Tasks**
  - **Parse Vina log files** for each variant.
  - **Collect best-scoring poses**.
  - **Compute metrics**, such as:
    - Binding energies  
    - Ligand RMSD  
    - Optional: interaction features (hydrogen bonds, contacts, etc.)
  - **Score and rank all variants** based on these metrics.

#### 7. Refinement (Optional)

**Objective**: Improve physical realism via molecular mechanics or MD simulations.

- **Tasks**
  - **Perform quick minimization or short MD** for the top-ranked poses.
  - **Recalculate stability and interaction energies** for refined complexes.
  - **Re-rank candidates** with these more accurate estimates.

#### 8. Selection for Synthesis

**Objective**: Choose peptide variants for wet-lab synthesis and biological testing.

- **Tasks**
  - **Filter and prioritize hits** using:
    - Docking scores  
    - Sequence properties  
    - Binding pose quality  
    - Stability predictions
  - **Export final candidates** (e.g., as ranked CSV + FASTA) for experimental teams.

#### High-Level Summary in One Sentence

We generate peptide variants, model and dock them to a protein target, score their binding potential, and select the most promising candidates for synthesis and experimental validation.

### Notes for Reviewers and Employers

- This repository is meant as a **conceptual, end-to-end scaffold**, not a fixed benchmark for a single target.
- All steps are implemented as **scripts with configurable parameters** (e.g., choice of target, peptide, mutation rules, docking settings), so results are **not hard-coded**.
- Example configurations (such as specific UniProt IDs or mutation schemes) are provided only to **demonstrate the workflow**; in practice, the same pipeline can be re-run on new targets and design hypotheses.
