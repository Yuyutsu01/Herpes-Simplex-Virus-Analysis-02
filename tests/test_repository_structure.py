"""
Test Suite: test_repository_structure.py
Purpose: Automated structural, linkage, and integrity validation for the reorganized repository.
Enforces that:
  1. Canonical scripts exist, are descriptive, and cover phases 01 through 16.
  2. Required raw datasets exist and match SHA-256 hashes.
  3. Critical processed artifacts and manifests exist and are non-empty.
  4. Manuscript figures and tables exist in figures/, tables/, and manuscript/.
  5. Internal markdown file links in docs/ and manuscript/ resolve to valid files.
  6. Documentation governance files exist (README, LICENSE, CITATION.cff, requirements.txt, .gitignore).
  7. No obvious secrets or sensitive keys are exposed.
  8. Scientific identifiers (UniProt accessions, homology groups, 687 benchmark rows) remain intact.
"""

import os
import re
import hashlib
from pathlib import Path
import pytest
import pandas as pd

BASE_DIR = Path(__file__).resolve().parent.parent

# 1. Test canonical scripts existence and phase coverage
def test_canonical_scripts_and_phase_coverage():
    """Verify scripts in scripts/ follow descriptive naming and cover phases 01 through 16."""
    scripts_dir = BASE_DIR / "scripts"
    assert scripts_dir.exists() and scripts_dir.is_dir()
    
    py_scripts = list(scripts_dir.glob("*.py"))
    assert len(py_scripts) > 0, "No scripts found in scripts/"
    
    phases_present = set()
    for script in py_scripts:
        # Match descriptive names: e.g. 01_fasta_inventory.py, 16_generate_main_figures.py
        match = re.match(r"^(\d+)_([a-zA-Z0-9_]+)\.py$", script.name)
        if match:
            phase_num = match.group(1)
            action_subject = match.group(2)
            assert len(action_subject) > 3, f"Script name not descriptive: {script.name}"
            phases_present.add(phase_num)

    # Assert canonical analytical phases are represented
    for expected in ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "15", "16"]:
        assert expected in phases_present, f"Missing script representation for phase {expected}"


# 2. Test raw datasets existence and cryptographic integrity
def test_raw_datasets_and_checksums():
    """Verify that raw FASTA files exist and match their authoritative SHA-256 checksums."""
    raw_dir = BASE_DIR / "data" / "raw"
    assert raw_dir.exists(), "Missing data/raw directory"
    
    checksum_file = raw_dir / "RAW_DATA_CHECKSUMS.sha256"
    assert checksum_file.exists(), "Missing data/raw/RAW_DATA_CHECKSUMS.sha256"
    
    with open(checksum_file, "r", encoding="utf-8") as f:
        lines = f.readlines()
        
    verified_files = 0
    for line in lines:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        expected_hash = parts[0]
        rel_path = parts[1].replace("\\", "/")
        
        target_file = BASE_DIR / rel_path
        assert target_file.exists(), f"Raw file listed in checksum manifest missing: {target_file}"
        
        # Calculate SHA-256 hash
        hasher = hashlib.sha256()
        with open(target_file, "rb") as f_in:
            while chunk := f_in.read(65536):
                hasher.update(chunk)
        actual_hash = hasher.hexdigest()
        assert actual_hash == expected_hash, f"Hash mismatch for {rel_path}: expected {expected_hash}, got {actual_hash}"
        verified_files += 1

    assert verified_files == 2, f"Expected 2 verified raw files, found {verified_files}"


# 3. Test critical processed datasets and manifests
def test_critical_processed_artifacts_exist():
    """Assert critical processed datasets, embeddings, and manifests exist and are non-empty."""
    processed_dir = BASE_DIR / "data" / "processed"
    m_supp_dir = BASE_DIR / "manuscript" / "supplementary_materials"
    
    critical_processed = [
        processed_dir / "curated_protein_candidates.csv",
        processed_dir / "valid_proteins_metadata.csv",
        processed_dir / "esm2_t12_representatives_embeddings.npy",
        processed_dir / "esm2_t12_representatives_metadata.csv",
        processed_dir / "08_temporal_label_adjudication.csv",
        processed_dir / "08_temporal_gene_evidence.csv",
        processed_dir / "08_gene_alias_dictionary.json",
        processed_dir / "08_1_strict_temporal_benchmark.csv",
        processed_dir / "09_group_split_manifest.csv",
        m_supp_dir / "table_s1_supervised_benchmark_687.csv",
        m_supp_dir / "table_s2_split_manifest_5fold.csv",
        m_supp_dir / "table_s3_family_oof_predictions.csv",
        m_supp_dir / "table_s4_sequence_oof_predictions.csv",
    ]
    for filepath in critical_processed:
        assert filepath.exists(), f"Missing critical processed artifact: {filepath}"
        assert filepath.stat().st_size > 0, f"Critical artifact is empty: {filepath}"


# 4. Test manuscript figures and tables
def test_manuscript_figures_and_tables_exist():
    """Verify that main publication figures and tables exist and are non-empty in primary & manuscript dirs."""
    figures_dir = BASE_DIR / "figures"
    m_figures_dir = BASE_DIR / "manuscript" / "figures"
    m_tables_dir = BASE_DIR / "manuscript" / "tables"
    m_supp_dir = BASE_DIR / "manuscript" / "supplementary_materials"
    
    # Check figures 1 through 8
    for i in range(1, 9):
        matched_pdfs = list(figures_dir.glob(f"figure_{i}_*.pdf"))
        assert len(matched_pdfs) > 0, f"Missing vector PDF for Figure {i} in figures/"
        matched_pngs = list(figures_dir.glob(f"figure_{i}_*.png"))
        assert len(matched_pngs) > 0, f"Missing PNG for Figure {i} in figures/"
        m_pdfs = list(m_figures_dir.glob(f"figure_{i}_*.pdf"))
        assert len(m_pdfs) > 0, f"Missing vector PDF for Figure {i} in manuscript/figures/"
        
    # Check main tables 1 through 7 in manuscript/tables
    for i in range(1, 8):
        m_tables = list(m_tables_dir.glob(f"table_{i}_*.csv"))
        assert len(m_tables) > 0, f"Missing CSV for Table {i} in manuscript/tables/"

    # Check supplementary tables S1 through S8
    for i in range(1, 9):
        matched_supp = list(m_supp_dir.glob(f"table_s{i}_*.csv"))
        assert len(matched_supp) > 0, f"Missing Supplementary Table S{i} in manuscript/supplementary_materials/"


# 5. Test documentation files existence
def test_core_documentation_files_exist():
    """Verify all required documentation hub files, guides, and ADRs exist."""
    docs_dir = BASE_DIR / "docs"
    required_docs = [
        docs_dir / "README.md",
        docs_dir / "PROJECT_FILE_GUIDE.md",
        docs_dir / "REPRODUCIBILITY.md",
        docs_dir / "DATA_DICTIONARY.md",
        docs_dir / "SCIENTIFIC_SCOPE_AND_LIMITATIONS.md",
        docs_dir / "SUBMISSION_PACKAGE.md",
    ]
    for doc in required_docs:
        assert doc.exists(), f"Missing core documentation: {doc}"
        assert doc.stat().st_size > 100, f"Documentation file suspiciously small: {doc}"

    # Check 10 phases exist
    phases_dir = docs_dir / "phases"
    for i in range(1, 11):
        phase_files = list(phases_dir.glob(f"PHASE_{i:02d}_*.md"))
        assert len(phase_files) == 1, f"Missing or duplicate document for Phase {i:02d}"

    # Check decision records exist
    decisions_dir = docs_dir / "decisions"
    for i in range(1, 6):
        decision_files = list(decisions_dir.glob(f"DECISION_{i:02d}_*.md"))
        assert len(decision_files) == 1, f"Missing or duplicate document for Decision {i:02d}"


# 6. Test root governance files
def test_root_governance_files_exist():
    """Verify README.md, LICENSE, CITATION.cff, requirements.txt, and .gitignore exist."""
    for filename in ["README.md", "LICENSE", "CITATION.cff", "requirements.txt", ".gitignore"]:
        p = BASE_DIR / filename
        assert p.exists(), f"Missing root governance file: {filename}"
        assert p.stat().st_size > 0, f"Root governance file empty: {filename}"


# 7. Test markdown links integrity
def test_internal_markdown_file_links_resolve():
    """Scan docs/ for file:/// links and verify target existence."""
    link_pattern = re.compile(r"\[.*?\]\(file:///([^\)\#]+)(?:#[^\)]*)?\)")
    docs_dir = BASE_DIR / "docs"
    broken_links = []
    
    for md_file in docs_dir.rglob("*.md"):
        content = md_file.read_text(encoding="utf-8", errors="ignore")
        matches = link_pattern.findall(content)
        for target in matches:
            target_clean = target.replace("\\", "/").rstrip("/")
            if ":" in target_clean:
                p = Path(target_clean)
            else:
                p = BASE_DIR / target_clean
                
            if not p.exists():
                broken_links.append((md_file.name, target))
                
    assert len(broken_links) == 0, f"Found broken internal markdown links: {broken_links[:10]}"


# 8. Test secret scanning (no sensitive keys/tokens)
def test_no_secrets_in_tracked_files():
    """Assert that no API keys, private credentials, or auth tokens exist in codebase."""
    secret_patterns = [
        re.compile(r"(?i)api[_-]?key\s*=\s*['\"][A-Za-z0-9_\-]{16,}['\"]"),
        re.compile(r"(?i)bearer\s+[A-Za-z0-9_\-\.]{20,}"),
        re.compile(r"-----BEGIN (RSA|EC|DSA|OPENSSH|PRIVATE) KEY-----"),
    ]
    
    for ext in ["*.py", "*.md", "*.json", "*.yaml", "*.yml"]:
        for fpath in BASE_DIR.rglob(ext):
            if any(part in fpath.parts for part in ["venv", ".git", ".pytest_cache", "__pycache__"]):
                continue
            text = fpath.read_text(encoding="utf-8", errors="ignore")
            for pattern in secret_patterns:
                match = pattern.search(text)
                assert not match, f"Potential sensitive secret found in {fpath}: {match.group(0)}"


# 9. Test scientific identifiers integrity
def test_scientific_identifiers_and_benchmark_integrity():
    """Assert benchmark dataset maintains 687 sequences, 63 families, and fold balance."""
    supp_dir = BASE_DIR / "manuscript" / "supplementary_materials"
    benchmark_path = supp_dir / "table_s1_supervised_benchmark_687.csv"
    assert benchmark_path.exists(), "Benchmark dataset table_s1 missing"
    
    df = pd.read_csv(benchmark_path)
    assert len(df) == 687, f"Expected 687 benchmark sequences, got {len(df)}"
    
    # Assert 63 unique homology groups
    families = df["effective_gene_symbol"].unique()
    assert len(families) == 63, f"Expected 63 homology groups, got {len(families)}"
    
    # Assert 4 kinetic classes are present
    classes = set(df["normalized_temporal_class"].unique())
    assert classes == {"IE_ALPHA", "E_BETA", "L1_GAMMA1", "L2_GAMMA2"}, f"Unexpected classes: {classes}"
    
    # Assert split manifest maintains 5 folds and 687 rows
    splits_path = supp_dir / "table_s2_split_manifest_5fold.csv"
    df_splits = pd.read_csv(splits_path)
    assert len(df_splits) == 687, f"Expected 687 rows in split manifest, got {len(df_splits)}"
    assert set(df_splits["fold_index"].unique()) == {0, 1, 2, 3, 4}, "Folds do not match 0-4"
    
    # Assert family out-of-fold predictions has 315 rows (63 families x 5 models)
    oof_fam_path = supp_dir / "table_s3_family_oof_predictions.csv"
    df_oof = pd.read_csv(oof_fam_path)
    assert len(df_oof) == 315, f"Expected 315 family predictions, got {len(df_oof)}"
