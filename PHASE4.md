# Phase 4 — Extensions & Optimization

## 4.1 Isotope Pattern Scoring

**Goal**: Score candidates by comparing theoretical isotope distribution against experimental data.

### Input
- Experimental isotope peaks: `Vec<(f64, f64)>` — (m/z, relative_intensity)
- Provided via CLI: `--exp-pattern 506.9958:100,507.9991:13.3,509.0025:3.5`
- Or from file: `--exp-pattern-file spectrum.csv`

### Scoring Methods

1. **Cosine similarity** (default)
   ```
   score = Σ(theo_i × exp_i) / (||theo|| × ||exp||)
   ```
   Match peaks by nearest nominal mass (M+0, M+1, M+2, ...).

2. **Chi-squared**
   ```
   χ² = Σ((exp_i - theo_i)² / theo_i)
   ```
   Lower is better. Good for penalizing large deviations.

3. **Sigma-fit** (Xcalibur-style)
   ```
   σ = sqrt(Σ((exp_i - theo_i)²) / n)
   ```

### Output
- Add `isotope_score` field to `FormulaCandidate`
- New table column "Score" when experimental pattern is provided
- Results sorted by score (best match first) instead of mass error

### Implementation: `src/scoring.rs`

```rust
pub fn cosine_similarity(
    theoretical: &[IsotopePeak],
    experimental: &[(f64, f64)],
) -> f64;

pub fn chi_squared(
    theoretical: &[IsotopePeak],
    experimental: &[(f64, f64)],
) -> f64;
```

---

## 4.2 Extended Element Support

**Goal**: Add F, Cl, Br, I, Si to the search engine with proper handling.

### Already Done
- Element data is defined in `elements.rs` (F, Cl, Br, I, Si)
- `parse_element_string()` handles multi-char symbols (Cl, Br, Si)

### Remaining Work
- **Halogen-aware pruning**: Cl and Br produce distinctive M+2 patterns.
  When experimental isotope pattern is available:
  - M+2 ≈ 32% of M+0 → likely 1 Cl
  - M+2 ≈ 97% of M+0 → likely 1 Br
  - Use this to constrain Cl/Br count before enumeration
- **Element ordering**: Halogens (mass 35-127) should be placed in outer
  recursion layers for better pruning. Current sort-by-mass already handles this.

---

## 4.3 Fine Isotope Structure

**Goal**: Distinguish contributions of different elements to M+1, M+2 peaks.

### Use Case
At sufficient resolving power (Orbitrap > 100k), M+1 peak splits into:
- ¹³C contribution at M + 1.00335
- ¹⁵N contribution at M + 0.99703
- ³³S contribution at M + 0.99939

### Implementation
- Keep mass/probability pairs separate during convolution (don't merge by nominal mass)
- Add `--fine-isotope` flag
- Output sub-peaks within each nominal mass group
- Useful for distinguishing formulas that have similar M+1 intensity but different composition

---

## 4.4 Database Matching

**Goal**: Cross-reference candidates against metabolite databases.

### Data Sources
- **HMDB** (Human Metabolome Database): ~220k metabolites with exact masses
- **KEGG COMPOUND**: ~18k biologically relevant compounds
- **PubChem**: broader coverage, millions of compounds
- **LipidMAPS**: specialized for lipids

### Implementation: `src/database.rs`

**Approach 1 — Embedded lookup table (fast, offline)**
- Ship a binary file with (exact_mass, formula, name, db_id) tuples
- Binary search by mass within tolerance
- ~5 MB for HMDB core metabolites

**Approach 2 — REST API queries (comprehensive, online)**
- Query HMDB/PubChem APIs by exact mass
- Cache results locally
- Slower but always up-to-date

### Output
- Add `database_hits` field to `FormulaCandidate`
- New columns: "DB Match", "Compound Name"
- Flag candidates that match known compounds

### CLI
```
mfcalc 180.06339 -p 2 --db hmdb
mfcalc 505.9885 -m neg -a "[M-H]-" --db hmdb,kegg
```

---

## 4.5 MS/MS Integration (Future)

**Goal**: Use fragmentation data to rank or filter candidates.

### Approach
- Accept MS/MS spectrum as input (precursor m/z + fragment list)
- For each candidate formula:
  - Generate plausible fragment formulas (combinatorial sub-formulas)
  - Match against observed fragments within tolerance
  - Score by number of explained fragments

### This is a large feature — defer to a separate project or integrate with
existing tools (CFM-ID, MetFrag, SIRIUS).

---

## 4.6 Performance: SIMD Optimization

**Goal**: Accelerate isotope convolution with SIMD intrinsics.

### Where It Helps
- Isotope distribution calculation is the bottleneck for batch + isotope mode
- Convolution inner loop is a multiply-accumulate pattern → ideal for SIMD

### Approach
- Use `std::simd` (nightly) or `packed_simd2` crate
- Vectorize the inner loop of `convolve()`
- Expected speedup: 2-4x on isotope calculation

### Priority: Low — current performance is already excellent (<1s for 931 masses)

---

## 4.7 Python Bindings (PyO3)

**Goal**: Expose core functionality to Python for integration with metabolomics pipelines.

### API Surface
```python
import mfcalc

# Simple search
results = mfcalc.search(mass=180.06339, ppm=2.0, elements="CHNOPS")
for r in results:
    print(r.formula, r.mass_error_ppm, r.dbe)

# With adduct
results = mfcalc.search(mz=505.9885, adduct="[M-H]-", ppm=2.0)

# Isotope distribution
pattern = mfcalc.isotope_pattern("C6H12O6")

# Batch
results = mfcalc.search_batch(masses=[180.06, 505.99, 744.08], ppm=2.0)
```

### Implementation
- Add `pyo3` feature flag in Cargo.toml
- Create `src/python.rs` with `#[pymodule]`
- Build with `maturin develop`
- Publish to PyPI as `mfcalc`

---

## 4.8 WebAssembly / GUI

**Goal**: Run in browser or as desktop app.

### WASM
- Compile core library to `wasm32-unknown-unknown`
- Wrap with JavaScript/TypeScript bindings
- Embed in a web UI (React/Svelte)

### Desktop (Tauri)
- Rust backend (already done) + web frontend
- Native performance with small binary size
- Cross-platform (macOS, Linux, Windows)

---

## Implementation Priority

| Feature | Impact | Effort | Priority |
|---------|--------|--------|----------|
| 4.1 Isotope scoring | High | Medium | **P0** |
| 4.2 Extended elements | Medium | Low | **P0** (mostly done) |
| 4.4 Database matching | High | Medium | **P1** |
| 4.3 Fine isotope | Medium | Medium | **P2** |
| 4.7 Python bindings | High | Medium | **P2** |
| 4.5 MS/MS integration | High | High | **P3** |
| 4.8 WASM/GUI | Medium | High | **P3** |
| 4.6 SIMD optimization | Low | Medium | **P4** |
