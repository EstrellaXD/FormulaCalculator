# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] - 2026-02-24

### Added
- Recursive DFS molecular formula search with H direct-calculation and suffix-max pruning
- Seven Golden Rules filtering (RDBE, element ratios loose/strict, nitrogen rule)
- Isotope distribution calculation via repeated-squaring convolution
- Adduct support: 10 common adducts for positive and negative ion modes
- Auto-adduct mode: `-m pos`/`-m neg` automatically tries all adducts and merges results
- CHNOPS + extended elements (F, Cl, Br, I, Si)
- CLI with table, JSON, CSV output formats
- Batch mode for processing mass lists from files
- Rayon-based parallel search (split by C count)
- 15 integration tests covering 9 reference compounds
