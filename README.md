# EGM

Analysis code for "Induction of rejuvenation reprogramming in aged mice and non-human primates by defined plasma factors" (Kim J, Hwang Y, Cho B, et al.).

## Layout

- `Figure1/` – `Figure5/`, `SuppFigure*/`: scripts for each figure. File names start with the panel they produce (e.g. `Fig3D_...` → Fig. 3D).
- `shared/`: upstream processing used by more than one figure (bulk RNA-seq, single-cell preprocessing and CellOracle, hdWGCNA, chromatin profiling, enrichment plots).

## Notes

- Code only; data are not included.
- File paths are relative; adjust them to the location of the data.
- Written in R, Python and bash.

## Affiliation

Institute for Stem Cells and Regenerative Medicine (ISR), Department of Chemistry, Dongguk University, Seoul 04620, Republic of Korea.

## Contact

Jongpil Kim (jpkim153@dongguk.edu)

## License

Copyright © 2026 Institute for Stem Cells and Regenerative Medicine (ISR), Dongguk University. Released under the MIT License; see [LICENSE](LICENSE).
