# RIKEN_LIPIDOMICS
copy of All datasets (Excel) in https://metabography.riken.jp/menta.cgi/lipidomics/download_data_set

This work is licensed under a CC BY-NC 4.0 license.

`mztab/` holds one mzTab-M 2.1 file (profile M+S) for each `[[datasets]]` entry in `sample_metadata.toml`. Metadata comes from the TOML file, and the small-molecule table comes from the matching file in `parquet/`. Regenerate them with [RmzTabM](https://github.com/kozo2/RmzTabM):

```sh
Rscript scripts/export_mztabm.R sample_metadata.toml mztab
Rscript scripts/add_sml_from_parquet.R mztab parquet
```
