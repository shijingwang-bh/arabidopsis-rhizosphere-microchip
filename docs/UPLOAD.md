# Upload this release to GitHub

Unpack the delivery ZIP. The folder `arabidopsis-rhizosphere-microchip` is the repository root. Submit its contents as repository files; a single uploaded ZIP is not an inspectable source repository.

The public repository is [shijingwang-bh/arabidopsis-rhizosphere-microchip](https://github.com/shijingwang-bh/arabidopsis-rhizosphere-microchip). For a fresh checkout, clone that repository. To initialize the unpacked release for an empty remote, use:

```sh
git init
git add .
git commit -m "Add manuscript data and code release 2026.11.25"
git branch -M main
git remote add origin https://github.com/shijingwang-bh/arabidopsis-rhizosphere-microchip.git
git push -u origin main
```

Use the existing `.gitattributes` to preserve measurement-file checksums across platforms. `validate_package.py` should pass after cloning or downloading. Commit `data/`, `reference/`, `code/`, `assets/`, metadata and documentation together. Reproduction output should be written to a new sibling directory and need not be committed.

The package is checked so every repository file is below 25 MiB. GitHub's documented browser-upload limit is 25 MiB per file; ordinary Git repositories warn above 50 MiB and reject files above 100 MiB. A Git client is convenient for uploading this multi-file package. See [GitHub's official large-file documentation](https://docs.github.com/en/repositories/working-with-files/managing-large-files/about-large-files-on-github). The external delivery ZIP is an interchange artifact and is not part of the repository contents.

The existing MIT license covers code. Existing data terms are retained separately in `DATA_REUSE_TERMS.txt`; they are not replaced by the code license. `CITATION.cff` records authors, release version and the repository URL. Add archive or manuscript DOIs when assigned. If any packaged files change, regenerate the checksum manifest for the new release and record the change; the present validator intentionally detects differences from this frozen snapshot.
