# ResiCritic

A hybrid command-line tool with file dialogs for analyzing mass spectrometry data from Proteome Discoverer (Thermo Fisher) or PEAKS (Bioinformatics Solutions).

## Building

It is possible to build a standalone executable of ResiCritic by using the provided `build.bat`/`build.sh` scripts or by running the following commands:
```sh
# Optional - create a virtual environment
python -m venv .venv
./.venv/Scripts/activate
# Note: "activate" script may appear in .venv/bin, and the appropriate script may have a different file extension depending on your operating system

# Install requirements
pip install -r requirements.txt

# (Optional) Create a standalone executable
pyinstaller --onefile --console --name resicritic __main__.py
```
The build scripts additionally install `uv` as a faster alternative to `pip`.

## How to run

An executable version of ResiCritic for 64-bit Windows is available in the [Releases](https://github.com/AlgebraManiacABC/ResiCritic/releases) section of the project,
or built by following the steps above. Alternatively, the script can be run in python after installing the necessary requirements (shown above):

```
python __main__.py
```

## How to use

ResiCritic operates inside a command line, but uses a file dialog for any file loading/saving. Here is a summary of operations:
1. ResiCritic prints welcome message with version number, then asks the user for an Excel file to analyze
2. User selects Excel document to load, then selects which sheet the data is located using arrow keys (enter to submit)
3. If the user has analyzed this document before, a `.pkl` file has been saved logging the progress, and the user may return to any step in the process.
4. Depending on the software and method, the Excel file may contain multiple protein analyses. If only one is found, a confirmation message appears.
5. By default, the program looks for "Acetyl" modifications. The user may enter any modification as a string (e.g. Methyl, Oxidation, etc.) as it appears in the Excel
6. The program then proceeds through each column with "PSMs" or "#Scans", and requests a FASTA file for each, provided by the user. If a FASTA has already been loaded, the user may select it from the list of loaded FASTAs rather than re-using the file dialog. In some cases, a column will contain data for multiple peptides. In such cases, multiple FASTAs may be selected.
7. Once all columns have been prepared, the program presents the user with confirmation of the FASTA selections. The user may at this time edit a single column's chosen FASTAs, completely replace a FASTA file for another, or restart the FASTA selection completely.
8. Once the user is satisfied and continues to analysis, the script computes the coverages and modification ratios for the given columns and FASTAs.
9. The user is prompted to save the coverages as a `.txt` file (or as a `.zip` file containing all coverages)
10. The user is prompted to save the modification ratio files as a `.tsv` (tab-separated values) file (or as a `.zip` file containing all coverages)
11. The program completes! Note the saved files, as well as the `.pkl` file which acts as a progress tracker and a sort of verbose log.

## Acknowledgements

This work is part of a project which received support from Environmental Molecular Sciences Laboratory project 61285, a Univ. of Missouri College of Agriculture Joy of Discovery award, DOE grant DE-SC0023142 to JJT, and Hatch Multi-State Project No. 1013013 from the National Institute of Food and Agriculture. Initially cryo-electron microscopy screening was supported by the Kansas University Medical Center Cryo-electron microscopy Facility and the Glacios equipment purchase via S10OD036339-01.
