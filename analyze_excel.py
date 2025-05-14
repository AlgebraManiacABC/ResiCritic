import os
import sys
import pandas as pd
import tempfile
import yaml
from yaml import SafeDumper, SafeLoader
import questionary as q
from questionary import Validator, ValidationError
import re
import tkinter as tk
from tkinter import filedialog, messagebox
from rich import print
from pathlib import Path
import math
import unicodedata
import zipfile
import shutil

VERSION_MAJOR = 1
VERSION_MINOR = 2
print(f"Welcome to the [bold green]Mass Spec Ratio Analyzer![/] (Version {VERSION_MAJOR}.{VERSION_MINOR})")
print(f"Developed by [bold purple]Luke Dunn[/] at the [bold orange3]University of Missouri[/]")

def select_file(title="Select file", filetypes=None):
    root = tk.Tk()
    root.withdraw()
    if filetypes:
        file_path = filedialog.askopenfilename(title=title,filetypes=filetypes)
    else:
        file_path = filedialog.askopenfilename(title=title)
    root.destroy()
    return file_path

def sanitize_filename(filename: str) -> str:
    # Normalize Unicode, convert to ASCII
    filename = unicodedata.normalize("NFKD", filename).encode('ascii', 'ignore').decode('ascii')
    # Remove badness
    filename = re.sub(r'[<>:"/\\|?*\x00-\x1F]', '_', filename)
    # Collapse and replace whitespace
    filename = re.sub(r'\s+', ' ', filename)
    # Strip whitespace and trailing dots
    filename = filename.strip(" .")
    # Limit length
    filename = filename[:255]
    # Done!
    return filename

# Expects the Excel file as the first and only argument
# Otherwise, have the user select the file.
filename = ""
if len(sys.argv) != 2:
    while not filename.lower().endswith('.xlsx'):
        messagebox.showinfo("Information", "Please select an Excel file to analyze.")
        filename = select_file("Open Excel File to Analyze",[("Excel files","*.xlsx")])
        if not filename or filename == "":
            sys.exit(-1)
else:
    filename = sys.argv[1]
    if not filename.lower().endswith('.xlsx'):
        print(f"'{filename}' was not a valid Excel file.")
        sys.exit(-1)

file = None
try:
    file = pd.ExcelFile(filename)
except FileNotFoundError:
    messagebox.showerror(f"{filename} does not exist!")
except ValueError as ve:
    messagebox.showerror(f"There was an error parsing the Excel file '{filename}': {ve}")
except Exception as e:
    messagebox.showerror(f"An unexpected error occurred: {e}")
finally:
    if not file:
        sys.exit(-1)

# Grab the first sheet (we are assuming the data is here)
df = file.parse(file.sheet_names[0])
file.close()

class Peptide:
    def __init__(self, seq: str, mod_list: list[str], counts: dict[str, int], mod_dict=None, fasta_match_indices=None):
        self.sequence = seq
        self.mod_list = mod_list
        self.counts = counts
        self.mod_dict = mod_dict if mod_dict is not None else {}
        self.fasta_match_indices = fasta_match_indices if fasta_match_indices is not None else {}
            # Will be dict of str (fasta path) to int (index of peptide in fasta)

    def __repr__(self):
        pep = f"Peptide: '{self.sequence}'\nMods (filtered): {self.mod_list}\nCounts: {self.counts}"
        if self.mod_dict:
            pep += f"\nMod Dictionary: {repr(self.mod_dict)}"
        if self.fasta_match_indices:
            pep += f"\nFasta Match Indices: {repr(self.fasta_match_indices)}"
        return pep

def peptide_yaml_representer(dumper, peptide: Peptide):
    return dumper.represent_mapping('!Peptide', {
        'sequence': peptide.sequence,
        'mod_list': peptide.mod_list,
        'counts': peptide.counts,
        'mod_dict': peptide.mod_dict,
        'fasta_match_indices': peptide.fasta_match_indices
    })

def peptide_yaml_constructor(loader, node):
    values = loader.construct_mapping(node, deep=True)
    peptide = Peptide(
        values['sequence'],
        values['mod_list'],
        values['counts'],
        values['mod_dict'],
        values['fasta_match_indices']
    )
    return peptide

yaml.add_representer(Peptide, peptide_yaml_representer, Dumper=SafeDumper)
yaml.add_constructor('!Peptide', peptide_yaml_constructor, Loader=SafeLoader)

class Ratio:
    '''
    mod_count: count found modified
    count: total count found
    ratio (method): mod_count/count
    '''
    def __init__(self, mod_count: int, count: int):
        self.mod_count = mod_count
        self.count = count
    
    # Custom ratio (i.e. - division) function to return `nan` on divByZero
    def ratio(self):
        return float(self.mod_count) / float(self.count) if self.count else 'nan'

    def __repr__(self):
        return f"{self.mod_count}:{self.count} ({self.ratio()})"

def ratio_yaml_representer(dumper, ratio: Ratio):
    return dumper.represent_mapping('!Ratio', {
        'mod_count': ratio.mod_count,
        'count': ratio.count
    })

def ratio_yaml_constructor(loader, node):
    values = loader.construct_mapping(node, deep=True)
    ratio = Ratio(values['mod_count'],values['count'])
    return ratio

yaml.add_representer(Ratio, ratio_yaml_representer, Dumper=SafeDumper)
yaml.add_constructor('!Ratio', ratio_yaml_constructor, Loader=SafeLoader)

class Result:
    '''
    total_hits: FASTA match (path str) -> (Dict of fasta_index (int) -> Dict[str to either Ratio,List[Peptide]])
    '''
    def __init__(self, total_hits=None, coverages=None):
        self.total_hits = total_hits if total_hits is not None else {}
        self.coverages = coverages if coverages is not None else {}

    def __repr__(self):
        return f"{self.total_hits}"

    def set_fasta_results(self, fasta_name: str, counts: list[int], mod_counts: list[int]):
        if fasta_name in self.total_hits.keys():
            print(f"Already saved results for {fasta_name}! Won't overwrite.")
            return
        zipped = zip(mod_counts,counts)
        self.total_hits[fasta_name] = {}
        for i, (x,y) in enumerate(zipped):
            self.total_hits[fasta_name][i+1] = {'ratio':Ratio(x,y),'peptides':[]}

    def append_peptide_to_fasta_results(self, peptide: Peptide, fasta_match: str, psm_col_name: str):
        # Retrieve current counts (dict[list[Ratio]])
        fasta_hit = None
        if fasta_match in self.total_hits.keys():
            fasta_hit = self.total_hits[fasta_match]
        else:
            fasta_hit = {}
            for i in range(len(fasta_sequences[fasta_match])):
                fasta_hit[i+1] = {'ratio':Ratio(0,0),'peptides':[]}
        # Get index of first char of peptide in fasta
        index = peptide.fasta_match_indices[fasta_match] + 1
        # Add (peptide.counts[col]) to the mod counts
        for mod in peptide.mod_dict[fasta_match]:
             num = int(mod[1:])
             #  print(f"{mod} -> {num}")
             fasta_hit[num]['ratio'].mod_count += peptide.counts[psm_col_name]
        # Add (peptide.counts[col]) to all overlapping indices
        for i in range(len(peptide.sequence)):
            fasta_hit[index+i]['ratio'].count += peptide.counts[psm_col_name]
            fasta_hit[index+i]['peptides'].append(peptide)
        self.total_hits[fasta_match] = fasta_hit

    def calculate_coverages(self):
        # Goal: return a dict of dict:
        # 'fasta': str: The complete fasta
        # 'cover': str: fasta, with spaces where not found
        # 'percent': float: percent coverage
        # 'total_found': list[int]: total number of scans found
        # 'total_modified': list[int]: total number of scans which were found modified
        # 'ratios': dict[int,float]: modification ratios (the true goal of this whole script)
        self.coverages = {}
        for fasta_match in self.total_hits.keys():
            coverage = {}
            coverage['fasta'] = fasta_sequences[fasta_match]
            chars = len(coverage['fasta'])
            coverage['total_found'] = [0] * chars
            coverage['total_modified'] = [0] * chars
            coverage['ratios'] = dict.fromkeys(self.total_hits[fasta_match].keys(),None)
            coverage['cover'] = ""
            for i,c in enumerate(coverage['fasta']):
                # v Dicts v
                ratio = self.total_hits[fasta_match][i+1]['ratio']
                coverage['ratios'][i+1] = ratio.ratio()
                # v Lists v
                coverage['cover'] += (c if ratio.count > 0 else '.')
                coverage['total_found'][i] = ratio.count
                coverage['total_modified'][i] = ratio.mod_count
            coverage['percent'] = (chars - coverage['cover'].count('.')) / chars
            self.coverages[fasta_match] = coverage
        return self.coverages

    def print_coverage(self, fasta_file: str, output=""):
        try:
            coverage = self.coverages[fasta_file]
        except KeyError as ke:
            print("No coverage for '{fasta_file}'.")
            return

        match output:
            case 'total_found':
                print(coverage[output])
            case 'total_modified':
                print(coverage[output])
            case 'ratios':
                print(coverage[output])
            case _:
                # 'pretty' output
                for i,c in enumerate(coverage['cover']):
                    if c == '.':
                        print(f"[strike gray69]{coverage['fasta'][i]}[/]",end='')
                    else:
                        print(f"[bold blue]{c}[/]",end='')
                print("")

def result_yaml_representer(dumper, result: Result):
    return dumper.represent_mapping('!Result', {
        'total_hits': result.total_hits,
        'coverages': result.coverages
    })

def result_yaml_constructor(loader, node):
    values = loader.construct_mapping(node, deep=True)
    result = Result(values['total_hits'],values['coverages'])
    return result

yaml.add_representer(Result, result_yaml_representer, Dumper=SafeDumper)
yaml.add_constructor('!Result', result_yaml_constructor, Loader=SafeLoader)

# We can use a dict to dump to a yaml log

# Other Globals:
psm_col_names = None #dataset_names
psm_col_count = None #dataset_count
psm_cols = None #dataset_scan_counts
peptide_col = None #peptide_sequences
mod_col = None #target_modifications
peptides = None #peptides
column_to_fastas = None #dataset_fasta_lists
dataset_results = None #dataset_results
coverage_maps = None #coverage_maps

def init_log():
    log_dict = {}
    log_dict['step'] = 'init'
    log_dict['dataset_names'] = None
    log_dict['dataset_count'] = None
    log_dict['dataset_scan_counts'] = None
    log_dict['peptide_sequences'] = None
    log_dict['target_modifications'] = None
    log_dict['peptides'] = None
    log_dict['dataset_fasta_lists'] = None
    log_dict['fasta_sequences'] = None
    log_dict['dataset_results'] = None
    log_dict['coverage_maps'] = None
    Path(log_filename).touch()
    atomic_write_yaml(log_dict, log_filename)
    return log_dict

'''
Steps
    'get-scans': get_psm_columns,
    'get-peptides': get_peptide_column,
    'get-target-modifications': get_target_modifications_column,
    'compile-peptides': compile_peptides,
    'select-fasta-files': select_fasta_files,
    'read-fasta-sequences': read_fasta_sequences,
    'align-peptides': align_peptides,
    'compute-results': compute_results,
    'get-coverages': get_coverages,
    'get-modification-lists': get_modification_lists,
'''

def restore_from_log(log_dict):
    global next_step, psm_col_names, psm_col_counts, psm_cols, peptide_col
    global mod_col, peptides, column_to_fastas, fasta_sequences, dataset_results
    global coverage_maps
    try:
        next_step = log_dict['step']
        if next_step == 'init' or next_step == 'get-scans':
            return
        psm_col_names = log_dict['dataset_names']
        psm_col_count = log_dict['dataset_count']
        psm_cols = log_dict['dataset_scan_counts']
        if next_step == 'get-peptides':
            return
        peptide_col = log_dict['peptide_sequences']
        if next_step == 'get-target-modifications':
            return
        mod_col = log_dict['target_modifications']
        if next_step == 'compile-peptides':
            return
        peptides = log_dict['peptides']
        if next_step == 'select-fasta-files':
            return
        column_to_fastas = log_dict['dataset_fasta_lists']
        if next_step == 'read-fasta-sequences':
            return
        fasta_sequences = log_dict['fasta_sequences']
        if next_step == 'align-peptides' or next_step == 'compute-results':
            return
        dataset_results = log_dict['dataset_results']
        if next_step == 'get-coverages':
            return
        coverage_maps = log_dict['coverage_maps']
    except KeyError as ke:
        print(f"There was an issue reading the log file! Please delete and start again.")
        sys.exit(-1)

def load_log(filename: str):
    filepath = Path(filename)
    if filepath.exists() and filepath.is_file():
        with open(filename,'r') as f:
            return yaml.safe_load(f) or init_log()
    return init_log()

def atomic_write_yaml(data, filepath: str):
    dir_name = Path(filepath).absolute().parent.resolve()
    with tempfile.NamedTemporaryFile('w', delete=False, dir=dir_name, prefix='._', suffix='.tmp') as tmp_file:
        if not tmp_file:
            print("ERROR: No tmp_file!")
        yaml.safe_dump(data, tmp_file)
        tmp_name = tmp_file.name
    os.replace(tmp_name, filepath)  # Atomic on most OSes

next_step = 'get-scans'
log_filename = f"{filename}.log.yaml"
log_dict = load_log(log_filename)

def log_event(event_name: str, event_data):
    global log_dict
    log_dict[event_name] = event_data
    atomic_write_yaml(log_dict, log_filename)


class NonEmptyValidator(Validator):
    def validate(self, document):
        if not document.text.strip():
            raise ValidationError(
                message="Input cannot be empty",
                cursor_position=len(document.text)
            )
        if not document.text.isalpha():
            raise ValidationError(
                message="Input must only include letters from the English alphabet",
                cursor_position=len(document.text)
            )

def get_psm_columns():
    # Firstly, we need to know if the file given is a combined or single data set.
    # Let's ask the user:
    choice = q.select(
        f"In '{filename}', are there results for a single sample, or are results combined?",
        choices=[
            "Single Sample",
            "Combined",
            "I don't know (Abort)"
        ]
    ).ask()

    if not choice or "Abort" in choice:
        print("Aborting.")
        sys.exit(-1)

    single = None
    if "Single" in choice:
        single = True
    elif "Combined" in choice:
        single = False

    if single == None:
        print("Unexpected choice. Aborting.")
        sys.exit(-1)

    # PSM/Scan counts:
    #   "# PSMs"
    #   "#Scan.*"
    #
    # Do note, however, that a Combined file may only list the total PSMs.
    #
    # If this is a combined file, we expect multiple PSM/Scan columns:
    global psm_col_names, psm_col_count
    psm_col_names = df.columns[df.columns.str.contains("# PSMs",case=True)]
    psm_col_count = len(psm_col_names.tolist())
    if psm_col_count == 0:
        psm_col_names = df.columns[df.columns.str.contains("#Scan",case=True)]
        psm_col_count = len(psm_col_names.tolist())
        if psm_col_count == 0:
            print(
            "Could not find peptide hit counts in file!",
            "Expected '# PSMs' or '#Scan [sample]'"
            )
            sys.exit(-1)
    if psm_col_count == 1 and not single:
        print("Found multiple peptide hit columns despite 'single' dataset (is it combined?):")
        print(psm_col_names.tolist())
        sys.exit(-1)
    if psm_col_count > 1 and single:
        print("Only found a single peptide hit column despite having a combined dataset:")
        print('\t'+psm_col_names)
        print("If this is the total count of peptide hits, please retrieve the individual counts.")
        sys.exit(-1)

    psm_col_names = psm_col_names.tolist()
    # Now we know the number of columns and their names
    log_event('dataset_count', psm_col_count)
    log_event('dataset_names', psm_col_names)

    # We now have the relevant column name(s).
    # We will now work column by column,
    #  where each column is done row by row.
    # We do need to align each PSM count column
    #  to the peptide column now
    global psm_cols
    psm_cols = []
    for name in psm_col_names:
        psm_cols.append(df[name].tolist())
    log_event('dataset_scan_counts', psm_cols)

    return 'get-peptides'

def get_peptide_column():
    # Finds the "peptide" column:
    peptide_col_raw = None
    try:
        peptide_col_raw = df['Annotated Sequence']
    except KeyError:
        try:
            peptide_col_raw = df['Peptide']
        except KeyError:
            print("Could not find 'Peptide'/'Annotated Sequence' column! Please double check input.")

    # Let's go ahead and clean up the peptide column rows
    #
    # Peptide sequence
    #   "Annotated Sequence" -> [X].ABCDEF.[Y] where ABCDEF is the peptide found
    #   "Peptide" -> X.ABC(+##.##)DEF.Y where ABCDEF is the peptide found
    #
    # Sadly, I have just found an example of A.BCDE (with no ending `.`).
    #   Now more work is necessary. No list comprehension for me!

    global peptide_col
    peptide_col = []
    for i,row in enumerate(peptide_col_raw):
        if row[1] == row[-2] == '.':
            # safe to remove ends
            row = row[2:-2]
        else:
            # Otherwise, which end to remove?
            if row[1] == '.':
                row = row[2:]
            else:
                row = row[:-2]
        # This line removes any parenthesized data:
        row = re.sub(r'\([^)]*\)', '', row)
        if any(c not in set('ACDEFGHIKLMNPQRSTVWY') for c in row):
            print(f"Row[{i}]: {peptide_col_raw[i]} -> {row}")
            sys.exit(-1)
        peptide_col.append(row)
    #  for i in range(len(peptide_col)):
        #  print(f"{peptide_col_raw[i]} -> {peptide_col[i]}")
    print(f"Found {len(peptide_col)} peptide sequences.")
    log_event('peptide_sequences',peptide_col)

    return 'get-target-modifications'

def split_mods(mod_list_str: str):
    if ';' not in mod_list_str:
        return [mod_list_str]

    if '[' not in mod_list_str:
        return mod_list_str.split(';')

    # This case is more difficult:
    # NxAcetyl [X1(##); ...; Yn(##)]
    # Also:
    # 1xAcetyl [Kn(##)]

    # Start by grabbing all '[...]' strings:
    brackets = re.findall(r'\[.*\]',mod_list_str)
    # If there are no ';' in any of these, we're done
    if not any(';' in bracket for bracket in brackets):
        return mod_list_str.split(';')

    # NOW we have to do work.
    # For each bracket with colons, let's split it into n strings:
    # NxAcetyl [X1; Y2; Z3] -> 1xAcetyl [X1]; 1xAcetyl [Y2]; 1xAcetyl [Z3]
    mods = []
    mod_list_str_updated = mod_list_str
    while mod_list_str_updated != "":
        index = mod_list_str_updated.find(';')
        if index < 0:
            mods.append(mod_list_str_updated)
            mod_list_str_updated = ""
            continue
        bracket_index = mod_list_str_updated.find(']')
        if bracket_index < 0:
            # Something strange is occurring. I do not recognize this pattern.
            raise RuntimeError(f"Unexpected string {mod_list_str}")
        if bracket_index < index:
            # In other words, if the mod entry ends before the semicolon
            splat = mod_list_str_updated.split(';',maxsplit=1)
            mods.append(splat[0])
            mod_list_str_updated = splat[1]
            continue
        # Otherwise, we have [ ; ]
        match = re.search(r'\[.*\]',mod_list_str_updated)
        if not match:
            raise RuntimeError(f"Unexpected string {mod_list_str}")
        bracket = match.group()
        semi_count = bracket.count(';')
        inner_list = bracket.strip('[]').split(';')
        mod_list_before_bracket = mod_list_str_updated.split('[')[0]
        if mod_list_before_bracket[0] in '123456789':
            mod_list_before_bracket = '1'+mod_list_before_bracket[1:]
        for inner in inner_list:
            mods.append(f"{mod_list_before_bracket}[{inner}]")
        mod_list_str_updated = mod_list_str_updated.split(';',maxsplit=semi_count)[-1]
    return [mod.strip() for mod in mods]

def remove_uncertain_peptides(mod_col: list[list[str]], pep_col: list[str], psm_cols: list[list[str]]):
    if len(mod_col) != len(pep_col):
        raise RuntimeError("Modification and Peptide lists were not even!")
    for i in reversed(range(len(mod_col))):
        if any('/' in mod for mod in mod_col[i]):
            # this is overkill.
            # We could end up removing a certain peptide as well.
            # Let's notify the user of this important destruction of data
            print(f"Peptide [blue]{pep_col[i]}[/] contains mods {mod_col[i]}")
            print("Due to uncertainty, this peptide is being discarded!!")
            mod_col.pop(i)
            pep_col.pop(i)
            for col in psm_cols:
                col.pop(i)

def simplify_mods(mod_col: list[list[str]]):
    new_mod_col = []
    for i in range(len(mod_col)):
        new_mod_list = []
        for mod in mod_col[i]:
            match = re.search(r'[A-Z]\d+',mod)
            if not match:
                raise RuntimeError(f"Unexpected modification string '{mod}'!")
            new_mod_list.append(match.group())
        new_mod_col.append(new_mod_list)
    return new_mod_col

def get_target_modifications_column():
    # Next, get the modifications:
    #   "AScore"
    #   "Modifications"

    mod_col_raw = None
    try:
        mod_col_raw = df['AScore']
    except KeyError:
        try:
            mod_col_raw = df['Modifications']
        except KeyError:
            print("Could not find 'AScore'/'Modifications' column! Please double check input.")

    # Let's first find out which modification(s) we're looking for:
    choice = q.select(
        "Which type of modification are you looking for?",
        choices=[
            "Acetyl",
            "Other"
        ]
    ).ask()

    if choice is None:
        print("Expected to analyze a modification. Aborting.")
        sys.exit(-1)

    if choice == "Other":
        choice = q.text("Which modification?:",validate=NonEmptyValidator()).ask()

    global target_mod
    target_mod = choice
    print(f"Analyzing all '{target_mod}' modifications!")
    log_event('target_mod',target_mod)

    # "AScore"
    #    eg.:
    #    "C1:Carbamidomethylation:1000.00;M32:Oxidation (M):1000.00"
    #    "C5:Carbamidomethylation:1000.00;K14:Acetylation (K):1000.00"
    #    "T1:Acetylation (STY):99.62;M2:Oxidation (M):1000.00;K3:Acetylation (K):1000.00"
    # "Modifications"
    #    eg.:
    #    "1xCarbamidomethyl [C5]; 1xAcetyl [K14(100)]"
    #    "1xAcetyl [K5(100)]"
    #    "1xAcetyl [T/K/S]"
    # First split by ';', then find 'Acetyl' or whatever target modification it is.

    global mod_col
    mod_col = ['' if str(mod) == 'nan' else str(mod) for mod in mod_col_raw]
    mod_col = [split_mods(mod) for mod in mod_col]
    # mod_col is now a 2D array
    mod_col = [[mod for mod in mod_list if target_mod in mod] for mod_list in mod_col]
    # At this point, we have either:
    #   1: "Kn:Acetylation (K):####.##"
    #   2: "Sn:Acetylation (STY):##.##"
    #   3: "1","2"
    #   4: "1xAcetyl [Kn(##)]"
    #   5: "1xAcetyl [T/K/S]"
    # We cannot use case 5, so the entire peptide hit should be DISCARDED
    #  (not counted toward either acetylated/unacetylated).

    global peptide_col, psm_cols
    remove_uncertain_peptides(mod_col=mod_col,pep_col=peptide_col,psm_cols=psm_cols)

    # Notably, all other instances contain a residue-number "Xn" substring.
    # We only have one target_mod. We can simplify each mod to "Xn"
    mod_col = simplify_mods(mod_col)

    # Rewrite psm_cols, peptide_col, mod_col:
    log_event('dataset_scan_counts',psm_cols)
    log_event('peptide_sequences',peptide_col)
    log_event('target_modifications',mod_col)

    return 'compile-peptides'

def compile_peptides():
    # Now we can compile all of this into a nice list of Peptides:
    global peptides
    peptides = []
    for i in range(len(mod_col)):
        counts = {}
        for psm_col,psm_col_name in zip(psm_cols,psm_col_names):
            counts[psm_col_name] = psm_col[i]
        peptides.append(Peptide(peptide_col[i],mod_col[i],counts))
    
    log_event('peptides',peptides)

    #  for peptide in peptides:
        #  print(peptide.sequence)
    return 'select-fasta-files'

# Recall: The overall goal is to count the total "finds" of each residue as acetylated or not.
# This means we need the FASTA(s)

def flatten_nested_list(nested: list[list[str]]):
    return [item for sublist in nested for item in sublist]

def unique(lst: list):
    return list(set(lst))

def unique_flattened_nested_str_list(nested: list[list[str]]):
    return [val for val in set(flatten_nested_list(nested)) if val != ""]

def select_a_fasta(fasta_filenames_flat: list[str], col_name: str):
    choice = q.select(
        f"Select a FASTA file for the column labeled '{col_name}':",
        choices=fasta_filenames_flat+["Choose new FASTA file","Ignore this column"]
    ).ask()

    if "Ignore" in choice:
        return ""

    if "Choose new" in choice:
        return select_file(f"Select FASTA file for {col_name}",[("FASTA files",".FASTA .fasta .fas"),("All files",".*")])
    return choice

def select_col_fastas(cur_column_fastas: dict[str,list[str]], cur_col_name: str):
    fasta_filenames = []
    while True:
        cur_fasta_filename_lists = list(cur_column_fastas.values())
        cur_fasta_filename_lists.append(fasta_filenames)
        new_fasta_filename = select_a_fasta(unique_flattened_nested_str_list(cur_fasta_filename_lists),cur_col_name)
        fasta_filenames.append(new_fasta_filename)
        if new_fasta_filename == "":
            break

        # There are cases where a second FASTA is necessary (due to combination experiments)
        choice = q.select(
            f"Do you need to select a second FASTA for '{cur_col_name}'?",
            choices=["Yes", "No"]
        ).ask()
        
        if choice == "No":
            break
    return fasta_filenames

def select_all_fastas(psm_col_names: list[str]):
    column_to_fastas = {}
    for name in psm_col_names:
        column_to_fastas[name] = select_col_fastas(column_to_fastas,name)

    return column_to_fastas

def print_pathname(path: str):
    p = Path(path)
    sep = os.path.sep
    p1 = f"{p.parent}{os.path.sep}"
    p2 = f"{p.name}"
    if path == "":
        print("[red]Ignored[/]")
    else:
        print("[cyan]"+p1+"[/]"+"[bold green]"+p2+"[/]")

def print_column_fastas(column_fastas: dict[str, list[str]]):
    for k,v in column_fastas.items():
        print(f"[bold blue]{k}[/]:")
        for path in v:
            print('\t',end='')
            print_pathname(path)

def select_fasta_files():
    global column_to_fastas
    if not column_to_fastas:
        column_to_fastas = select_all_fastas(psm_col_names)
        log_event('dataset_fasta_lists',column_to_fastas)

    while True:
        print("[bold yellow]Please review the following peptide data sets for accuracy:[/]")
        print_column_fastas(column_to_fastas)

        choice = q.select(
            "Are the listed FASTAs accurate?",
            choices=[
                "Yes, continue to analysis",
                "No, I'd like to make a change",
                "Forget it all!"
            ]
        ).ask()

        if "Forget" in choice:
            sys.exit(0)
        if "Yes" in choice:
            break

        # Change necessary
        choice = q.select(
            "Restart FASTA selection or pick a file to change?",
            choices=[
                "Restart",
                "Replace a FASTA file (all occasions of file)",
                "Edit a single column"
            ]
        ).ask()

        if choice == "Restart":
            column_to_fastas = select_all_fastas(psm_col_names)
            log_event('dataset_fasta_lists',column_to_fastas)
            continue

        if "all occasions" in choice:
            # Pick a FASTA file to replace:
            choice = q.select(
                "Which FASTA to replace?",
                choices=unique_flattened_nested_str_list(list(column_to_fastas.values()))+["Nevermind, we're all good"]
            ).ask()

            if "Nevermind" in choice:
                # Assuming no FASTA will ever have "Nevermind" as part of the name...
                continue

            # Need to update all instances of the given FASTA file.
            #  indices = [(i,j) for j, fasta in enumerate(fasta_sublist) for i, fasta_sublist in enumerate(fasta_filenames) if fasta == choice]
            indices = [(name, i) for name, fasta_list in column_to_fastas.items()
                                for i, fasta in enumerate(fasta_list)
                                if fasta == choice]

            new_fasta_filename = select_file()
            for k,i in indices:
                column_to_fastas[k][i] = new_fasta_filename
            log_event('dataset_fasta_lists',column_to_fastas)
        else:
            # Pick a column name to edit
            choice = q.select(
                "Which column to edit?",
                choices=list(column_to_fastas.keys())+["Nevermind, we're all good"]
            ).ask()

            if "Nevermind" in choice:
                continue

            column_to_fastas[choice] = select_col_fastas(column_to_fastas,choice)
            log_event('dataset_fasta_lists',column_to_fastas)
    
    log_event('dataset_fasta_lists',column_to_fastas)

    return 'read-fasta-sequences'


# Continuing to analysis!

def read_fasta_sequences():
    # Starting by exporting FASTA files to Dictionary:
    global fasta_sequences
    fasta_sequences = {}
    for fasta_filename in unique_flattened_nested_str_list(list(column_to_fastas.values())):
        with open(fasta_filename,'r') as fasta_file:
            fasta_lines = fasta_file.readlines()
            for i in reversed(range(len(fasta_lines))):
                if fasta_lines[i][0] == '>':
                    fasta_lines.pop(i)
            fasta_sequences[fasta_filename] = "".join(fasta_lines)
            fasta_file.close()

    print("Here are the FASTAs you selected:")
    for fasta,seq in fasta_sequences.items():
        print(f"\t'{fasta}':")
        print(f"\t\t{seq}")

    log_event('fasta_sequences',fasta_sequences)
    return 'align-peptides'


# Next steps:
#   Start by aligning peptides to FASTAs:
#   Go peptide by peptide, finding which FASTA(s?) they are a part of
#   Renumber modifications

def align_peptides():
    global peptides
    for n,peptide in enumerate(peptides):
        if not peptide.sequence or peptide.sequence == "":
            print(f"Peptide {n} was empty! Skipping.")
            continue
        for fasta,seq in fasta_sequences.items():
            index = seq.find(peptide.sequence)
            if index < 0:
                continue
            # Otherwise, this sequence contains the peptide.
            peptide.fasta_match_indices[fasta] = index
            peptide.mod_dict[fasta] = []
            for i,mod in enumerate(peptide.mod_list):
                num = int(mod[1:])
                peptide.mod_dict[fasta].append(f"{mod[0]}{num+index}")
        if len(peptide.mod_dict.keys()) > 1:
            print(f"Peptide {n} ({peptide.sequence}) was found in multiple FASTAs:")
            print(f"\t{peptide.fasta_match_indices}")

    log_event('peptides',peptides)
    return 'compute-results'

# For a given sample, I want:
#  * FASTA matches
#  * Coverage (method)
#  * PSM count sums
#  * Acetylation count sums
#  * Acetylation ratios (method)

def compute_results():
    # Go column by column:
    #   Go peptide by peptide:
    #       Increment column counters
    global dataset_results
    dataset_results = {}
    for col_name,fasta_list in column_to_fastas.items():
        if not fasta_list or fasta_list == [""]:
            continue
        dataset_results[col_name] = Result()
        for peptide in peptides:
            fasta_matches = set(fasta_list).intersection(set(peptide.mod_dict.keys()))
            #  if len(fasta_matches) > 1:
                #  print("How did we get here?")
                #  print(f"\tcol_name: {col_name}")
                #  print(f"\tfasta_list: {fasta_list}")
                #  print(f"\tpeptide.mod_dict.keys(): {peptide.mod_dict.keys()}")
            #  print(list(fasta_matches))
            if not fasta_matches:
                # This peptide is unrelated to this column.
                continue # to the next peptide
            # This peptide is part of at least one of the given fastas
            #  print(peptide)
            for fasta_match in fasta_matches:
                dataset_results[col_name].append_peptide_to_fasta_results(peptide,fasta_match,col_name)

    log_event('dataset_results',dataset_results)
    return 'get-coverages'

def get_coverages():
    global coverage_maps
    for dataset,result in dataset_results.items():
        print(f"[bold yellow]{dataset} coverages:[/]")
        result.calculate_coverages()
        for fasta_file in result.total_hits.keys():
            print(f"\t'{fasta_file}' ({result.coverages[fasta_file]['percent']*100.0:.2f}%):")
            result.print_coverage(fasta_file)
    # This data can be put into a file, if desired:
    choice = q.select(
        "Do you wish to save these coverages to a file?",
        choices = ["Yes", "No"]
    ).ask()

    if choice == "Yes":
        while True:
            coverage_filename = filedialog.asksaveasfilename(title="Choose a new file for the coverage maps")
            if not coverage_filename:
                break
            coverage_filepath = Path(coverage_filename)
            with open(coverage_filepath,'w') as cfile:
                for dataset,result in dataset_results.items():
                    cfile.write(f"Coverages for {dataset}:\n")
                    for fasta_file in result.total_hits.keys():
                        cfile.write(f"\t'{fasta_file}' ({result.coverages[fasta_file]['percent']*100.0:.2f}%):\n")
                        cfile.write(f"\t{result.coverages[fasta_file]['cover']}\n")
                    cfile.write("\n")
            coverage_maps = coverage_filepath
            break
    return 'get-modification-lists'

def get_modification_lists():
    print("The modification ratios are ready to be saved.")
    choice = q.select(
        f"You have {len(dataset_results)} dataset(s) ready to be saved. Will you save them now?",
        choices = ["Yes", "No"]
    ).ask()
    if choice == "Yes":
        while True:
            filename = filedialog.asksaveasfilename(title="Save as .zip",defaultextension=".zip")
            if not filename:
                print("Results not yet saved. Run this program again to save them!")
                break
            filepath = Path(filename)
            with tempfile.TemporaryDirectory() as temp_dir:
                temp_path = Path(temp_dir)
                temp_filepaths = []
                for dataset,result in dataset_results.items():
                    temp_filename_base = f"{dataset}"
                    for fasta_file in result.total_hits.keys():
                        temp_filename = f"{temp_filename_base}_{Path(fasta_file).stem}.mods"
                        temp_filepath = temp_path / temp_filename # nifty
                        with open(temp_filepath, 'w') as temp_file:
                            temp_file.write("Residue#\tResidue\tModified\tTotal\tRatio\n")
                            for i in sorted(result.total_hits[fasta_file]):
                                ratio = result.total_hits[fasta_file][i]['ratio']
                                temp_file.write(f"{i}\t{fasta_sequences[fasta_file][i-1]}\t")
                                temp_file.write(f"{ratio.mod_count}\t{ratio.count}\t{ratio.ratio()}\n")
                        temp_filepaths.append(temp_filepath)
                if len(temp_filepaths) < 1:
                    raise RuntimeError("There were no files to save! You might want to contact the developer.")
                elif len(temp_filepaths) == 1:
                    # Just save as .mods
                    shutil.copy(temp_filepaths[0], filepath.with_suffix('.mods'))
                else:
                    with zipfile.ZipFile(filepath.with_suffix('.zip'), 'w', zipfile.ZIP_DEFLATED) as zipf:
                        for temp_filepath in temp_filepaths:
                            zipf.write(temp_filepath, arcname=temp_filepath.name)
            break

    return 'Done!'

# Steps:
step_functions = {
    'get-scans': get_psm_columns,
    'get-peptides': get_peptide_column,
    'get-target-modifications': get_target_modifications_column,
    'compile-peptides': compile_peptides,
    'select-fasta-files': select_fasta_files,
    'read-fasta-sequences': read_fasta_sequences,
    'align-peptides': align_peptides,
    'compute-results': compute_results,
    'get-coverages': get_coverages,
    'get-modification-lists': get_modification_lists,
}

# Now we can resume execution to where we left it, if desired:
if log_dict['step'] in step_functions:
    choice = q.select(
        f"It looks like you've analyzed '{Path(filename).name}' before. Do you want to resume analysis?",
        choices = [
            f"Resume analysis of {Path(filename).name} (will resume at '{log_dict['step']}')",
            "Resume at a previous step",
            "Restart analysis"
        ]
    ).ask()

    if "Resume" in choice:
        restore_from_log(log_dict)
        if "previous step" in choice:
            prev_steps = ['init']
            for step in step_functions.keys():
                prev_steps.append(step)
                if step == log_dict['step']:
                    break
            choice = q.select(
                "Which step to return to?",
                choices = prev_steps
            ).ask()
            next_step = choice
    elif "Restart" in choice:
        log_dict = init_log()
elif log_dict['step'] == 'init':
    next_step = 'get-scans'

while next_step in step_functions.keys():
    log_event('step', next_step)
    func = step_functions[next_step]
    next_step = func()

messagebox.showinfo("Success!",f"{next_step}")

