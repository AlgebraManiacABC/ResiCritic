import re
from pathlib import Path
import tui
from peptide import Peptide
import util
import files
from exit_codes import EXIT_SUCCESS, EXIT_FAILURE

def read_experiments(context: dict):
    experiment_names = read_experiment_names(context)
    if not experiment_names:
        return EXIT_FAILURE

    if len(experiment_names) == 1:
        confirm = tui.confirm(f"Only one dataset was found in {Path(context['filename']).name}. Do you wish to continue?")
        if not confirm:
            tui.say_error("Aborting analysis. Please check your file for columns labeled '# PSMs' or '#Scan [sample]'.")
            return EXIT_FAILURE

    # At this point, we can be reasonably sure those pieces of data are acceptable
    context['experiment_count'] = len(experiment_names)
    context['experiment_names'] = experiment_names

    peptide_sequences = read_peptide_sequences(context)
    if not peptide_sequences:
        return EXIT_FAILURE

    modification_lists_raw = read_modification_lists(context)

    experiment_columns = []
    try:
        for name in experiment_names:
            experiment_columns.append(context['dataframe'][name].tolist())
    except KeyError:
        tui.say_error("Unexpected error getting experimental data!")
        return EXIT_FAILURE

    if context['experiment_count'] == 1:
        experiment_names = [Path(context['filename']).stem]
        context['experiment_names'] = experiment_names

    target_mod = tui.choose_between("Which type of modification are you looking for?",
                                    ["Acetyl","Other"])
    if target_mod == "Other":
        target_mod = tui.answer_alpha_only("Which modification?")

    if not target_mod:
        tui.say_error(f"Expected a string modification, but received {target_mod}! Aborting.")
        return EXIT_FAILURE

    tui.say(f"Analyzing all '{target_mod}' modifications!")
    context['target_modification'] = target_mod

    peptide_list = []
    for i,peptide_sequence in enumerate(peptide_sequences):
        try:
            cleaned_modification_list = clean_modification_list(modification_lists_raw[i],target_mod)
        except RuntimeError as e:
            tui.say(f"Peptide [blue]{peptide_sequence}[/]:")
            tui.say_error(e)
            return EXIT_FAILURE
        if cleaned_modification_list is None:
            tui.say(f"Peptide [blue]{peptide_sequence}[/] contains modification string '{modification_lists_raw[i]}'.")
            tui.say_error("Due to uncertainty, this peptide is being discarded!!")
            continue
        scan_counts = {}
        for j,name in enumerate(experiment_names):
            scan_counts[name] = experiment_columns[j][i]
        peptide_list.append(Peptide(peptide_sequence,cleaned_modification_list,scan_counts))

    context['peptides'] = peptide_list
    context['step'] = 'select-fasta-files'

    return EXIT_SUCCESS

def read_experiment_names(context: dict):
    experiment_names = context['dataframe'].columns[context['dataframe'].columns.str.contains("# PSMs",case=True)].tolist()
    if len(experiment_names) == 0:
        experiment_names = context['dataframe'].columns[context['dataframe'].columns.str.contains("#Scan",case=True)].tolist()
    if len(experiment_names) == 0:
        tui.say_error(f"Could not find peptide hit counts in {Path(context['filename']).name}!" +
                      "Expected '# PSMs' or '#Scan [sample]'")
        return None
    return experiment_names

def read_peptide_sequences(context: dict) -> list[str]:
    try:
        peptide_list_raw = context['dataframe']['Annotated Sequence']
    except KeyError:
        try:
            peptide_list_raw = context['dataframe']['Peptide']
        except KeyError:
            tui.say_error("Could not find 'Peptide' or 'Annotated Sequence' column! Please double check input.")
            return None

    peptide_sequences = []
    for row in peptide_list_raw:
        cleaned_sequence = clean_peptide_sequence(row)
        if not cleaned_sequence:
            return None
        peptide_sequences.append(cleaned_sequence)
    return peptide_sequences

def read_modification_lists(context: dict):
    try:
        modification_lists = context['dataframe']['AScore']
    except KeyError:
        try:
            modification_lists = context['dataframe']['Modifications']
        except KeyError:
            tui.say_error("Could not find 'AScore' or 'Modifications' column in '{Path(context['filename']).name}'!")
            return []
    return modification_lists

def clean_peptide_sequence(dirty: str) -> str:
    # Remove brackets
    dirty = re.sub(r'[\[\]]', '', dirty)
    # Remove trailing ends, if there
    if dirty[1] == '.':
        dirty = dirty[2:]
    if dirty[-2] == '.':
        dirty = dirty[:-2]
    # Remove parentheticals
    clean = re.sub(r'\([^)]*\)', '', dirty)
    # Check that we only have accepted residues
    if any(c not in set('ACDEFGHIKLMNPQRSTVWY') for c in clean):
        tui.say_error(f"An unexpected character was found in the 'cleaned' peptide sequence '{clean}'.")
        tui.say("Firstly, ensure all residues are capitalized. If the error persists and the sequence looks correct, contact the developer.")
        return ""

    return clean

def modification_uncertain(mod: str):
    return '/' in mod

def clean_modification_list(dirty, target_mod: str):
    dirty = '' if str(dirty) == 'nan' else str(dirty)
    dirty_split = split_modification_str_to_list(dirty)
    clean_list = [mod for mod in dirty_split if target_mod in mod]
    # At this point, we have either:
    #   1: "Kn:Acetylation (K):####.##"
    #   2: "Sn:Acetylation (STY):##.##"
    #   3: "1xAcetyl [Kn(##)]"
    #   4: "1xAcetyl [T/K/S]"
    # We cannot use case 4, so the entire peptide hit should be DISCARDED
    #  (not counted toward either acetylated/unacetylated).
    if any(modification_uncertain(mod) for mod in clean_list):
        return None

    for i,mod in enumerate(clean_list):
        # We can further simplify our list by
        #  grabbing only the "Xn" string (e.g., K107)
        match = re.search(r'[A-Z]\d+',mod)
        if not match:
            raise RuntimeError(f"Unexpected modification string '{dirty}'!")
        clean_list[i] = match.group()

    return clean_list

def split_modification_str_to_list(mods: str) -> list:
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
    if ';' not in mods:
        # Only a single modification.
        # Convert to list
        return [mods]
    if '[' not in mods:
        # "AScore" version. Split by semicolon
        return mods.split(';')

    inside_brackets = re.findall(r'\[.*\]',mods)
    if not any(';' in bracket_text for bracket_text in inside_brackets):
        # No semicolon inside brackets; simply split by semicolon
        return mods.split(';')

    # All that remains now are the following cases:
    # NxAcetyl [X1; Y2; Z3]
    # Go semicolon by semicolon
    mods_split = []
    mods_updated = mods
    while mods_updated != "":
        semicolon_index = mods_updated.find(';')
        if semicolon_index < 0:
            # No more semicolons
            mods_split.append(mods_updated)
            mods_updated = ""
            continue
        right_bracket_index = mods_updated.find(']')
        if right_bracket_index < 0:
            # Unrecognized modification string
            raise RuntimeError(f"Unexpected modification string '{mods}' - please consult developer")
        if right_bracket_index < semicolon_index:
            # The modification entry ends before the next semicolon ("] ;")
            splat = mods_updated.split(';',maxsplit=1)
            mods_split.append(splat[0])
            mods_updated = splat[1]
            continue
        # Otherwise we have [ ; ]
        match = re.search(r'\[.*\]', mods_updated)
        if not match:
            # Unrecognized modification string
            raise RuntimeError(f"Unexpected modification string '{mods}' - please consult developer")
        bracket_string = match.group()
        semi_count = bracket_string.count(';')
        inner_list = bracket_string.strip('[]').split(';')
        mod_prefix = mods_updated.split('[',maxsplit=1)[0]
        if mod_prefix[0] in '123456789':
            # #xTargetmod
            # Since we are splitting into one each:
            # 1xTargetmod
            mod_prefix = '1'+mod_prefix[1:]
        for inner in inner_list:
            mods_split.append(f"{mod_prefix}[{inner}]")
        mods_updated = mods_updated.split(';',maxsplit=semi_count)[-1]
    return [mod.strip() for mod in mods_split]

def remove_fasta_comments(fasta_lines: list[str]) -> list:
    for i in reversed(range(len(fasta_lines))):
        if fasta_lines[i][0] == '>':
            fasta_lines.pop(i)
    return fasta_lines

def read_multi_fasta_sequence(fasta_filename: str, fasta_lines: list[str]) -> list[str]:
    desc_indices = [i for i,line in enumerate(fasta_lines) if line[0] == '>']
    multi_sequence = []
    current_desc = ""
    current_fasta = ""
    fasta_count = 0
    for i,line in fasta_lines:
        if i in desc_indices:
            # new fasta start
            if current_fasta:
                # end current fasta
                fasta_count += 1
                multi_sequence.append(current_fasta)
                current_fasta = ""
                current_desc = ""
            if not current_desc:
                current_desc = line[1:].strip()
            else:
                current_desc += ' ' + line[1:].strip()
        else:
            # continue on current fasta
            current_fasta += line.strip()
    # After the loop:
    if current_desc and not current_fasta:
        # Trailing descriptor/comment.
        # Warn user, but ignore
        tui.say(f"Warning: '{fasta_filename}' ends with a '>' line (are you missing a sequence?)")
    elif current_desc and current_fasta:
        fasta_count += 1
        multi_sequence.append(current_fasta)

    return multi_sequence

def read_fasta_sequences(context: dict):
    fasta_sequences = {}
    for fasta_filename in util.unique_flattened_nested_str_list(list(context['experiment_fasta_files'].values())):
        fasta_lines = files.read_text(fasta_filename).splitlines()
        # If there are more than one lines beginning with '>' with text between and after each,
        #   then we treat this fasta file as a "multi-fasta file".
        #   Double check with user first, of course.
        desc_indices = [i for i,line in enumerate(fasta_lines) if line[0] == '>']
        if len(desc_indices) <= 1:
            # Only one (or no) descriptor/comment line
            fasta_sequences[fasta_filename] = ''.join(remove_fasta_comments(fasta_lines))
            continue
        # Otherwise, there are multiple lines with '>'.
        multiple = tui.confirm(f"It looks like {Path(fasta_filename).name} contains multiple sequences. Does it?")
        if not multiple:
            fasta_sequences[fasta_filename] = ''.join(remove_fasta_comments(fasta_lines))
            continue
        # Otherwise, we need to treat the fasta as multiple...
        # This also means we need to place them all into the dictionary value as a list
        fasta_sequences[fasta_filename] = read_multi_fasta_sequence(fasta_filename,fasta_lines)

    context['fasta_sequences'] = fasta_sequences
    context['step'] = 'compute-results'

    return EXIT_SUCCESS
