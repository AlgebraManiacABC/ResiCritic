import math
from pathlib import Path
from peptide import Peptide
from ratio import Ratio
import tui
from exit_codes import EXIT_SUCCESS, EXIT_FAILURE

class Result:
    '''
    Comprehensive Dictionary of Experimental Data (for a single FASTA and single experiment)
    complete: Whether more peptides can be added (False) or not (True)
    residues: main dictionary - sequence index keys, dict values:
        'id'    -> 1-Letter code as char
        'ratio' -> Ratio(modified,total)
        'cover' -> 'id' if 'ratio' != 'nan' else '.'
        'peptides' -> List of Peptides which overlap at this residue
    '''
    def __init__(self, experiment_name: str, fasta_id: str, sequence: str, start_num: int=1):
        if not isinstance(sequence,str):
            raise TypeError(f"Expected a string sequence; got {type(sequence)}.")
        self.experiment_name = experiment_name
        self.fasta_id = fasta_id
        self.fasta_sequence = sequence
        self.residues = {}
        self.complete = False
        self.coverage_percent = None
        for i,resi in enumerate(sequence):
            self.residues[i+start_num] = dict.fromkeys(['id','ratio','cover','peptides'])
            self.residues[i+start_num]['id'] = resi
            self.residues[i+start_num]['ratio'] = Ratio(0,0)
            self.residues[i+start_num]['peptides'] = []

    def __str__(self):
        if not self.complete:
            return "Incomplete (use Result.complete_results())"
        return "TODO"

    def add_peptide_to_result(self, peptide: Peptide) -> int:
        '''
        For the given experiment and FASTA, add the peptide's
            counts to each overlapping residue. Also considers
            whether the residue was modified.
        '''
        if self.complete:
            tui.say_error("Results are already complete; cannot add new peptide!")
            return EXIT_FAILURE
        try:
            # Index of peptide's first residue in whole sequence
            zero_index = self.fasta_sequence.find(peptide.sequence)
            if zero_index < 0:
                return EXIT_SUCCESS

            for modification in peptide.mod_list:
                resi_num = int(modification[1:]) + zero_index
                self.residues[resi_num]['ratio'].numer += peptide.counts[self.experiment_name]

            for i in range(len(peptide.sequence)):
                self.residues[zero_index+1+i]['ratio'].denom += peptide.counts[self.experiment_name]
                self.residues[zero_index+1+i]['peptides'].append((zero_index+1,peptide))
        except (AttributeError, KeyError, IndexError, ValueError) as e:
            tui.say_error(f"Unexpected error while adding peptide to result! {e}")
            return EXIT_FAILURE
        return EXIT_SUCCESS

    def complete_results(self) -> None:
        '''
        Completes the Result by calculating the coverage map
        '''
        if self.complete:
            return
        self.coverage_percent = Ratio(0,len(self.residues))
        for _,resi in self.residues.items():
            if resi['ratio'].denom > 0:
                self.coverage_percent.numer += 1
                resi['cover'] = resi['id']
            else:
                resi['cover'] = '.'
        self.complete = True
        return

    def get_coverage_aligned(self, width: int=80) -> str:
        if not self.complete:
            return "---"
        # Using the following "look":
        # FASTA 1 ABCDEFGHIJ (width)
        # cover   ABC.EFGH.J
        coverage_string = ''.join([self.residues[i]['cover'] for i in sorted(self.residues)])

        # How many pairs of lines?
        line_pairs = math.ceil(len(coverage_string) / width)

        aligned_string = ""
        gap_width = len(str(len(coverage_string)))
        gap = ''.join([' ']*gap_width)
        for i in range(line_pairs):
            first_index = (i * width) + 1
            last_index = min(first_index + width - 1, len(coverage_string))
            aligned_string += f"{'FASTA':>{len(self.experiment_name)}s} {first_index:{gap_width}d} {self.fasta_sequence[first_index-1:first_index+width]} {last_index}\n"
            aligned_string += f"{self.experiment_name} {gap} {coverage_string[first_index-1:first_index+width]} {last_index}\n\n"
        return aligned_string

def compute_results(context: dict):
    tui.say("Computing results...")
    all_results = dict.fromkeys(context['experiment_names'])
    for experiment_name, fasta_file_list in context['experiment_fasta_files'].items():
        if not fasta_file_list or fasta_file_list == ['']:
            tui.say(f"Skipping {experiment_name}.")
            continue
        tui.say(f"Analyzing experiment [bold blue]{experiment_name}[/] with fastas:")
        all_results[experiment_name] = {}
        for fasta_filename in fasta_file_list:
            tui.say(f"\t- '{Path(fasta_filename).name}'")
            fasta_seq = context['fasta_sequences'][fasta_filename]
            if isinstance(fasta_seq,list):
                # Multi-sequence FASTA
                for i,fasta_sequence in enumerate(fasta_seq):
                    fasta_id = f"{fasta_filename}_{i}"
                    all_results[experiment_name][fasta_id] = Result(experiment_name,fasta_id,fasta_sequence)
                    for peptide in context['peptides']:
                        err = all_results[experiment_name][fasta_id].add_peptide_to_result(peptide)
                        if err:
                            return EXIT_FAILURE
                    all_results[experiment_name][fasta_id].complete_results()
            else:
                # Single-sequence FASTA
                all_results[experiment_name][fasta_filename] = Result(experiment_name,fasta_filename,fasta_seq)
                for peptide in context['peptides']:
                    err = all_results[experiment_name][fasta_filename].add_peptide_to_result(peptide)
                    if err:
                        return EXIT_FAILURE
                all_results[experiment_name][fasta_filename].complete_results()
    context['results'] = all_results
    context['step'] = 'export-coverages'
    return EXIT_SUCCESS
