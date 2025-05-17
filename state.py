from parser import read_experiments, read_fasta_sequences
from result import compute_results
from output import export_coverages, export_modification_lists

def get_fasta_files(context: dict):
    from tui import select_fasta_files
    return select_fasta_files(context)

# Steps:
state_functions = {
    'read-experiments': read_experiments,
    'select-fasta-files': get_fasta_files,
    'read-fasta-sequences': read_fasta_sequences,
    'compute-results': compute_results,
    'export-coverages': export_coverages,
    'export-modification-lists': export_modification_lists,
}

def run(context: dict):
    func = state_functions[context['step']]
    return func(context)
