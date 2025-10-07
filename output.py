import tempfile
import zipfile
import shutil
from pathlib import Path
import tui
import files
from result import Result
from exit_codes import EXIT_SUCCESS, EXIT_FAILURE

def export_text(export_what: str, temp_filepath: Path, result: Result):
    if export_what == 'coverage':
        temp_filepath.write_text(result.get_coverage_aligned())
        return

    if export_what == 'modifications':
        with open(temp_filepath, 'w') as temp_file:
            temp_file.write("Residue#\tResidue\tModified\tTotal\tRatio\n")
            for i in sorted(result.residues):
                ratio = result.residues[i]['ratio']
                temp_file.write(f"{i}\t{result.residues[i]['id']}\t")
                temp_file.write(f"{ratio.numer}\t{ratio.denom}\t{ratio.ratio()}\n")

def export_to_file(context: dict, export_what: str, filename: str):
    filepath = Path(filename)
    file_ext = ''
    match export_what:
        case 'coverage':
            file_ext = '.txt'
        case 'modifications':
            file_ext = '.tsv'
        case _:
            file_ext = '.txt'
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        temp_filepaths = []
        for experiment_name, results in context['results'].items():
            if not results:
                continue
            temp_filename_base = f"{experiment_name}"
            for fasta_id, result in results.items():
                temp_filename = f"{temp_filename_base}_{Path(fasta_id).stem}{file_ext}"
                temp_filepath = temp_path / temp_filename # Nifty Path operation
                export_text(export_what,temp_filepath,result)
                temp_filepaths.append(temp_filepath)
        if len(temp_filepaths) < 1:
            tui.say_error("There were no files to save! You may want to contact the developer.")
            return
        if len(temp_filepaths) == 1:
            # Just save as file_ext
            shutil.copy(temp_filepaths[0], filepath.with_suffix(file_ext))
            return
        # Otherwise, multiple files
        with zipfile.ZipFile(filepath.with_suffix('.zip'), 'w', zipfile.ZIP_DEFLATED) as zipf:
            for temp_filepath in temp_filepaths:
                zipf.write(temp_filepath, arcname=temp_filepath.name)

def export_coverages(context: dict):
    save_coverages = tui.confirm("Do you wish to save the coverages to file?")
    if save_coverages:
        coverage_filename = files.open_new_file("Choose a new file for the coverage maps")
        if not coverage_filename:
            tui.say("Skipping coverage map export - you may rerun the script to save your coverages!")
            return EXIT_SUCCESS
        export_to_file(context,'coverage',coverage_filename)
    context['step'] = 'export-modification-lists'
    return EXIT_SUCCESS

def export_modification_lists(context: dict):
    num_results = len([result.values() for result in context['results'].values() if result])
    save_modifications = tui.confirm(f"Do you wish to save the modification lists to a file? There are {num_results} experiment(s) to save.")
    if save_modifications:
        modification_list_filename = files.open_new_file("Choose a new file for the modification list(s)")
        if not modification_list_filename:
            tui.say("Skipping modification list export - you may rerun the script to save your data!")
            return EXIT_SUCCESS
        export_to_file(context,'modifications',modification_list_filename)
    context['step'] = 'done'
    return EXIT_SUCCESS
