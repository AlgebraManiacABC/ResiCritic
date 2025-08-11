from pathlib import Path
from tkinter import messagebox
from rich import print
import questionary as q
from questionary import Validator, ValidationError
import util
import files
import state
from exit_codes import EXIT_SUCCESS, EXIT_FAILURE
from version import __version__

class AlphaNonEmptyValidator(Validator):
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

def say(that: str = ""):
    print(that)

def say_error(that: str = ""):
    print(f"[bold red]{that}")

def print_welcome():
    say(f"Welcome to [bold green]ResiCritic![/] (Version {__version__})")
    say("Developed by [bold purple]Luke Dunn[/] at the [bold orange3]University of Missouri[/] (dunnl@umsystem.edu)")

def done():
    say("All done!")
    messagebox.showinfo("Success!","All done!")

def confirm(to_confirm: str):
    r = q.select(to_confirm, ["Yes","No"]).ask()
    if r == "Yes":
        return True
    return False

def answer_alpha_only(question: str) -> str:
    return q.text(question,validate=AlphaNonEmptyValidator()).ask()

def choose_between(ask: str="Please choose an option:",choices=None):
    r = q.select(ask,choices=choices).ask()
    return r

def already_analyzed_pick_step(filename: str, context: dict) -> str:
    choice = q.select(
        f"It looks like you've analyzed '{Path(filename).name}' before. Do you want to resume analysis?",
        choices = [
            f"Resume analysis of '{Path(filename).name}' (will resume at '{context['step']}')",
            "Resume at a previous step (overwrite any progress since the given step)",
            "Restart analysis completely (overwriting all prior progress - previous saved files are unaffected)"
        ]
    ).ask()

    if "Restart" in choice:
        return ""

    if "previous step" in choice:
        prev_steps = ['init']
        for step in state.state_functions:
            prev_steps.append(step)
            if step == context['step']:
                break
        prev_step = choose_between("Which step to return to?",prev_steps)
        return prev_step

    return context['step']

def select_single_fasta_file(fasta_filenames_flat: list[str], experiment_name: str):
    fasta_filenames_flat.sort()
    choice = q.select(
        f"Select a FASTA file for the experiment labeled '{experiment_name}':",
        choices=fasta_filenames_flat+["Choose new FASTA file","Ignore this column"]
    ).ask()

    if "Ignore" in choice:
        return ""

    if "Choose new" in choice:
        return files.open_file(f"Select FASTA file for {experiment_name}",[("FASTA files",".FASTA .fasta .fas"),("All files",".*")])
    return choice

def select_experiment_fastas(experiment_fasta_filenames: dict[str,list[str]], experiment_name: str):
    fasta_filenames = []
    while True:
        all_current_fastas_list = fasta_filenames.copy()
        for v in experiment_fasta_filenames.values():
            if v:
                all_current_fastas_list += v
        all_current_fastas = set(all_current_fastas_list)
        new_fasta_filename = select_single_fasta_file(util.remove_empty(all_current_fastas), experiment_name)
        fasta_filenames.append(new_fasta_filename)
        if new_fasta_filename == "":
            # Ignoring this experiment, or done selecting fastas
            break

        # Some experiments require multiple protein sequences
        # Allow the user to select multiple FASTAs if they don't have a combined one
        choice = q.select(
            f"Do you need to select a second FASTA for '{experiment_name}'?",
            choices=["Yes", "No"]
        ).ask()

        if choice == "No":
            break
    return fasta_filenames

def select_all_fastas(experiment_names: list[str]):
    experiment_fasta_files = {}
    for name in experiment_names:
        experiment_fasta_files[name] = select_experiment_fastas(experiment_fasta_files,name)
    return experiment_fasta_files

def print_fastas_by_experiment(experiment_fasta_files):
    for exp,fastas in experiment_fasta_files.items():
        print(f"[bold blue]{exp}[/]:")
        for fasta in fastas:
            print('\t',end='')
            if fasta == "":
                print("[red]Ignored[/]")
            else:
                print("[bold green]"+Path(fasta).name+"[/]")

def select_fasta_files(context: dict):
    try:
        experiment_fasta_files = context['experiment_fasta_files']
    except KeyError:
        experiment_fasta_files = select_all_fastas(context['experiment_names'])
        context['experiment_fasta_files'] = experiment_fasta_files

    while True:
        say("[bold yellow]Please review the following peptide data sets for accuracy:[/]")
        print_fastas_by_experiment(experiment_fasta_files)

        choice = q.select(
            "Are the listed FASTAs accurate?",
            choices=[
                "Yes, continue to analysis",
                "No, I'd like to make a change",
                "Forget it all!"
            ]
        ).ask()

        if "Forget" in choice:
            return EXIT_FAILURE
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
            experiment_fasta_files = select_all_fastas(context['experiment_names'])
            context['experiment_fasta_files'] = experiment_fasta_files
            continue

        if "Replace" in choice:
            # Pick a FASTA file to replace:
            choice = q.select(
                "Which FASTA to replace?",
                choices=util.unique_flattened_nested_str_list(list(experiment_fasta_files.values()))+["Nevermind, we're all good"]
            ).ask()

            if "Nevermind" in choice:
                # Assuming no FASTA will ever have "Nevermind" as part of the name...
                continue

            # Need to update all instances of the given FASTA file.
            #  indices = [(i,j) for j, fasta in enumerate(fasta_sublist) for i, fasta_sublist in enumerate(fasta_filenames) if fasta == choice]
            indices = [(name, i) for name, fasta_list in experiment_fasta_files.items()
                                for i, fasta in enumerate(fasta_list)
                                if fasta == choice]


            new_fasta_filename = files.open_file(f"Select FASTA file for {choice}",[("FASTA files",".FASTA .fasta .fas"),("All files",".*")])
            for k,i in indices:
                experiment_fasta_files[k][i] = new_fasta_filename
            context['experiment_fasta_files'] = experiment_fasta_files
            continue

        if "Edit" in choice:
            # Pick a column name to edit
            choice = q.select(
                "Which column to edit?",
                choices=list(experiment_fasta_files)+["Nevermind, we're all good"]
            ).ask()

            if "Nevermind" in choice:
                continue

            new_fasta_filelist = select_experiment_fastas(experiment_fasta_files,choice)
            experiment_fasta_files[choice] = new_fasta_filelist

    context['experiment_fasta_files'] = experiment_fasta_files
    context['step'] = 'read-fasta-sequences'

    return EXIT_SUCCESS
