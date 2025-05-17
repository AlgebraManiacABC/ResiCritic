from pathlib import Path
from tkinter import filedialog, messagebox
import pandas as pd
import tui
from version import __version__

def open_file(title="Select file",filetypes=None) -> str:
    filetypes_final = filetypes.copy()
    filetypes_final.append(("All files","*.*"))
    file_path = filedialog.askopenfilename(title=title,filetypes=filetypes_final)
    return file_path

def open_new_file(title="Save as...",extension=None) -> str:
    return filedialog.asksaveasfilename(title=title,extension=extension)

def get_excel_filename(argv: list[str]) -> str:
    try:
        filename = argv[1]
    except IndexError:
        filename = ""
        while not filename:
            filename = open_file(filetypes=[("Excel files","*.xlsx *.xls")])
            if not filename:
                give_up = messagebox.askyesno("Quit?","No Excel file was selected. Abort?")
                if give_up:
                    break
    return filename

def parse_excel(filename: str="") -> pd.DataFrame:
    try:
        file = pd.ExcelFile(filename)
    except FileNotFoundError:
        messagebox.showerror(f"{filename} does not exist!")
        return None
    except ValueError as ve:
        messagebox.showerror(f"There was an error parsing the Excel file '{filename}': {ve}")
        return None
    except Exception as e:
        messagebox.showerror(f"An unexpected error occurred: {e}")
        return None

    try:
        sheet_name = tui.choose_between(choices=file.sheet_names+["[Abort]"],ask="Which sheet to analyze?")
        if sheet_name == "[Abort]":
            return None
        return file.parse(sheet_name=sheet_name)
    finally:
        file.close()

def read_text(filename: str):
    return Path(filename).read_text()
