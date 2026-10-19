# imports
# ----------------------------------------------------------------------------------------------------------------------
import json
from pathlib import Path
# ----------------------------------------------------------------------------------------------------------------------


# JSON
########################################################################################################################
def get_json_dict(json_file_path: (str | Path)) -> dict:
    """
    Parameters
    ----------
    json_file_path : str | Path
        Absolute-path of JSON file.

    Returns
    -------
    dict
        JSON file as dictionary.
    """

    if not Path(json_file_path).exists():
        raise FileExistsError("JSON File could not be found!")
    with open(json_file_path, 'r') as in_file:
        return json.load(in_file)


def save_as_json_dict(dict_saving: dict, json_file_path: str, name_file: str) -> None:
    """
    Parameters
    ----------
    dict_saving : dict
        Dictionary that should be written to JSON file.
    json_file_path : str
        Target directory for saving.
    name_file : str
        Name of file; written without the '.json'.
    """

    if Path(json_file_path).exists() and isinstance(dict_saving, dict):
        target_file = Path(json_file_path) / f"{name_file}.json"
        if not target_file.exists():
            target_file.touch()
            print(f"File {target_file} has been created...")
        with open(target_file, "w") as out_file:
            json.dump(dict_saving, out_file)
    else:
        raise ValueError("Inputted path or dictionary are incorrect!")

