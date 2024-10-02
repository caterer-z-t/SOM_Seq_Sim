import pandas as pd
import os


class FileLoader:
    def __init__(self, file_path):
        """Initializes the FileLoader object with the provided file path.

        Args:
            file_path (str): The path to the file to be loaded.

        Raises:
            FileNotFoundError: If the file does not exist.
            ValueError: If the file extension is unsupported.
        """
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"The file {file_path} does not exist.")

        self.file_path = file_path
        self.file_extension = os.path.splitext(file_path)[1].lower()

    def load_file(self):
        """Loads the file based on its extension and converts it to a pandas DataFrame.

        Returns:
            pd.DataFrame: A pandas DataFrame with the file's contents.

        Raises:
            ValueError: If the file extension is not supported.
        """
        if self.file_extension == ".csv":
            return self._load_csv()
        elif self.file_extension in [".xls", ".xlsx"]:
            return self._load_excel()
        elif self.file_extension == ".json":
            return self._load_json()
        elif self.file_extension == [".txt", ".tsv"]:
            return self._load_txt()
        else:
            raise ValueError(f"Unsupported file format: {self.file_extension}")

    def _load_csv(self):
        """Loads a CSV file into a DataFrame."""
        try:
            return pd.read_csv(self.file_path)
        except Exception as e:
            print(f"Error loading CSV file: {e}")
            return None

    def _load_excel(self):
        """Loads an Excel file into a DataFrame."""
        try:
            return pd.read_excel(self.file_path)
        except Exception as e:
            print(f"Error loading Excel file: {e}")
            return None

    def _load_json(self):
        """Loads a JSON file into a DataFrame."""
        try:
            return pd.read_json(self.file_path)
        except Exception as e:
            print(f"Error loading JSON file: {e}")
            return None

    def _load_txt(self, delimiter="\t"):
        """Loads a text file (e.g., tab-delimited or space-delimited) into a DataFrame.

        Args:
            delimiter (str, optional): The delimiter used in the text file. Defaults to "\t" (tab).
        """
        try:
            return pd.read_csv(self.file_path, delimiter=delimiter)
        except Exception as e:
            print(f"Error loading text file: {e}")
            return None
