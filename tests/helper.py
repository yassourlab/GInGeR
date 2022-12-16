from pathlib import Path


def get_test_data_dir() -> str:
    return f"{Path(__file__).resolve().parent}/data"


def get_test_output_dir() -> str:
    return f"{Path(__file__).resolve().parent}/tests_outputs"
