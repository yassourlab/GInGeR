from pathlib import Path

def get_filedir() -> str:
    currentdir = Path(__file__).resolve().parent
    return f"{currentdir}/test_files"