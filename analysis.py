import pandas as pd

def mean_expression_from_csv(csv_path: str) -> float:
    df = pd.read_csv(csv_path)
    if "expression" not in df.columns:
        raise KeyError("CSV must contain an 'expression' column")
    return float(df["expression"].mean())

if __name__ == "__main__":
    import sys

    if len(sys.argv) != 2:
        raise SystemExit("Usage: python analysis.py <path-to-csv>")

    print(mean_expression_from_csv(sys.argv[1]))