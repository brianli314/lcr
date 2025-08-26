import csv
import sys
from collections import defaultdict, Counter

def read_tsv(path):
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        records = []
        for row in reader:
            records.append((row["Name"], int(row["Start"]), int(row["End"])))
        return records

def compare_tsv(file1, file2):
    records1 = read_tsv(file1)
    records2 = read_tsv(file2)

    # Group by name
    groups1 = defaultdict(list)
    groups2 = defaultdict(list)
    for name, start, end in records1:
        groups1[name].append((start, end))
    for name, start, end in records2:
        groups2[name].append((start, end))

    all_names = set(groups1.keys()) | set(groups2.keys())

    differences = False
    for name in sorted(all_names):
        coords1 = groups1.get(name, [])
        coords2 = groups2.get(name, [])

        if not coords1:
            print(f"{name} present only in {file2}")
            differences = True
            continue
        if not coords2:
            print(f"{name} present only in {file1}")
            differences = True
            continue

        # Count frequencies of (start, end) pairs
        c1 = Counter(coords1)
        c2 = Counter(coords2)

        if c1 == c2:
            continue  # perfectly matches

        differences = True
        print(f"Differences for {name}:")
        for pair in set(c1.keys()) | set(c2.keys()):
            n1 = c1.get(pair, 0)
            n2 = c2.get(pair, 0)
            if n1 != n2:
                print(f"  Coordinates {pair}: {file1} has {n1}, {file2} has {n2}")

    if not differences:
        print("The two TSV files match exactly (ignoring ordering).")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: python compare_tsv.py file1.tsv file2.tsv")
        sys.exit(1)
    compare_tsv(sys.argv[1], sys.argv[2])
