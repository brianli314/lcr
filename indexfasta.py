import sys
file_path = sys.argv[1]
start = int(sys.argv[2])
end = int(sys.argv[3])

with open(file_path, 'r') as f:
    string = ""
    for line in f:
        if line.startswith(">"):
            continue
        else:
            string = string + line.strip("\n")

    print(f"{len(string)} bp")
    print(string[start:end])
