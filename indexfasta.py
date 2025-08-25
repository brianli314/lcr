with open("test.fasta", 'r') as f:
    string = ""
    for line in f:
        if line.startswith(">"):
            print("header")
            continue
        else:
            string = string + line.strip("\n")

    print(string[0:10])
