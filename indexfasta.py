with open("random_1mbp.fasta", 'r') as f:
    string = ""
    for line in f:
        if line.startswith(">"):
            print("header")
            continue
        else:
            string = string + line.strip("\n")

    print(string[886503:886530])
    print(string[412273:412302])

    print(string[511502:511521])
    print(string[511501:511530])

    print(string[568128:568158])
    print(string[568159:568174])

    print(string[568128:568174])
