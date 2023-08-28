import re
acc_regex = re.compile("(?P<accession>[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})(?P<isotype>-\d+)?")
reg_positon_residue = re.compile("_(\w)(\d+)")


def replace_special_with_dot(s: str):
    return re.sub(r'[^a-zA-Z0-9.]', '.', s)



# a function to read fasta file into a dictionary where the key is the uniprot accession id extracted from the fasta
def read_fasta(fasta_file: str) -> dict[str, str]:
    fasta_dict = {}

    current_acc = ""
    for line in fasta_file.split("\n"):
        line = line.strip()
        if line.startswith('>'):
            match = acc_regex.search(line.replace(">", ""))
            if match:
                acc = match.groupdict(default="")["accession"]
                iso = match.groupdict(default="")["isotype"]
                fasta_dict[acc+iso] = ""
                current_acc = acc+iso

            else:
                fasta_dict[line.replace(">", "")] = ""
                current_acc = line.replace(">", "")

        else:
            fasta_dict[current_acc] += line
    return fasta_dict