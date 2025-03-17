import argparse

def filter_pdb(input_file, output_file, remove_connect):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith("HETATM") and "HEM" not in line:
                continue
            if remove_connect and line.startswith("CONECT"):
                continue
            outfile.write(line)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Filter a PDB file by removing non-HEM HETATM lines and optionally removing CONECT lines."
    )
    parser.add_argument("input_file", help="Path to the input PDB file")
    parser.add_argument("output_file", help="Path for the output filtered PDB file")
    parser.add_argument("--remove_connect", action="store_true",
                        help="If set, remove all CONECT lines")
    
    args = parser.parse_args()
    filter_pdb(args.input_file, args.output_file, args.remove_connect)