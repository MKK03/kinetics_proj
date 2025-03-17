#!/usr/bin/env python3
import argparse

def parse_pdb_compute_com_and_max_serial(lines):
    """
    Compute the geometric center (center of mass) of all ATOM records in the PDB file.
    Also determine the maximum atom serial number.
    """
    sum_x = sum_y = sum_z = 0.0
    count = 0
    max_serial = 0

    for line in lines:
        if line.startswith("ATOM"):
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                sum_x += x
                sum_y += y
                sum_z += z
                count += 1
            except ValueError:
                continue

            try:
                serial = int(line[6:11].strip())
                if serial > max_serial:
                    max_serial = serial
            except ValueError:
                continue
        elif line.startswith("HETATM"):
            # also update max_serial from HETATM lines
            try:
                serial = int(line[6:11].strip())
                if serial > max_serial:
                    max_serial = serial
            except ValueError:
                continue

    if count == 0:
        raise ValueError("No ATOM records found to compute center of mass.")
    com = (sum_x / count, sum_y / count, sum_z / count)
    return com, max_serial

def create_o2_atoms(com, distance, max_serial, o2_resseq=999, chain_id="L"):
    """
    Create two HETATM records for the O2 molecule.
    - The midpoint of the O2 molecule is placed at COM + (distance, 0, 0).
    - The O2 molecule is assumed to be linear with an O–O bond of 1.21 Å.
      Here, we place one oxygen above and one below the midpoint along the y-axis.
    """
    # displacement vector from COM (along x axis)
    com_x, com_y, com_z = com
    mid_x = com_x + distance
    mid_y = com_y
    mid_z = com_z

    # O2 bond parameters
    o2_bond = 1.21      # Å, typical O-O bond length
    half_bond = o2_bond / 2.0

    # Define positions for the two oxygen atoms.
    # We place them along y axis relative to the midpoint.
    O1_pos = (mid_x, mid_y + half_bond, mid_z)
    O2_pos = (mid_x, mid_y - half_bond, mid_z)

    # Format HETATM records.
    # PDB fixed width format for HETATM:
    # Columns: Record, atom serial, atom name, altLoc, resName, chainID, resSeq, iCode,
    #          x, y, z, occupancy, tempFactor, element, charge
    # We'll use a simple format string.
    hetatm_format = (
        "HETATM{atom_num:5d} {atom_name:^4s} {res_name:>3s} {chain_id:1s}"
        "{res_seq:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           {element:>2s}\n"
    )

    o2_res_name = "O2"
    # Create two lines
    line1 = hetatm_format.format(atom_num=max_serial + 1,
                                 atom_name="O1",
                                 res_name=o2_res_name,
                                 chain_id=chain_id,
                                 res_seq=o2_resseq,
                                 x=O1_pos[0], y=O1_pos[1], z=O1_pos[2],
                                 element="O")
    line2 = hetatm_format.format(atom_num=max_serial + 2,
                                 atom_name="O2",
                                 res_name=o2_res_name,
                                 chain_id=chain_id,
                                 res_seq=o2_resseq,
                                 x=O2_pos[0], y=O2_pos[1], z=O2_pos[2],
                                 element="O")
    return line1, line2

def insert_o2(input_file, output_file, distance):
    """
    Reads the input PDB file, computes the center of mass from ATOM lines,
    and inserts an O2 molecule (as two HETATM lines) at a distance along the x-axis.
    The new lines are inserted just before the "END" record if present.
    """
    with open(input_file, "r") as f:
        lines = f.readlines()

    # Compute COM and get the max atom serial number
    com, max_serial = parse_pdb_compute_com_and_max_serial(lines)
    print("Center of Mass (computed from ATOM records):", com)
    print("Max atom serial number:", max_serial)

    # Create new O2 HETATM lines.
    o2_line1, o2_line2 = create_o2_atoms(com, distance, max_serial)

    # Insert the new lines before the END record if it exists.
    end_index = None
    for i, line in enumerate(lines):
        if line.strip() == "END":
            end_index = i
            break

    if end_index is not None:
        new_lines = lines[:end_index] + [o2_line1, o2_line2] + lines[end_index:]
    else:
        new_lines = lines + [o2_line1, o2_line2]

    with open(output_file, "w") as f:
        f.writelines(new_lines)
    print(f"New PDB with O2 inserted saved as {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Insert an O2 molecule (as two HETATM records) into a PDB file. "
                    "The O2 molecule is placed on a plane at a user-specified distance from the protein's center of mass."
    )
    parser.add_argument("input_file", help="Input PDB file path")
    parser.add_argument("output_file", help="Output PDB file path with O2 inserted")
    parser.add_argument("--distance", type=float, default=50.0,
                        help="Distance (in Å) from the protein's center of mass to the plane where O2 is placed (default: 50.0 Å)")
    args = parser.parse_args()

    insert_o2(args.input_file, args.output_file, args.distance)