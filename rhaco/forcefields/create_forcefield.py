import copy
from rhaco.definitions import ATOM_MASSES


def update_EAM_forcefield(file_name, elements_list):
    # Load the forcefield
    with open(file_name, "r") as original_FF:
        FF_lines = original_FF.readlines()

    # Number of Elements Line
    n_elements = FF_lines[3]
    # Update the element line
    new_n_elements = len(elements_list)
    FF_lines[3] = "{0} {1}\n".format(len(elements_list), " ".join(elements_list))

    # Number of Points Line
    n_points = FF_lines[4]
    split_n_points = n_points.split()
    nrho = int(split_n_points[0])
    drho = float(split_n_points[1])
    nr = int(split_n_points[2])
    dr = float(split_n_points[3])
    cutoff = float(split_n_points[4])

    # Complete header
    file_header = FF_lines[:5]

    # First element line: Atomic Number, Atomic Mass (AMU), Lattice Constant (Ang), Lattice Type (fcc, bcc, etc.)
    first_element_line = [FF_lines[5]]

    first_FF_line = FF_lines[6]
    numbers_per_line = len(first_FF_line.split())
    # The embedding function is first, there should be nrho numbers 
    # (nrho / numbers_per_line lines starting at 6)
    FF_starts_line = 6
    emb_func_fin_before = FF_starts_line + (nrho // numbers_per_line)
    first_element_emb_func = FF_lines[FF_starts_line:emb_func_fin_before]
    # Next is the density function, which has nr values
    # (nr / numbers_per_line lines starting at (nrho / numbers_per_line + 6))
    dens_func_fin_before = emb_func_fin_before + (nr // numbers_per_line)
    first_element_dens_func = FF_lines[emb_func_fin_before:dens_func_fin_before]
    # Then come the potential values in the pattern: 1-1, 1-2,...1-N.
    # Again there are nr values
    FF_ends_line = dens_func_fin_before + (nr // numbers_per_line)
    first_element_potentials = FF_lines[dens_func_fin_before:FF_ends_line]
    print(first_element_emb_func[0], "...", first_element_emb_func[-1])
    print(first_element_dens_func[0], "...", first_element_dens_func[-1])
    print(first_element_potentials[0], "...", first_element_potentials[-1])

    # Now create the additional element lines
    new_element_lines = create_element_lines(
        first_element_line, elements_list
    )

    # Now we want to add in a ton of zeroes corresponding to the other interactions
    blank_line = [" ".join(["{:24.16E}" * numbers_per_line, "\n"]).format(
        *[0.0] * numbers_per_line
    )]

    # LAMMPS documentation http://www.afs.enea.it/software/lammps/doc15/pair_eam.html
    # and cross-checking with *.eam.alloy files available through NIST:
    # The alloy file format for 3 elements looks like this:
    # Element 1 indentifier line
    # Element 1 emb fn
    # Element 1 dens fn
    # Element 2 indentifier line
    # Element 2 emb fn
    # Element 2 dens fn
    # Element 3 indentifier line
    # Element 3 emb fn
    # Element 3 dens fn
    # All potentials ([1, 1], [1, 2], [1, 3], [2, 2], [2, 3], [3, 3])
    # For n elements, the number of potential sets we have is the nth triangle number
    # = (n^{2} + n) / 2

    # For the potentials:
    # We already have E1-E1 interactions in first_element_potentials
    potentials = copy.deepcopy(first_element_potentials)
    # Add all the other interactions to forcefield
    for _ in range(int(((new_n_elements ** 2) + new_n_elements) / 2) - 1):
        potentials += blank_line * (nr // numbers_per_line)

    # For the embedding function and the density function
    embedding_function_lines = []
    density_function_lines = []
    for _ in range(new_n_elements):
        embedding_function_lines.append(blank_line * (nrho // numbers_per_line))
        density_function_lines.append(blank_line * (nrho // numbers_per_line))

    # Now reformat the EAM components
    print("Reformatting the EAM file components...")
    reformatted_EAM = reformat_EAM_components(
        file_header, new_element_lines, embedding_function_lines,
        density_function_lines,
    )


def create_element_lines(first_element_line, elements_list):
    # OLD CODE:
    # # Al
    # second_element_line = ["   ".join(["13", "{:24.16E}".format(26.981539), "{:24.16E}".format(0.0), "fcc\n"])]
    # # O
    # third_element_line = ["   ".join(["8", "{:24.16E}".format(15.999), "{:24.16E}".format(0.0), "fcc\n"])]
    additional_element_lines = first_element_line
    for element in elements_list:
        additional_element_lines.append(
            ["   ".join(
                [
                    "13", "{:24.16E}".format(ATOM_MASSES[element]),
                    "{:24.16E}".format(0.0), "fcc\n"
                ]
            )]
        )
    return additional_element_lines


def reformat_EAM_components(
    file_header, element_lines, embed_func_lines, dens_func_lines,
    new_numbers_per_line=4,
):
    new_FF_lines = file_header
    for element_ID, element_line in enumerate(element_lines):
        element_floats = [
            list(map(float, embed_func_lines[element_ID][:-1].split())),
            list(map(float, dens_func_lines[element_ID][:-1].split())),
        ]
        
        
        ***THIS IS WHERE I AM UP TO ***
        new_line_format = "".join(["{:20.12E}" * new_numbers_per_line, "\n"])
        new_FF_lines += element_names[elementID]
        # Create a set of iterators of dynamic length to iterate through the number list
        # effectively creating a rolling `window' of new_numbers_per_line to format
        iterators = [[iter(element_floats[i::new_numbers_per_line]) for i in range(new_numbers_per_line)]]
        line_of_floats = [list(zip(*iterator)) for iterator in iterators][0]
        for values in line_of_floats:
            new_FF_lines += new_line_format.format(*values)




if __name__ == "__main__":
    update_EAM_forcefield("Ag_Zhou04.eam.alloy", ["Ag", "Al", "O"])




    new_FF_lines = file_header
    new_numbers_per_line = 4
    element_names = {0: first_element_line,
                     1: second_element_line,
                     2: third_element_line}
    # Based on the *.eam.alloy format, here is the correct ordering:
    element_properties = {0: [first_element_emb_func, first_element_dens_func],
                          1: [second_element_emb_func, second_element_dens_func],
                          2: [third_element_emb_func, third_element_dens_func, potentials],
                         }
    for elementID in range(3):
        element_floats = []
        for element_props in element_properties[elementID]:
            for line in element_props:
                element_floats += list(map(float, line[:-1].split()))

    # print("Concatenating the eam file components...")
    # # Now we have all the components needed to write the file:
    # new_FF_lines = file_header + first_element_line + first_element_emb_func + first_element_dens_func + first_element_potentials + second_element_line + second_element_emb_func + second_element_dens_func + second_element_potentials + third_element_line + third_element_emb_func + third_element_dens_func + third_element_potentials

    print("Writing the new eam file...")
    #file_name = "corundum.eam.fs"
    file_name = "corundum_Zhou04.eam.alloy"
    with open(file_name, "w+") as new_FF_file:
        new_FF_file.writelines(new_FF_lines)
    print("EAM file written to", file_name)
