if __name__ == "__main__":
    # Create a hybrid 3-element forcefield based on the silver forcefield
    with open("./Ag.eam.fs", "r") as original_FF:
        FF_lines = original_FF.readlines()
    # Number of Elements Line
    n_elements = FF_lines[3]
    # Number of Points Line
    n_points = FF_lines[4]
    split_n_points = n_points.split()
    nrho = int(split_n_points[0])
    drho = float(split_n_points[1])
    nr = int(split_n_points[2])
    dr = float(split_n_points[3])
    cutoff = float(split_n_points[4])
    # Update the element line
    new_n_elements = 3
    elements = ["Ag", "Al", "O"]
    FF_lines[3] = "3 Ag Al O\n"
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
    #FF_lines[3] = " ".join([new_n_elements] + [elements])

    # Al
    second_element_line = ["   ".join(["13", "{:24.16E}".format(26.981539), "{:24.16E}".format(0.0), "fcc\n"])]
    # O
    third_element_line = ["   ".join(["8", "{:24.16E}".format(15.999), "{:24.16E}".format(0.0), "fcc\n"])]

    # Now we want to add in a ton of zeroes corresponding to the other interactions
    blank_line = [" ".join(["{:24.16E}" * 5, "\n"]).format(*[0.0]*5)]

    # Element 1:
    # Add E1-E2 and E1-E3 interactions to forcefield
    # Firstly, Element 1 interacting with Element 2
    first_element_potentials += blank_line * (nr // numbers_per_line)
    # Then, Element 1 interacting with Element 3
    first_element_potentials += blank_line * (nr // numbers_per_line)

    # Element 2:
    # Firstly, we need the embedding function
    second_element_emb_func = blank_line * (nrho // numbers_per_line)
    # Then, the density function
    second_element_dens_func = blank_line * (nrho // numbers_per_line)
    # Then, the forcefield (E2-E2, E2-E3, E2-E1)
    second_element_potentials = blank_line * (nr // numbers_per_line) * 3

    # Element 3:
    # Firstly, we need the embedding function
    third_element_emb_func = blank_line * (nrho // numbers_per_line)
    # Then, the density function
    third_element_dens_func = blank_line * (nrho // numbers_per_line)
    # Then, the forcefield (E2-E2, E2-E3, E2-E1)
    third_element_potentials = blank_line * (nr // numbers_per_line) * 3


    new_FF_lines = file_header
    print("Reformatting the eam file components...")
    new_numbers_per_line = 4
    element_names = {0: first_element_line,
                     1: second_element_line,
                     2: third_element_line}
    # element_properties = {0: [first_element_emb_func, first_element_dens_func, first_element_potentials],
    #                       1: [second_element_emb_func, second_element_dens_func, second_element_potentials],
    #                       2: [third_element_emb_func, third_element_dens_func, third_element_potentials]}
    # NOTE: Perhaps density function only needed for single element potentials?
    element_properties = {0: [first_element_emb_func, first_element_potentials],
                          1: [second_element_emb_func, second_element_potentials],
                          2: [third_element_emb_func, third_element_potentials],
                         }
    for elementID in range(3):
        element_floats = []
        for element_props in element_properties[elementID]:
            for line in element_props:
                element_floats += list(map(float, line[:-1].split()))
        new_line_format = "".join(["{:20.12E}" * new_numbers_per_line, "\n"])
        new_FF_lines += element_names[elementID]
        # Create a set of iterators of dynamic length to iterate through the number list
        # effectively creating a rolling `window' of new_numbers_per_line to format
        iterators = [[iter(element_floats[i::new_numbers_per_line]) for i in range(new_numbers_per_line)]]
        line_of_floats = [list(zip(*iterator)) for iterator in iterators][0]
        for values in line_of_floats:
            new_FF_lines += new_line_format.format(*values)

    # print("Concatenating the eam file components...")
    # # Now we have all the components needed to write the file:
    # new_FF_lines = file_header + first_element_line + first_element_emb_func + first_element_dens_func + first_element_potentials + second_element_line + second_element_emb_func + second_element_dens_func + second_element_potentials + third_element_line + third_element_emb_func + third_element_dens_func + third_element_potentials

    print("Writing the new eam file...")
    with open("corundum.eam.fs", "w+") as new_FF_file:
        new_FF_file.writelines(new_FF_lines)
    print("EAM file written to corundum.eam.fs")
