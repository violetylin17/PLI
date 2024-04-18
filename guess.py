def guess_atom_name(original_name):
    new_name = original_name.strip()

    for index, char in enumerate(new_name):
        if char.isdigit():
            new_name = new_name[0:index]
            break
    if new_name in ["BR", "CL", "MG"]:
        new_name = new_name[0:2]
    else:
        new_name = new_name[0]

    return new_name
