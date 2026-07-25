import re
import os

# Read input file
params = {}

with open("input_PFM_E.in", "r") as f:
    for line in f:
        line = line.strip()

        if "=" in line:
            key, value = line.split("=", 1)
            params[key.strip()] = value.strip().replace(";", "")

# Read values
MESH_X = params["MESH_X"]
MESH_Y = params["MESH_Y"]
MESH_Z = params["MESH_Z"]

DELTA_X = params["DELTA_X"]
DELTA_t = params["DELTA_t"]

STARTTIME = params["STARTTIME"]
NTIMESTEPS = params["NTIMESTEPS"]
SAVET = params["SAVET"]

BARRIER_HEIGHT = params["W"]
GRADIENT_ENERGY_COEFF = params["k"]
INTERFACE_MOBILITY = params["L_sigma"]
KINETIC_COEFF = params["L_eta"]
ALPHA = params["alpha"]
DIFFUSIVITY = params["D"]
BULK_ELECTROLYTE_CONCENTRATION = params["C_0"]
SITE_DENSITY_ELECTRODE = params["C_m_s"]
SITE_DENSITY_ELECTROLYTE = params["C_m_l"]
ELECTRODE_CONDUCTIVITY = params["Sigma_s"]
ELECTROLYTE_CONDUCTIVITY = params["Sigma_l"]
EPSILON_ELECTRODE = params["epsilon_s"]
EPSILON_ELECTROLYTE = params["epsilon_l"]
NUM_MOLES = params["n"]
FARADAYS_CONSTANT = params["F"]
GAS_CONSTANT = params["R"]
TEMPERATURE = params["T"]

X_BOX = params["X_BOX"]
Y_BOX = params["Y_BOX"]
Z_BOX = params["Z_BOX"]

ETA_DEFAULT = params["ETA_DEFAULT"]
MU_DEFAULT = params["MU_DEFAULT"]
PHI_DEFAULT = params["PHI_DEFAULT"]
C_DEFAULT = params["C_DEFAULT"]

ETA_INIT = params["ETA_INIT"]
MU_INIT = params["MU_INIT"]
PHI_INIT = params["PHI_INIT"]
C_INIT = params["C_INIT"]

# Write parameters to generatedInput file
with open("generatedInput", "w") as f:
    f.write(f"MESH_X {MESH_X};\n")
    f.write(f"MESH_Y {MESH_Y};\n")
    f.write(f"MESH_Z {MESH_Z};\n")

    f.write(f"DELTA_X {DELTA_X};\n")
    f.write(f"DELTA_t {DELTA_t};\n")

    f.write(f"STARTTIME {STARTTIME};\n")
    f.write(f"NTIMESTEPS {NTIMESTEPS};\n")
    f.write(f"SAVET {SAVET};\n")

    f.write(f"BARRIER_HEIGHT {BARRIER_HEIGHT};\n")
    f.write(f"GRADIENT_ENERGY_COEFF {GRADIENT_ENERGY_COEFF};\n")
    f.write(f"INTERFACE_MOBILITY {INTERFACE_MOBILITY};\n")
    f.write(f"KINETIC_COEFF {KINETIC_COEFF};\n")
    f.write(f"ALPHA {ALPHA};\n")
    f.write(f"DIFFUSIVITY {DIFFUSIVITY};\n")
    f.write(f"BULK_ELECTROLYTE_CONCENTRATION {BULK_ELECTROLYTE_CONCENTRATION};\n")
    f.write(f"SITE_DENSITY_ELECTRODE {SITE_DENSITY_ELECTRODE};\n")
    f.write(f"SITE_DENSITY_ELECTROLYTE {SITE_DENSITY_ELECTROLYTE};\n")
    f.write(f"ELECTRODE_CONDUCTIVITY {ELECTRODE_CONDUCTIVITY};\n")
    f.write(f"ELECTROLYTE_CONDUCTIVITY {ELECTROLYTE_CONDUCTIVITY};\n")
    f.write(f"EPSILON_ELECTRODE {EPSILON_ELECTRODE};\n")
    f.write(f"EPSILON_ELECTROLYTE {EPSILON_ELECTROLYTE};\n")
    f.write(f"NUM_MOLES {NUM_MOLES};\n")
    f.write(f"FARADAYS_CONSTANT {FARADAYS_CONSTANT};\n")
    f.write(f"GAS_CONSTANT {GAS_CONSTANT};\n")
    f.write(f"TEMPERATURE {TEMPERATURE};\n")

print("Parameters Updated.")

# OpenFOAM boundary conditions
# OpenFOAM patch names
patch_names = [
    "floor",
    "boxtoCell",
    "ceiling",
    "electrode",
    "electrolyte"
]

# Mapping
bc_type = {
    0: "symmetryPlane",
    1: "zeroGradient",
    2: "fixedValue",
    3: "cyclic"
}

# Read Input.in
boundary = {}
boundary_value = {}

with open("input_PFM_E.in","r") as f:

    for line in f:

        line=line.strip()

        if line.startswith("BOUNDARY ="):

            s=line[line.find("{")+1:line.find("}")]
            data=[x.strip() for x in s.split(",")]

            field=data[0]
            boundary[field]=list(map(int,data[1:]))

        elif line.startswith("BOUNDARY_VALUE ="):

            s=line[line.find("{")+1:line.find("}")]
            data=[x.strip() for x in s.split(",")]

            field=data[0]
            boundary_value[field]=list(map(float,data[1:]))

for field in boundary.keys():

    filename = "0/" + field + ".orig"

    if not os.path.exists(filename):
        print(filename, "not found")
        continue

    with open(filename, "r") as f:
        lines = f.readlines()

    new_lines = []

    inside_patch = False
    patch_index = -1
    value_written = False

    i = 0

    while i < len(lines):

        line = lines[i]
        stripped = line.strip()

        # Detect patch
        if stripped in patch_names:

            inside_patch = True
            patch_index = patch_names.index(stripped)
            value_written = False

            new_lines.append(line)
            i += 1
            continue

        # Replace BC type
        if inside_patch and stripped.startswith("type"):

            new_type = bc_type[boundary[field][patch_index]]

            new_lines.append(
                f"        type            {new_type};\n"
            )

            i += 1
            continue

        # Existing value line
        if inside_patch and stripped.startswith("value"):

            if bc_type[boundary[field][patch_index]] == "fixedValue":

                value = boundary_value[field][patch_index]

                new_lines.append(
                    f"        value           uniform {value};\n"
                )

                value_written = True

            # Skip old value line
            i += 1
            continue

        # Before closing brace
        if inside_patch and stripped == "}":

            if (bc_type[boundary[field][patch_index]] == "fixedValue"
                    and not value_written):

                value = boundary_value[field][patch_index]

                new_lines.append(
                    f"        value           uniform {value};\n"
                )

            inside_patch = False
            patch_index = -1
            value_written = False

            new_lines.append(line)

            i += 1
            continue

        new_lines.append(line)

        i += 1

    with open(filename, "w") as f:
        f.writelines(new_lines)

    print(field, "boundary updated.")


# setFieldsDict (only boxToCell)

with open("system/setFieldsDict", "r") as f:
    text = f.read()

# Update defaultFieldValues
text = re.sub(
    r'(volScalarFieldValue\s+eta\s+)[-+0-9.eE]+',
    rf'\g<1>{ETA_DEFAULT}',
    text,
    count=1
)

text = re.sub(
    r'(volScalarFieldValue\s+mu\s+)[-+0-9.eE]+',
    rf'\g<1>{MU_DEFAULT}',
    text,
    count=1
)

text = re.sub(
    r'(volScalarFieldValue\s+phi\s+)[-+0-9.eE]+',
    rf'\g<1>{PHI_DEFAULT}',
    text,
    count=1
)

text = re.sub(
    r'(volScalarFieldValue\s+c\s+)[-+0-9.eE]+',
    rf'\g<1>{C_DEFAULT}',
    text,
    count=1
)

# Find the boxToCell block
pattern = r'(boxToCell\s*\{.*?\})'
match = re.search(pattern, text, re.DOTALL)

if match:

    box_block = match.group(1)

    # Replace box coordinates
    box_block = re.sub(
    r'box\s*\(\s*0\s+0\s+0\s*\)\s*\(\s*[-+0-9.eE]+\s+[-+0-9.eE]+\s+[-+0-9.eE]+\s*\);',
    f'box (0 0 0) ({X_BOX} {Y_BOX} {Z_BOX});',
    box_block
    )

    # Replace only values inside boxToCell
    box_block = re.sub(
        r'(volScalarFieldValue\s+eta\s+)[-+0-9.eE]+',
        rf'\g<1>{ETA_INIT}',
        box_block
    )

    box_block = re.sub(
        r'(volScalarFieldValue\s+mu\s+)[-+0-9.eE]+',
        rf'\g<1>{MU_INIT}',
        box_block
    )

    box_block = re.sub(
        r'(volScalarFieldValue\s+phi\s+)[-+0-9.eE]+',
        rf'\g<1>{PHI_INIT}',
        box_block
    )

    box_block = re.sub(
        r'(volScalarFieldValue\s+c\s+)[-+0-9.eE]+',
        rf'\g<1>{C_INIT}',
        box_block
    )

    # Put modified block back
    text = text[:match.start()] + box_block + text[match.end():]

with open("system/setFieldsDict", "w") as f:
    f.write(text)

print("setFieldsDict updated.")