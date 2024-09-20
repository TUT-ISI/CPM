import itertools
import subprocess

# Define the range of values for VerticalSpeed and number concentration for each mode
vertical_speed_values = [0.2,0.3,0.4,0.5, 0.6, 0.7]
#vertical_speed_values = [0.9,1.3,1.6,1.9, 2.0, 2.7]
mode1_number_concentration_values = [100]
mode2_number_concentration_values = [100,200,300, 400, 500, 600]

# Iterate over each combination for Mode 1, Mode 2, and vertical speed
for vertical_speed, mode1_number_concentration, mode2_number_concentration in itertools.product(vertical_speed_values, mode1_number_concentration_values, mode2_number_concentration_values):
    # Convert number concentrations to string format
    mode1_number_concentration_str = "{:.2E}".format(mode1_number_concentration)
    mode2_number_concentration_str = "{:.2E}".format(mode2_number_concentration)
    #mode1_number_concentration_str = float(mode1_number_concentration)
    #mode2_number_concentration_str = float(mode2_number_concentration)

    # Read the environ template file
    with open("input/environ_template.dat", "r") as f:
        environ_template = f.read()

    # Replace the placeholders in environ_template for vertical speed
    environ_content_mode1 = environ_template.format(vertical_speed)

    # Write the modified content to a temporary file for Mode 1
    with open("input/environ.dat", "w") as f:
        f.write(environ_content_mode1)

    # Read the aerosol template file
    with open("input/aerosol_template.dat", "r") as f:
        aerosol_template = f.read()

    # Replace the placeholders in aerosol_template for Mode 1 and Mode 2
    aerosol_content = aerosol_template.format(mode1_number_concentration_str, mode2_number_concentration_str)

    # Write the modified content to the same file for both modes
    with open("input/aerosol.dat", "w") as f:
        f.write(aerosol_content)

    # Run the Fortran code for each combination of Mode 1 and Mode 2
    subprocess.run(["./cmodel"], check=True)
