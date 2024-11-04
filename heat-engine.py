import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.animation import FuncAnimation
import numpy as np

# Constants
g = 9.81  # Acceleration due to gravity in m/s²
mass = 0.05  # Mass of the weight in kg (50 g)
piston_area = 0.001  # Area of the piston in m² (arbitrary but fixed for simplicity)
start_height = 0.005  # Initial height of the weight in meters (5 mm)
pressure = 101325 + mass/piston_area  # Atmospheric pressure in Pascals (1 atm)

gas_properties = {
    "air": {
        "cv": 717,      # Specific Cv for air in J/kg·K
        "cp": 1005,     # Specific Cp for air in J/kg·K
        "density": 1.225 # Density of air in kg/m³
    },
    "helium": {
        "cv": 3116,     # Specific Cv for helium in J/kg·K
        "cp": 5193,     # Specific Cp for helium in J/kg·K
        "density": 0.1785 # Density of helium in kg/m³
    },
    "argon": {
        "cv": 312,      # Specific Cv for argon in J/kg·K
        "cp": 520,      # Specific Cp for argon in J/kg·K
        "density": 1.784 # Density of argon in kg/m³
    },
    "oxygen": {
        "cv": 658,      # Specific Cv for oxygen in J/kg·K
        "cp": 918,      # Specific Cp for oxygen in J/kg·K
        "density": 1.429 # Density of oxygen in kg/m³
    },
    "nitrogen": {
        "cv": 743,     # Specific Cv for nitrogen in J/kg·K
        "cp": 1039,     # Specific Cp for nitrogen in J/kg·K
        "density": 1.250 # Density of nitrogen in kg/m³
    }
}

# Animation Function
def animate_piston(T_initial, T_final, distance_moved, start_height=5):
    # Setup the figure
    fig, ax = plt.subplots()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 50)  # Cylinder height from 0 to 50 mm
    
    # Draw the cylinder (rectangle)
    cylinder = patches.Rectangle((0.2, 0), 0.6, start_height, linewidth=2, edgecolor='black', facecolor='white')
    ax.add_patch(cylinder)
    
    # Create the piston as a horizontal rectangle
    piston_height = start_height  # Initial piston height in mm
    piston = patches.Rectangle((0.2, piston_height), 0.6, 1, linewidth=2, edgecolor='black', facecolor='gray')
    ax.add_patch(piston)
    
    # Text annotations for temperature
    text_loc = 45
    if distance_moved > 40:
        text_loc = 55
    T_initial_c = T_initial - 273.15
    T_final_c = T_final - 273.15
    temp_text = ax.text(0.1, text_loc, f"Temp: {T_initial_c:.2f} C", fontsize=10)
    
    # Function to update the piston and gas color in each frame
    def update(frame):
        # Update the piston height
        new_height = piston_height + frame * (distance_moved / 100)
        piston.set_y(new_height)
        cylinder.set_height(new_height)
        
        # Update the gas color based on the temperature change
        temp_fraction = frame / 100  # Fraction of total temperature change
        current_temp = T_initial_c + temp_fraction * (T_final_c - T_initial_c)
        
        # Interpolate color from blue (cold) to red (hot)
        gas_color = (1, 0, 0, temp_fraction)  # RGBA where A controls transparency
        cylinder.set_facecolor(gas_color)
        
        # Update the temperature text
        temp_text.set_text(f"Temp: {current_temp:.2f} C")
    
    # Create the animation
    ani = FuncAnimation(fig, update, frames=np.linspace(0, 100, 100), repeat=False)
    
    # Show the animation
    plt.show()

# Heat engine calculations
def heat_engine(system_gas, T_cold, T_hot, V_initial, P, A_piston):
    # Determine gas properties
    Cp = gas_properties[system_gas]["cp"]
    Cv = gas_properties[system_gas]["cv"]
    rho = gas_properties[system_gas]["density"]
    m_gas = rho * V_initial

    # Calculate the volume change
    V_2 = V_initial * (T_hot / T_cold)
    delta_h = ((V_2 - V_initial) / A_piston) * 1000  # Convert to mm
    height_limited = False
    if delta_h > 45:
        delta_h = 45
        V_2 = delta_h/1000 * A_piston + V_initial
        height_limited = True
        W_out = P * (V_2 - V_initial)
        T_inter = T_cold * (V_2 / V_initial)
        delta_U = m_gas * Cp * (T_inter - T_cold)
        delta_U += m_gas * Cv * (T_hot - T_inter)
        Q_in = W_out + delta_U
    else:
        W_out = P * (V_2 - V_initial)
        delta_U = m_gas * Cp * (T_hot - T_cold)
        Q_in = W_out + delta_U

    # Calculate the efficiency of the heat engine and change in height
    eta = W_out / Q_in
    
    return W_out, Q_in, delta_U, eta, delta_h, height_limited

# Print the results
def print_report(W_out, Q_in, delta_U, eta, delta_h, height_limit, hide_efficiency=True):
    print("\nHeat Engine Simulation Results:")
    print(f"Work done by the engine: {W_out:.2f} J")
    print(f"Heat added to the engine: {Q_in:.2f} J")
    print(f"Change in internal energy of the gas: {delta_U:.2f} J")
    if not hide_efficiency:
        print(f"Efficiency of the engine: {eta:.2f}")
    if height_limit:
        print("The change in height for the piston was limited to 45 mm.")
    else:
        print(f"Change in height of the weight: {delta_h:.2f} mm\n")

# Example inputs (students can modify these values)
gas_options = ['air', 'helium', 'argon', 'oxygen', 'nitrogen']
system_gas = input("Enter the type of gas (air, helium, argon, oxygen, nitrogen): ").lower()
if system_gas not in gas_options:
    print("Invalid gas type. Please choose from: air, helium, argon, oxygen, nitrogen.")
    exit()
T_cold = float(input("Enter the cold bath temperature (C): ")) + 273.15
if T_cold < 273.15:
    print("Enter a temperature greater than 0 C.")
    exit()
T_hot = float(input("Enter the hot bath temperature (C): ")) + 273.15
if T_hot > 423.15:
    print("Enter a temperature below 150 C.")
    exit()

# Perform the simulation
work_done, heat_added, internal_energy_change, efficiency, distance_moved_mm, h_lim = heat_engine(system_gas, T_cold, T_hot, start_height * piston_area + 0.0001, pressure, piston_area)
print_report(work_done, heat_added, internal_energy_change, efficiency, distance_moved_mm, h_lim)

# Run the animation
animate_piston(T_cold, T_hot, distance_moved_mm)
