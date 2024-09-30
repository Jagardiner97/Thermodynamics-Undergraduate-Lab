import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.animation import FuncAnimation
import numpy as np

# Constants
g = 9.81  # Acceleration due to gravity in m/s²
mass = 0.05  # Mass of the weight in kg (50 g)
piston_area = 0.01  # Area of the piston in m² (arbitrary but fixed for simplicity)

# Specific heat capacities at constant pressure (Cp) for different gases (in J/kg·K)
gas_cv = {
    "air": 717,      # Specific Cv for air in J/kg·K
    "helium": 5193,  # Specific Cv for helium in J/kg·K
    "argon": 520,    # Specific Cv for argon in J/kg·K
    "oxygen": 918,   # Specific Cv for oxygen in J/kg·K
    "nitrogen": 1040 # Specific Cv for nitrogen in J/kg·K
}

gas_cp = {
    "air": 1005,      # Specific Cp for air in J/kg·K
    "helium": 5193,   # Specific Cp for helium in J/kg·K
    "argon": 520,     # Specific Cp for argon in J/kg·K
    "oxygen": 918,    # Specific Cp for oxygen in J/kg·K
    "nitrogen": 1040  # Specific Cp for nitrogen in J/kg·K
}

gas_density = {
    "air": 1.225,      # Density of air in kg/m³
    "helium": 0.1785,  # Density of helium in kg/m³
    "argon": 1.784,    # Density of argon in kg/m³
    "oxygen": 1.429,   # Density of oxygen in kg/m³
    "nitrogen": 1.250  # Density of nitrogen in kg/m³
}

def calculate_work_and_distance(P, delta_V):
    """Calculate work done by gas and distance the piston moves."""
    # Work done by the gas (W = P * ΔV)
    work_done = P * delta_V
    
    # Calculate the distance moved by the piston (ΔV = area * distance)
    distance_moved = delta_V / piston_area  # Distance in meters
    
    return work_done, distance_moved

def first_law_thermodynamics(system_gas, T_environment, V_initial, Q_in, P):
    if system_gas not in gas_data:
        raise ValueError(f"Gas '{system_gas}' is not available. Choose from: {list(gas_data.keys())}")
    
    Cv = gas_cv[system_gas]  # Get specific Cv for the selected gas (J/kg·K)
    Cp = gas_cp[system_gas]  # Get specific Cp for the selected gas (J/kg·K)
    density = gas_density[system_gas]  # Get density of the selected gas (kg/m³)
    k = Cp / Cv  # Ratio of specific heat capacities (Cp/Cv)
    
    # Assume we have a constant amount of gas (say 1 kg for simplicity)
    gas_mass = density * V_initial # Mass of the gas in the cylinder (kg)
    
    # Calculate the change in internal energy (ΔU = Cv * ΔT * gas_mass)
    delta_U = Q_in  # Since no external work is done, ΔU = Q_in
    
    # Calculate the change in temperature using ΔU = Cv * ΔT * gas_mass (constant volume process)
    delta_T = delta_U / (Cv * gas_mass)
    
    # Final temperature
    T_final = T_environment + delta_T
    
    # Calculate work done by the gas (W = Q_in at constant pressure)
    delta_V = Q_in / P  # Work done by gas W = P * ΔV, so ΔV = W / P
    
    # Calculate distance moved by the piston
    work_done, distance_moved = calculate_work_and_distance(P, delta_V)
    
    # Final volume
    V_final = V_initial + delta_V
    
    # Convert distance to mm for better readability
    distance_moved_mm = distance_moved * 1000
    
    return T_final, V_final, distance_moved_mm, delta_U

# Animation Function
def animate_piston(T_initial, T_final, distance_moved):
    # Setup the figure
    fig, ax = plt.subplots()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 20)  # Cylinder height from 0 to 20 mm
    
    # Draw the cylinder (rectangle)
    cylinder = patches.Rectangle((0.2, 0), 0.6, 10, linewidth=2, edgecolor='black', facecolor='white')
    ax.add_patch(cylinder)
    
    # Create the piston as a horizontal rectangle
    piston_height = 10  # Initial piston height in mm
    piston = patches.Rectangle((0.2, piston_height), 0.6, 1, linewidth=2, edgecolor='black', facecolor='gray')
    ax.add_patch(piston)
    
    # Text annotations for temperature
    temp_text = ax.text(0.1, 18, f"Temp: {T_initial:.2f} K", fontsize=10)
    
    # Function to update the piston and gas color in each frame
    def update(frame):
        # Update the piston height
        new_height = piston_height + frame * (distance_moved / 100)
        piston.set_y(new_height)
        
        # Update the gas color based on the temperature change
        temp_fraction = frame / 100  # Fraction of total temperature change
        current_temp = T_initial + temp_fraction * (T_final - T_initial)
        
        # Interpolate color from blue (cold) to red (hot)
        gas_color = (1, 0, 0, temp_fraction)  # RGBA where A controls transparency
        cylinder.set_facecolor(gas_color)
        
        # Update the temperature text
        temp_text.set_text(f"Temp: {current_temp:.2f} K")
    
    # Create the animation
    ani = FuncAnimation(fig, update, frames=np.linspace(0, 100, 100), repeat=False)
    
    # Show the animation
    plt.show()

# Example inputs (students can modify these values)
system_gas = input("Enter the type of gas (air, helium, argon, oxygen, nitrogen): ").lower()
T_environment = float(input("Enter the environmental temperature (K): "))
V_initial = float(input("Enter the initial volume of the gas (m³): "))
Q_in = float(input("Enter the heat added to the system (J): "))
P = float(input("Enter the pressure of the gas (Pa): "))

# Perform the simulation
T_final, V_final, distance_moved_mm, delta_U = first_law_thermodynamics(system_gas, T_environment, V_initial, Q_in, P)

# Run the animation
animate_piston(T_environment, T_final, distance_moved_mm)
