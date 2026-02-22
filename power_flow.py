import numpy as np

# System parameters
# Line data [From Bus, To Bus, Length (km)]
lines = np.array([
    [1, 2, 20],
    [1, 4, 40],
    [1, 3, 10],
    [2, 4, 10],
    [3, 4, 20]
])

# Constants
f = 50  # frequency in Hz
w = 2 * np.pi * f

# Line impedance parameters
R = 109e-3  # Resistance in ohm/km
X = 1.2e-3  # Reactance in H/km
C = 9.5e-9  # Capacitance in F/km

# Bus data
# [Bus, Voltage Mag (p.u), Angle (degree), P (MW), Q (MVAR), P-Q spec type]
buses = np.array([
    [1, 1.05, 0.0,    0,  0],    # Slack bus
    [2, 1.02, 0.0, -300,  0],    # PV bus
    [3, 1.0,  0.0,  200, 50],    # PQ bus
    [4, 1.0,  0.0,  100, 20]     # PQ bus
])

# Number of buses and lines
num_buses = len(buses)
num_lines = len(lines)

# Create the admittance matrix
Ybus = np.zeros((num_buses, num_buses), dtype=complex)

for line in lines:
    fb, tb, length = line
    fb = int(fb) - 1
    tb = int(tb) - 1
    z = complex(R * length, X * length * w)
    y = 1 / z
    c = C * length * w

    Ybus[fb, fb] += y + 1j * c / 2
    Ybus[tb, tb] += y + 1j * c / 2
    Ybus[fb, tb] -= y
    Ybus[tb, fb] -= y

real_admittance = np.real(Ybus)
imaginary_admittance = np.imag(Ybus)

# Power equation
def power_eq(is_active_power, subindex, nodal_voltage, phase):
    subindex -= 1
    power = 0.
    sum = 0.
    if is_active_power == True:
        for k in range (len(nodal_voltage)):
            delta_theta = phase[subindex] - phase[k]
            sum += nodal_voltage[k]*(real_admittance[subindex, k] * np.cos(delta_theta) + imaginary_admittance[subindex, k] * np.sin(delta_theta)) 
    else:
        for k in range (len(nodal_voltage)):
            delta_theta = phase[subindex] - phase[k]
            sum += nodal_voltage[k] * (real_admittance[subindex, k] * np.sin(delta_theta) - imaginary_admittance[subindex, k] * np.cos(delta_theta))
    power = nodal_voltage[subindex] * sum
    return power

# Partial phase calculation
def partial_phase(is_active_power, function, subindex, nodal_voltage, phase):
    function -= 1
    subindex -= 1
    delta_theta = phase[function] - phase[subindex]
    value = 0.0

    if is_active_power:
        if function == subindex:
            value = -(power_eq(False, subindex + 1, nodal_voltage, phase) + nodal_voltage[subindex]**2 * imaginary_admittance[subindex, subindex])
        else:
            value = nodal_voltage[function] * nodal_voltage[subindex] * (
                real_admittance[function, subindex] * np.sin(delta_theta) -
                imaginary_admittance[function, subindex] * np.cos(delta_theta)
            )
    else:
        if function == subindex:
            value = power_eq(True, subindex + 1, nodal_voltage, phase) - nodal_voltage[subindex]**2 * real_admittance[subindex, subindex]
        else:
            value = -nodal_voltage[function] * nodal_voltage[subindex] * (
                real_admittance[function, subindex] * np.cos(delta_theta) +
                imaginary_admittance[function, subindex] * np.sin(delta_theta)
            )
    return value

# Partial voltage calculation
def partial_voltage(is_active_power, function, subindex, nodal_voltage, phase):
    function -= 1
    subindex -= 1
    delta_theta = phase[function] - phase[subindex]
    value = 0.0

    if is_active_power:
        if function == subindex:
            value = power_eq(True, subindex + 1, nodal_voltage, phase) / nodal_voltage[subindex] + nodal_voltage[subindex] * real_admittance[subindex, subindex]
        else:
            value = nodal_voltage[function] * (
                real_admittance[function, subindex] * np.cos(delta_theta) +
                imaginary_admittance[function, subindex] * np.sin(delta_theta)
            )
    else:
        if function == subindex:
            value = power_eq(False, subindex + 1, nodal_voltage, phase) / nodal_voltage[subindex] - nodal_voltage[subindex] * imaginary_admittance[subindex, subindex]
        else:
            value = nodal_voltage[function] * (
                real_admittance[function, subindex] * np.sin(delta_theta) -
                imaginary_admittance[function, subindex] * np.cos(delta_theta)
            )
    return value

# Solver function
def solver(power_target, nodal_voltage, phase):
    size = len(power_target)
    power_n = np.zeros(size)
    power_n[0] = power_eq(True, 2, nodal_voltage, phase)
    power_n[1] = power_eq(True, 3, nodal_voltage, phase)
    power_n[2] = power_eq(True, 4, nodal_voltage, phase)
    power_n[3] = power_eq(False, 3, nodal_voltage, phase)
    power_n[4] = power_eq(False, 4, nodal_voltage, phase)

    f_nonlinear_n = power_n - power_target
    j_matrix = np.zeros((size, size))

    for i in range(3):
        for j in range(3):
            j_matrix[i, j] = partial_phase(True, i + 2, j + 2, nodal_voltage, phase)
        j_matrix[i, 3] = partial_phase(False, i + 2, 3, nodal_voltage, phase)
        j_matrix[i, 4] = partial_phase(False, i + 2, 4, nodal_voltage, phase)

    for i in range(3, 5):
        for j in range(3):
            j_matrix[i, j] = partial_voltage(True, i, j + 2, nodal_voltage, phase)
        j_matrix[i, 3] = partial_voltage(False, i, 3, nodal_voltage, phase)
        j_matrix[i, 4] = partial_voltage(False, i, 4, nodal_voltage, phase)

    delta_x = np.linalg.solve(np.transpose(j_matrix), -f_nonlinear_n)
    return delta_x, power_n

# Map function to ensure value is within 2*pi
def map_2pi(x):
    return np.sign(x) * (x % (2 * np.pi))

def main(turbine_power_generation=300, compressor_power=0):

    power_target=np.array([turbine_power_generation, compressor_power - 200, -100.0, -50.0, -20.0])

    # Operation conditions at nodes
    nodal_voltage = np.array([1.05, 1.02, -0.0, -0.0]) * 110  # kV
    phase = np.array([0.0, -0.0, -0.0, -0.0])  # rad

    # Initialization for the NR-solver
    tolerance = 1e-6
    delta_x = 10 * np.ones(len(power_target))
    n_iter = 0
    nodal_voltage[2] = 1.0 * 110
    nodal_voltage[3] = 1.0 * 110

    while np.linalg.norm(delta_x) > tolerance:
        delta_x, power_n = solver(power_target, nodal_voltage, phase)
        phase[1:4] += delta_x[0:3]
        nodal_voltage[2:4] += delta_x[3:5]
        n_iter += 1

    line_currents = []
    for line in lines:
        fb, tb, length = line
        fb = int(fb) - 1
        tb = int(tb) - 1
        vi = nodal_voltage[fb] * np.exp(1j * phase[fb])
        vj = nodal_voltage[tb] * np.exp(1j * phase[tb])
        y = Ybus[fb, tb]
        current = (vi - vj) * y
        line_currents.append(current)

    print(f'Simulation converges in {n_iter} iterations!')
    print(f"Norm of delta_x = {np.linalg.norm(delta_x)}")
    print(f"Phases = {phase}")
    print(f"Nodal Voltages = {nodal_voltage}")
    print(f"Line currents:\n L1: {line_currents[0]:.2f} L2: {line_currents[1]:.2f} L3: {line_currents[2]:.2f} L4: {line_currents[3]:.2f} L5: {line_currents[4]:.2f}")

if __name__ == '__main__':
    main()