import numpy as np
from scipy.constants import bar, atm, zero_Celsius, pi, R, g

# Define constant variables
T_ST = 15 + zero_Celsius  # standard temperature (15 degree Celcius in K)
P_ST = 1 * atm  # standard pressure in Pa

# Natural Gas (mixture) properties
RD = 0.58 
Z = 0.88
RHO = 0.7  # kg/sm^3
CP = 2142  # kJ/(kg*K)

# Average molar mass of air
M_AIR = 28.97  # g/mol

# Determine the C1 and C2
C1 = pi * (R * 1000 / 16 / M_AIR)**0.5 # R is in J/(mol*K), we need it in J/(kmol*K)
C2 = 2 * g * M_AIR / (R * 1000)

# Function to determine the flow direction
def flow_direction(p1, p2, e):
    if (p1**2 - p2**2 - e) >= 0:
        return 1
    else:
        return -1

# Function to calculate the height difference effect
def height_effect(h1, h2, p1, p2, t1, t2, d=RD, Z=Z):
    p_avg = 2/3 * ((p1 + p2) - (p1 * p2) / (p1 + p2))
    t_avg = 1/2 * (t1 + t2)
    return C2 * d * (h2 - h1) * p_avg**2 / (t_avg*Z)

# Create a helper function to calculate C_ij of a pipe which won't change during one iteration
def c_pipe(t1, t2, L, D=0.5, f=0.01, eta=0.85, Z=Z, d=RD):
    t_avg = (t1 + t2)/2
    return C1 * T_ST/P_ST * (D**2.5) * eta * (1 / (L*t_avg*f*d*Z))**0.5


def volumetric_flow_rate(p1, p2, t1, t2, h1, h2, L, D=0.5, d=RD, Z=Z, f=0.01, eta=0.85):
    p_avg = 2/3 * ((p1+p2) - (p1*p2) / (p1 + p2))
    t_avg = 1/2 * (t1 + t2)
    e = height_effect(h1, h2, p1, p2, t1, t2)
    f_direction = flow_direction(p1, p2, e)
    return f_direction * c_pipe(t1, t2, L, D) * (abs(p1**2 - p2**2 - e))**0.5

def pressure_compression_station(p1, r=1.1):
    return p1 * r

def power_compression_station(m_flow, p1, p2, t1, n=1.3, rho=RHO, cp=CP):
    return (cp * t1) * ((p2/p1)**((n-1)/n) - 1) * m_flow

def main(f_target=np.array([60., 14.0])):
    # Operation conditions at nodes
    p1 = 50 * bar # 50 bar
    t1 = 300
    t2 = 288.15
    t3 = 288.15
    t4 = 288.15
    
    # Node heights
    h1 = 30
    h2 = 30
    h3 = 30
    h4 = 0

    # Pipe lengths
    l1 = 200e3
    l2 = 400e3
    l3 = 100e3

    # Pipe diameters
    # d1 = 0.5
    # d2 = 0.5
    # d3 = 0.5

    c12 = c_pipe(t1, t2, l1)
    c34 = c_pipe(t3, t4, l2)
    c14 = c_pipe(t1, t4, l3)
    
    # Initialization for the NR-solver
    err = 1
    tol = 0.01
    p2 = 0.98 * p1
    p3 = pressure_compression_station(p2)
    p4 = 0.98 * p3
    x0 = np.array([p2, p4])
    p = x0
    n_iter = 0
    while err > tol:
        p2 = p[0]
        p3 = pressure_compression_station(p2)
        p4 = p[1]

        # F(X)
        q12 = volumetric_flow_rate(p1, p2, t1, t2, h1, h2, l1)
        q43 = volumetric_flow_rate(p4, p3, t4, t3, h4, h3, l2)
        q14 = volumetric_flow_rate(p1, p4, t1, t4, h1, h4, l3)
        q34 = volumetric_flow_rate(p3, p4, t3, t4, h3, h4, l2)

        # Jacobian matrix
        e12 = height_effect(h1, h2, p1, p2, t1, t2)
        e34 = height_effect(h3, h4, p3, p4, t3, t4)
        e14 = height_effect(h1, h4, p1, p4, t1, t4)
        dq12_dp2 = c12 * p2 / (abs(p1**2 - p2**2 - e12))**0.5
        dq43_dp3 = c34 * p3 / (abs(p3**2 - p4**2 - e34))**0.5
        dq43_dp4 = c34 * (-p4) / (abs(p3**2 - p4**2 - e34))**0.5
        dq34_dp3 = c34 * (-p3) / (abs(p3**2 - p4**2 - e34))**0.5
        dq14_dp4 = c14 * p4 / (abs(p1**2 - p4**2 - e14))**0.5
        dq34_dp4 = c34 * p4 / (abs(p3**2 - p4**2 - e34))**0.5 
        j_mat = np.zeros((2, 2))
        j_mat[0][0]= dq12_dp2 + dq43_dp3
        j_mat[0][1] = dq43_dp4
        j_mat[1][0] = dq34_dp3
        j_mat[1][1] = dq14_dp4 + dq34_dp4

        # Balance on node 2 and 3
        d23 = q12 + q43
        # Balance on node 4
        d4 = q14 + q34

        # Update the unkonwn variables X
        j_mat_inv = np.linalg.inv(j_mat)
        f = np.array([d23, d4])
        delta_f = f_target - f
        p = p - np.dot(j_mat_inv, delta_f)
        err = max(abs(delta_f))
        n_iter += 1

    print(f'Simulation converges in {n_iter} iterations!')
    print(f'P1: {p1}, P2: {p2}, P3: {p3}, P4: {p4}')
    print(f'Flow pipe 1: {q12}, flow pipe 2: {q34} flow pipe 3: {q14}')
          

    m_flow_compressor = (q12 - f_target[0]*0.5) * RHO
    power_consumption_compressor = power_compression_station(m_flow_compressor, p2, p3, t2) * (1e-6) # in MW
    
    return power_consumption_compressor


if __name__ == '__main__':
    main()