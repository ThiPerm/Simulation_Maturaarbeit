import matplotlib.pyplot as plt
import numpy as np


def closest_value(lst, target):
    closest = lst[0]  # erstes Element als Nächstes definieren
    # Alle anderen Elemente werden mit dem Nächsten verglichen und falls kleiner eintauschen
    for ind, element in enumerate(lst):
        if abs(element - target) < abs(closest - target):
            closest = element  # Überschreibe closest, wenn das element näher ist
            closest_ind = ind
    return closest_ind


def simulate():
    positions = []  # Liste, in der die Position jeden Zeitschritt gespeichert wird
    velocities = []  # Liste, in der die Geschwindigkeit jeden Zeitschritt gespeichert wird
    times = [] # Liste, in der die Zeit gespeichert wird, damit der Index mit anderen Listen verglichen werden kann

    position = np.array(init_position, dtype=float)
    velocity = np.array(init_velocity, dtype=float)
    impulse = np.array([init_velocity[0] * mass, init_velocity[1] * mass, init_velocity[2] * mass])

    for time in np.arange(0, total_time, time_step):
        # Lokale Variablen
        g = 9.81
        rho_air = 1.15
        diameter = 0.065
        area = diameter * diameter / 4 * np.pi
        magnus_coefficient = 1

        # Gravitationskraft berechnen
        gravitational_force = np.array([0, 0, -g * mass])

        # Luftwiderstand berechnen
        net_velocity_float = (velocity[0] ** 2 + velocity[1] ** 2 + velocity[2] ** 2) ** (1 / 2)
        total_air_resistance = -1 / 2 * cW * rho_air * area * net_velocity_float
        air_resistance = np.array([total_air_resistance * velocity[0], total_air_resistance * velocity[1],
                                   total_air_resistance * velocity[2]])

        # Magnus kraft berechnen
        magnus_force = -magnus_coefficient * 4 / 3 * np.pi * rho_air * (diameter / 2) ** 3 * np.cross(velocity,
                                                                                                      ang_velocity)

        # resultierende Kraft berechnen
        net_force = gravitational_force + air_resistance + magnus_force
        # Impuls
        impulse += net_force * time_step

        # Veränderung der Geschwindigkeit durch den Impuls
        velocity = impulse / mass

        # Veränderung der Position durch die Geschwindigkeit pro Zeitschritt
        position += velocity * time_step

        positions.append(position.copy())
        velocities.append(velocity.copy())
        times.append(time.copy())

        # Stoppt Simulation, wenn z (zz) kleiner als -0.5 ist
        xx, yy, zz = np.transpose(position)
        if zz < -0.5:
            return np.array(positions), np.array(velocities), np.array(times)

    return np.array(positions), np.array(velocities), np.array(times)

def results():
    # Simulate the motion
    positions, velocities, times = simulate()

    x, y, z = np.transpose(positions)
    end_index = closest_value(z, 0)
    x_end = x[end_index]
    time_end = times[end_index]
    print(f"The Ball touches the ground at x = {x_end}m after t = {time_end}s")
    return positions, velocities, times, x, y, z


def results_arr():
    x_ends = []
    times_end = []
    for ind, element in enumerate(init_pos_2d):
        init_position = init_pos_2d[ind]
        init_velocity = init_vel_2d[ind]
        positions, velocities, times = simulate()

        x, y, z = np.transpose(positions)
        end_index = closest_value(z, 0)
        x_end = x[end_index]
        x_ends.append(x_end.copy())
        time_end = times[end_index]
        times_end.append(time_end.copy())
    np.savetxt("end_positions.txt", x_ends)
    np.savetxt("end_times.txt", times_end)



# Anfangskonditionen
mass = 0.0577
cW = 0.472
init_position = [-11.88, 0, 0.58]
init_velocity = [17.74, 0, 7.41]
init_vel_2d = np.array([
    [14.04, 0, 9.61],
    [14.30, 0, 9.48],
    [13.73, 0, 9.60],
    [13.89, 0, 9.57],
    [14.12, 0, 9.80],
    [14.01, 0, 9.58],
    [13.84, 0, 9.40],
    [13.83, 0, 9.54],
    [14.05, 0, 9.89],
    [14.07, 0, 9.26],
    [13.92, 0, 9.53],
    [13.94, 0, 9.37],
    [15.89, 0, 7.77],
    [15.77, 0, 7.74],
    [15.64, 0, 7.68],
    [15.59, 0, 7.56],
    [15.49, 0, 7.61],
    [15.51, 0, 7.46],
    [15.70, 0, 7.88],
    [15.82, 0, 7.91],
    [15.96, 0, 7.73],
    [15.53, 0, 7.58],
    [15.72, 0, 7.62],
    [15.78, 0, 7.77],
    [17.74, 0, 7.41],
    [17.79, 0, 7.67],
    [17.61, 0, 7.54],
    [17.83, 0, 7.45],
    [18.00, 0, 7.60],
    [17.79, 0, 7.53],
    [17.39, 0, 7.49],
    [17.50, 0, 7.66],
    [17.66, 0, 7.42],
    [17.60, 0, 7.37],
    [17.56, 0, 7.53],
    [17.70, 0, 7.63],
])
init_pos_2d = np.array([
    [-11.88, 0, 0.55],
    [-11.88, 0, 0.55],
    [-11.88, 0, 0.55],
    [-11.88, 0, 0.55],
    [-11.88, 0, 0.55],
    [-11.88, 0, 0.55],
    [-11.88, 0, 0.55],
    [-11.88, 0, 0.55],
    [-11.88, 0, 0.55],
    [-11.88, 0, 0.55],
    [-11.88, 0, 0.55],
    [-11.88, 0, 0.55],
    [-11.88, 0, 0.57],
    [-11.88, 0, 0.57],
    [-11.88, 0, 0.57],
    [-11.88, 0, 0.57],
    [-11.88, 0, 0.57],
    [-11.88, 0, 0.57],
    [-11.88, 0, 0.57],
    [-11.88, 0, 0.57],
    [-11.88, 0, 0.57],
    [-11.88, 0, 0.57],
    [-11.88, 0, 0.57],
    [-11.88, 0, 0.57],
    [-11.88, 0, 0.58],
    [-11.88, 0, 0.58],
    [-11.88, 0, 0.58],
    [-11.88, 0, 0.58],
    [-11.88, 0, 0.58],
    [-11.88, 0, 0.58],
    [-11.88, 0, 0.58],
    [-11.88, 0, 0.58],
    [-11.88, 0, 0.58],
    [-11.88, 0, 0.58],
    [-11.88, 0, 0.58],
    [-11.88, 0, 0.58],
])
ang_velocity = [0, 0, 0]
time_step = 0.001
total_time = 3

positions, velocities, times, x, y, z = results()
#results_arr()



# Kreiert einen 3d Graphen und die gleichen Graphen wie in der Vensim Simulation mithilfe von matplotlib
# von ChatGPT 3.5 geschrieben

#plt.figure(figsize=(12, 8))

# 3D-Diagramm initialisieren
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Punktdiagramm im 3D-Raum erstellen
ax.scatter(x, y, z, c=z, cmap='viridis')

# Beschriftungen hinzufügen
ax.set_xlabel('X-Achse')
ax.set_ylabel('Y-Achse')
ax.set_zlabel('Z-Achse')

#plt.subplot(3, 1, 1)
#plt.plot(x, y)
#plt.title("XY")

#plt.subplot(3, 1, 2)
#plt.plot(x, z)
#plt.title("XZ")

#plt.subplot(3, 1, 3)
#plt.plot(y, z)
#plt.title("YZ")

#plt.tight_layout()
#plt.show()
