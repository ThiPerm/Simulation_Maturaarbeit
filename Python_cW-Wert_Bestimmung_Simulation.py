import numpy as np

def find_cWs(cW_start, cW_end, cW_accuracy, init_pos_2d, init_vel_2d, x_target_arr):
    if len(x_target_arr) == len(init_pos_2d) == len(init_vel_2d):
        cWs = []
        for index, element in enumerate(x_target_arr):
            init_position = init_pos_2d[index]
            init_velocity = init_vel_2d[index]
            x_target = x_target_arr[index]
            cW = find_cW(cW_start, cW_end, cW_accuracy, x_target)
            cWs.append(cW)
            print(index)
        np.savetxt("cWs.txt", cWs)
    else:
        print(f"Die Listen sind haben nicht die gleiche Länge: \n"
              f"x_target_arr: {len(x_target_arr)} Elemente \n"
              f"init_pos_2d: {len(init_pos_2d)} Elemente \n"
              f"init_vel_2d: {len(init_vel_2d)} Elemente \n")


def closest_value(lst, target):
    closest = lst[0]  # erstes Element als Nächstes definieren
    # Alle anderen Elemente werden mit dem Nächsten verglichen und falls kleiner eintauschen
    closest_ind = len(lst)-1
    for ind, element in enumerate(lst):
        if abs(element - target) < abs(closest - target):
            closest = element  # Update closest if the current value is closer to the target
            closest_ind = ind
    return closest_ind


def find_cW(cW_start, cW_end, cW_accuracy, x_target):
    cW_current = cW_start
    sim_x_ends = []
    cW_currents = []
    while cW_current < cW_end:
        cW_currents.append(cW_current)
        pos, vel, time = simulate(init_position, init_velocity, time_step, total_time, ang_velocity, cW_current)
        pos_x, pos_y, pos_z = np.transpose(pos)
        sim_x_end = pos_x[closest_value(pos_z, 0)]
        sim_x_ends.append(sim_x_end)
        cW_current = cW_current + cW_accuracy
    cW_closest = cW_currents[closest_value(sim_x_ends, x_target)]
    print(cW_closest)
    return cW_closest


def simulate(init_position, init_velocity, ang_velocity, magnus_coefficient):
    positions = []  # Liste, in der die Position jeden Zeitschritt gespeichert wird
    velocities = []  # Liste, in der die Geschwindigkeit jeden Zeitschritt gespeichert wird
    times = []
    mass = 0.0577

    position = np.array(init_position, dtype=float)
    velocity = np.array(init_velocity, dtype=float)
    impulse = np.array([init_velocity[0] * mass, init_velocity[1] * mass, init_velocity[2] * mass])
    for time in np.arange(0, total_time, time_step):
        # Lokale Variablen
        g = 9.81
        rho_air = 1.15
        diameter = 0.065
        area = diameter * diameter / 4 * np.pi

        net_velocity = (velocity[0] ** 2 + velocity[1] ** 2 + velocity[2] ** 2) ** (1 / 2)
        # Gravitationskraft berechnen
        gravitational_force = np.array([0, 0, -g * mass])

        # Luftwiderstand berechnen
        total_air_resistance = -1 / 2 * cW * rho_air * area * net_velocity
        air_resistance = np.array([total_air_resistance * velocity[0], total_air_resistance * velocity[1],
                                   total_air_resistance * velocity[2]])

        # Magnus kraft berechnen
        magnus_force = -magnus_coefficient * 4 / 3 * np.pi * rho_air * (diameter / 2) ** 3 * np.cross(velocity,
                                                                                                      ang_velocity)

        # resultierende Kraft berechnen
        net_force = gravitational_force + air_resistance + magnus_force

        # Veränderung des Impulses durch den Vektor der Kraft je Zeitschritt
        impulse = impulse + (net_force * time_step)

        # Veränderung der Geschwindigkeit durch den Impuls
        velocity = impulse / mass

        # Veränderung der Position durch die Geschwindigkeit pro Zeitschritt
        position += velocity * time_step

        positions.append(position.copy())
        velocities.append(velocity.copy())
        times.append(time.copy())

        xx, yy, zz = np.transpose(position)
        if zz < -0.5:
            return np.array(positions), np.array(velocities), np.array(times)

    return np.array(positions), np.array(velocities), np.array(times)


# Startkondition bei mehreren abfragen gleichzeitig
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
x_target_arr = np.array([
    9.765,
    9.806,
    9.683,
    9.601,
    9.909,
    9.641,
    9.453,
    9.82,
    10.311,
    9.555,
    9.719,
    9.506,
    9.213,
    9.227,
    9.085,
    9.053,
    8.768,
    8.72,
    9.116,
    9.243,
    9.18,
    8.895,
    8.895,
    9.077,
    11.755,
    11.616,
    11.968,
    11.899,
    11.749,
    11.094,
    11.081,
    11.464,
    11.098,
    10.964,
    10.724,
    11.074,
])

# Startkonditionen für eine abfrage
init_position = [-11.88, 0, 0.55]
init_velocity = [14.04, 0, 9.61]
ang_velocity = [0, 0, 0]
x_target = 9.765  # erzieltes X (durch Messung / Video)

# Allgemeine Startkonditionen
time_step = 0.001
total_time = 10  # unwichtig, da die Simulation bei zu tiefer Z-Koordinate abbricht

# Definition der zu überprüfenden cW-Werte
cW_start = 0.27  # Startwert für cW
cW_end = 0.61  # Endwert für cW
cW_accuracy = 0.001  # Vergrösserungsschritte von cW

# Eine der beiden folgenden Funktionen wählen (die andere "auskommentieren")
# Funktion, um den exakten cW-Wert für einen Schuss zu finden
find_cW(cW_start, cW_end, cW_accuracy, x_target)

# Funktion, um die exakten cW-Werte für alle Schüsse aus dem Array zu finden
find_cWs(cW_start, cW_end, cW_accuracy, init_pos_2d, init_vel_2d, x_target_arr)
