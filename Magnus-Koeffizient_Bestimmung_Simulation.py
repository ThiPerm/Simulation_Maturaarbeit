import numpy as np

def find_magnus_coefficients():
    if len(x_target_arr) == len(init_pos_2d) == len(init_vel_2d) == len(ang_vel_2d):
        mcs = []
        for index, element in enumerate(x_target_arr):
            init_position_a = init_pos_2d[index]
            init_velocity_a = init_vel_2d[index]
            ang_velocity_a = ang_vel_2d[index]
            x_target_a = x_target_arr[index]
            mc = find_magnus_coefficient(init_position_a, init_velocity_a, ang_velocity_a, x_target_a)
            mcs.append(mc)
        np.savetxt("Magnus-Koeffizienten.txt", mcs)
    else:
        print(f"Die Listen sind haben nicht die gleiche Länge: \n"
              f"x_target_arr: {len(x_target_arr)} Elemente \n"
              f"init_pos_2d: {len(init_pos_2d)} Elemente \n"
              f"init_vel_2d: {len(init_vel_2d)} Elemente \n"
              f"ang_vel_2d: {len(ang_velocity_b)} Elemente \n")



def closest_value(lst, target):
    closest = lst[0]  # erstes Element als Nächstes definieren
    # Alle anderen Elemente werden mit dem Nächsten verglichen und falls kleiner eintauschen
    closest_ind = len(lst)-1
    for ind, element in enumerate(lst):
        if abs(element - target) < abs(closest - target):
            closest = element  # Update closest if the current value is closer to the target
            closest_ind = ind
    return closest_ind


def find_magnus_coefficient(init_position, init_velocity, ang_velocity, x_target):
    mc_current = mc_start
    sim_x_ends = []
    mc_currents = []
    while mc_current < mc_end:
        mc_currents.append(mc_current)
        pos, vel, time = simulate(init_position, init_velocity, ang_velocity, mc_current)
        pos_x, pos_y, pos_z = np.transpose(pos)
        sim_x_end = pos_x[closest_value(pos_z, 0)]
        sim_x_ends.append(sim_x_end)
        mc_current = mc_current + mc_accuracy
    mc_closest = mc_currents[closest_value(sim_x_ends, x_target)]
    print(f"Der beste Wert für den Koeffizienten ist {mc_closest}")
    return mc_closest


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


# Startkonditionen für eine abfrage
init_position_b = [-11.88, 0, 0.55]
init_velocity_b = [14.04, 0, 9.61]
ang_velocity_b = [0, 10, 0]
x_target_b = 9.765   # erzieltes X (durch Messung / Video)

# Startkondition bei mehreren abfragen gleichzeitig
init_vel_2d = np.array([[14.04, 0, 9.61], [14.30, 0, 9.48]])
init_pos_2d = np.array([[-11.88, 0, 0.55], [-11.88, 0, 0.55]])
ang_vel_2d = np.array([[0, 0, 0], [0, 0, 0]])
x_target_arr = np.array([9.765, 7.892])

# Allgemeine Startkonditionen
time_step = 0.001
total_time = 10
cW = 0.472

# Definition der zu überprüfenden Koeffizienten
mc_start = 0  # Startwert für Magnus-Koeffizient
mc_end = 1  # Endwert für Magnus-Koeffizient
mc_accuracy = 0.01  # Vergrösserungsschritte vom Magnus-Koeffizient

# Eine der beiden folgenden Funktionen wählen (die andere "auskommentieren")
# Funktion, um den exakten Koeffizienten für einen Schuss zu finden
find_magnus_coefficient(init_position_b, init_velocity_b, ang_velocity_b, x_target_b)

# Funktion, um die exakten Koeffizienten für alle Schüsse aus dem Array zu finden
find_magnus_coefficients()
