import matplotlib.pyplot as plt
import numpy as np


def plot_x_velocity_vs_wind_speed(speeds, x_velocities, x_min = 0, x_max = 20, y_min = -1.0, y_max = 0.0):
    fig, ax = plt.subplots()
    ax.scatter(speeds, x_velocities, marker='o', alpha=0.5)
    ax.set_xlabel('U (m/s)')
    ax.set_xlim(x_min, x_max)
    ax.set_ylabel('<uw> (m/s)')
    ax.set_ylim(y_min, y_max)
    plt.show()
    
    
    
def plot_wave_height_vs_wind_speed(speeds, wave_heights,  x_min = 0, x_max = 15, y_min = 0.0, y_max = 4.0):
    fig, ax = plt.subplots()
    ax.scatter(speeds, wave_heights, marker='o', alpha=0.5)
    ax.set_xlabel('U (m/s)')
    ax.set_xlim(x_min, x_max)
    ax.set_ylabel('$\sigma_H$ (m)')
    ax.set_ylim(y_min, y_max)
    plt.show()

