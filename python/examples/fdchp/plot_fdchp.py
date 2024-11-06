import matplotlib.pyplot as plt
import numpy as np


def plot_x_velocity_vs_wind_speed(speeds, x_velocities):
    fig, ax = plt.subplots()
    ax.scatter(speeds, x_velocities, marker='o', alpha=0.5)
    ax.set_xlabel('U (m/s)')
    ax.set_xlim(0, 20)
    ax.set_ylabel('<uw> (m/s)')
    ax.set_ylim(-1.0, 0.0)
    plt.show()
    
    
    
def plot_wave_height_vs_wind_speed(speeds, wave_heights):
    fig, ax = plt.subplots()
    ax.scatter(speeds, wave_heights, marker='o', alpha=0.5)
    ax.set_xlabel('U (m/s)')
    ax.set_xlim(0, 15)
    ax.set_ylabel('$\sigma_H$ (m)')
    ax.set_ylim(0, 4)
    plt.show()

