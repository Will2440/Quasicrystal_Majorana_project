# import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib.animation as animation

# # Define the Mexican hat potential
# def mexican_hat_potential(x, y, b, a=1):
#     return a * ((x**2 + y**2) - b**2)**2

# # Create a meshgrid for 2D plotting
# x = np.linspace(-3, 3, 400)
# y = np.linspace(-3, 3, 400)
# X, Y = np.meshgrid(x, y)

# # Parameters for animation
# b_values = np.linspace(0, 2, 50)  # Varying b from 0 (trivial) to 2 (non-trivial)

# # Create figure for the animation
# fig, ax = plt.subplots(subplot_kw={"projection": "3d"}, figsize=(8, 6))
# ax.set_zlim(0, 10)

# def update(frame):
#     ax.clear()  # Clear the plot for updating
#     Z = mexican_hat_potential(X, Y, b_values[frame])
#     ax.plot_surface(X, Y, Z, cmap="viridis", edgecolor="none")
#     ax.set_title(f'Mexican Hat Potential with b = {b_values[frame]:.2f}')
#     ax.set_zlim(0, 10)  # Ensure consistent z-axis limits

# # Create the animation
# ani = animation.FuncAnimation(fig, update, frames=len(b_values), interval=100)

# # Save the animation as a GIF
# ani.save('mexican_hat_potential.gif', writer='pillow')

# # plt.show()



# import numpy as np
# import matplotlib.pyplot as plt

# def mexican_hat_potential(x, y, b, a=1):
#     return a * ((x**2 + y**2) - b**2)**2

# x = np.linspace(-3, 3, 400)
# y = np.linspace(-3, 3, 400)
# X, Y = np.meshgrid(x, y)

# # # b = 100.0
# b_values = np.linspace(2, -2, 5)

# for i, value in enumerate(b_values):
#     Z = mexican_hat_potential(X, Y, value)

#     fig, ax = plt.subplots(subplot_kw={"projection": "3d"}, figsize=(8, 6))
#     ax.set_zlim(0, 20)
#     surface = ax.plot_surface(X, Y, Z, cmap="viridis", edgecolor="none", alpha=0.7)

#     ax.set_xlabel('X axis')
#     ax.set_ylabel('Y axis')
#     ax.set_zlabel('Potential V(x, y)')
#     ax.set_title(f'Mexican Hat Potential with b = {value}')

# plt.show()










import numpy as np
import matplotlib.pyplot as plt
import imageio
import os

def mexican_hat_potential(x, y, b, a=1):
    return a * ((x**2 + y**2) - b**2)**2

def plot_mexican_hat_gif(b_values, gif_name='mexican_hat.gif'):
    """
    Creates a GIF that cycles through the Mexican hat potential for each b value.

    Args:
        - b_values: List or array of b values for the Mexican hat potential.
        - gif_name: Str: Name of the output GIF file.
    """
    # Create a meshgrid for the X and Y values
    x = np.linspace(-3, 3, 400)
    y = np.linspace(-3, 3, 400)
    X, Y = np.meshgrid(x, y)

    # Create a list to hold filenames of individual plot images
    filenames = []

    for i, value in enumerate(b_values):
        Z = mexican_hat_potential(X, Y, value)

        # Create a new figure for each subplot
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"}, figsize=(8, 6))
        ax.set_zlim(0, 20)
        surface = ax.plot_surface(X, Y, Z, cmap="viridis", edgecolor="none", alpha=0.7)

        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('Potential V(x, y)')
        ax.set_title(f'Mexican Hat Potential with b = {value:.2f}')

        # Save the current figure to a temporary file
        filename = f'frame_{i}.png'
        plt.savefig(filename)
        plt.close(fig)  # Close the figure to avoid display
        filenames.append(filename)

    # Create a GIF from the saved frames
    with imageio.get_writer(gif_name, mode='I', duration=0.5) as writer:
        for filename in filenames:
            image = imageio.imread(filename)
            writer.append_data(image)

    # Optionally, clean up the temporary files
    for filename in filenames:
        os.remove(filename)

    print(f'GIF saved as {gif_name}')

# Example usage
b_values = np.linspace(2, -2, 5)
plot_mexican_hat_gif(b_values)
