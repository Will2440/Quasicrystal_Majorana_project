# from PIL import Image
# import os

# def create_gif(image_paths, output_path, duration=500, loop=0, deletion=False):
#     """
#     Creates a GIF from a list of PNG images.

#     Parameters:
#     - image_paths: List of file paths to the PNG images.
#     - output_path: Path to save the generated GIF (including .gif extension).
#     - duration: Duration for each frame in milliseconds. Default is 500 ms.
#     - loop: Number of loops for the GIF. Default is 0 (infinite loop).

#     Returns:
#     - None
#     """
#     if not image_paths:
#         raise ValueError("The list of image paths is empty.")
    
#     images = [Image.open(img) for img in image_paths]
    
#     images = [img.convert("RGB") if img.mode != "RGB" else img for img in images]
    
#     # Save as gif
#     images[0].save(
#         output_path,
#         save_all=True,
#         append_images=images[1:],  # Add the rest of the images
#         duration=duration,        # Duration per frame
#         loop=loop                 # Number of loops
#     )
#     print(f"GIF saved at {output_path}")

#     if deletion:
#         for img_path in image_paths:
#             try:
#                 os.remove(img_path)
#                 print(f"Deleted: {img_path}")
#             except OSError as e:
#                 print(f"Error deleting {img_path}: {e}")


# # Example usage
# if __name__ == "__main__":

#     png_files = []
#     for i in range(401):
#         # png_files.append(f"/Users/Will/Documents/FINAL_PROJECT/simulations/mp_sitewise_{i+1}.png")
#         png_files.append(f"/Users/Will/Documents/FINAL_PROJECT/simulations/Series_FFT_GQC_N377_rho{i+1}.png")
    
#     output_gif = "/Users/Will/Documents/FINAL_PROJECT/simulations/images_and_graphs/dimerisation/Series_FFT_GQC_N377_rho.gif"
    
#     # Ensure all files exist
#     png_files = [f for f in png_files if os.path.exists(f)]
    
#     create_gif(png_files, output_gif, duration=500, loop=10, deletion=False)


from PIL import Image
import os
import numpy as np

def create_gif(image_paths, output_path, duration=500, loop=0, deletion=False, bounce=False):
    """
    Creates a GIF from a list of PNG images.

    Parameters:
    - image_paths: List of file paths to the PNG images.
    - output_path: Path to save the generated GIF (including .gif extension).
    - duration: Duration for each frame in milliseconds. Default is 500 ms.
    - loop: Number of loops for the GIF. Default is 0 (infinite loop).
    - deletion: Boolean to delete images after creating the GIF. Default is False.
    - bounce: Boolean to enable bouncing effect. Default is False.

    Returns:
    - None
    """
    if not image_paths:
        raise ValueError("The list of image paths is empty.")
    
    images = [Image.open(img) for img in image_paths]
    
    images = [img.convert("RGB") if img.mode != "RGB" else img for img in images]
    
    # Handle bouncing effect
    if bounce:
        images += images[::-1]

    # Save as GIF
    images[0].save(
        output_path,
        save_all=True,
        append_images=images[1:],  # Add the rest of the images
        duration=duration,        # Duration per frame
        loop=loop                 # Number of loops
    )
    print(f"GIF saved at {output_path}")

    if deletion:
        for img_path in image_paths:
            try:
                os.remove(img_path)
                print(f"Deleted: {img_path}")
            except OSError as e:
                print(f"Error deleting {img_path}: {e}")

# Example usage
if __name__ == "__main__":
    png_files = []
    # for i in range(101):
    #     png_files.append(f"/Users/Will/Documents/FINAL_PROJECT/simulations/P_complexity_GQC_{i+1}.png")

    root = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/raw_data/np/all_crystal_grad_testruns/mu_vs_rho_mp_heatmaps/PQC_N(50-50-1)_t1(1.0-1.0-101__t2(0.0-10.0-101)_mu(0.0-10.0-101)_Delta(0.0-2.0-21)/"
    Delta_range = np.linspace(0.0, 2.0, 21)
    for Delta in Delta_range:
        Delta_safe = str(Delta).replace(".", "p")
        png_files.append(f"{root}/_PQC_sig3p0_N50_Delta{Delta_safe}_mp.png")
    
    output_gif = f"{root}_sig3p0_mp.gif"
    
    # Ensure all files exist
    png_files = [f for f in png_files if os.path.exists(f)]
    
    create_gif(png_files, output_gif, duration=500, loop=10, deletion=False, bounce=False)
