# from PIL import Image, ImageDraw, ImageFont
# import os
# import numpy as np
# import re
# import io
# import tempfile
# import subprocess


# # try:
# #     import cairosvg
# # except ImportError:
# #     print("Error: 'cairosvg' library is missing.")
# #     print("Please install it running: pip install cairosvg")
# #     # On macOS you may also need: brew install cairo
# #     cairosvg = None

# def create_gif(images, output_path, duration=500, loop=0):
#     """
#     Creates a GIF from a list of PIL Image objects.
#     """
#     if not images:
#         raise ValueError("The list of images is empty.")
    
#     # Save as GIF
#     images[0].save(
#         output_path,
#         save_all=True,
#         append_images=images[1:],  # Add the rest of the images
#         duration=duration,         # Duration per frame
#         loop=loop,                 # Number of loops
#         disposal=2                 # Restore to background color (helps with transparency)
#     )
#     print(f"GIF saved at {output_path}")

# def extract_phi(filename):
#     # Extracts the phi value using regex
#     match = re.search(r"phi([\d\.]+)_", filename)
#     if match:
#         return float(match.group(1))
#     return -1.0


# def process_frame(filepath, phi_val):
#     """
#     Converts SVG to PNG using Inkscape and draws the phi value.
#     """
#     # Create a temporary PNG file
#     with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp_file:
#         tmp_png_path = tmp_file.name
    
#     # Use Inkscape to convert SVG to PNG (full path to CLI)
#     try:
#         subprocess.run([
#             '/Applications/Inkscape.app/Contents/MacOS/inkscape',  # Full path to CLI
#             filepath,
#             '--export-filename', tmp_png_path,
#             '--export-type=png',
#             '--export-dpi', '96'  # Adjust DPI for resolution (96 is standard)
#         ], check=True, capture_output=True)
#     except subprocess.CalledProcessError as e:
#         print(f"Error converting {filepath}: {e}")
#         os.unlink(tmp_png_path)
#         return None
    
#     # Load the PNG with PIL
#     try:
#         img = Image.open(tmp_png_path).convert("RGB")
#     except Exception as e:
#         print(f"Error loading PNG: {e}")
#         os.unlink(tmp_png_path)
#         return None
    
#     # Clean up temp file
#     os.unlink(tmp_png_path)
    
#     # Setup drawing
#     draw = ImageDraw.Draw(img)
    
#     # Try to load a nice font, otherwise use default
#     try:
#         font = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 36)
#     except IOError:
#         font = ImageFont.load_default()
    
#     # Text to draw
#     text = f"phi = {phi_val:.3f}"
    
#     # Calculate text position (bottom right)
#     bbox = draw.textbbox((0, 0), text, font=font)
#     text_width = bbox[2] - bbox[0]
#     text_height = bbox[3] - bbox[1]
    
#     width, height = img.size
#     margin = 20
#     x = width - text_width - margin
#     y = height - text_height - margin
    
#     # Draw text with black color
#     draw.text((x, y), text, font=font, fill=(0, 0, 0))
    
#     return img

# # Example usage
# if __name__ == "__main__":

#     root = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/results/bp_results/MBData_final/hof_style_slopes_N1000_phason_0.0-1-0.0_nbins1000_npb1_N(500-500-1)_t1(1p0-1p0-1)_t2(2p0-2p0-1)_mu(0p0-5p0-146627)_Delta(0p05-0p05-1)/mp_gap_comp/mptarg-1.0_mptol0.01_idos0.05,0.004,200_N_rangenothing_Delta_rangenothing_tn_rangenothing_phi_rangenothing_phason_range[0.0]/bulk_vs_qc_gaps"
    
#     # 1. Collect all valid files and their phi values
#     files_with_phi = []
    
#     print("Scanning directory...")
#     for filename in os.listdir(root):
#         if filename.endswith(".svg") and "bulk_vs_qc_gaps_robust_from_eigs_robust_direct_logscalefalse" in filename:
#             phi = extract_phi(filename)
#             if phi != -1.0:
#                 files_with_phi.append((phi, os.path.join(root, filename)))
    
#     # 2. Sort by phi
#     files_with_phi.sort(key=lambda x: x[0])
    
#     if not files_with_phi:
#         print("No matching files found.")
#         exit()

#     print(f"Found {len(files_with_phi)} files. Processing frames...")

#     # 3. Process images
#     processed_images = []
#     for i, (phi, filepath) in enumerate(files_with_phi):
#         print(f"[{i+1}/{len(files_with_phi)}] Processing phi={phi:.4f}...")
#         img = process_frame(filepath, phi)
#         if img:
#             processed_images.append(img)
    
#     # 4. Save GIF
#     output_gif = os.path.join(root, "bulk_vs_qc_gaps_ordered_by_phi.gif")
#     create_gif(processed_images, output_gif, duration=500, loop=0)



# # ------------------------------------------------------------------------------------------------------------------
# # ------------------------------------------------------------------------------------------------------------------
## # For all eigs vs mu in dir, ordered by slope extracted from filename
# from PIL import Image, ImageDraw, ImageFont
# import os
# import numpy as np
# import re
# import io
# import tempfile
# import subprocess

# def create_gif(images, output_path, duration=500, loop=0):
#     """
#     Creates a GIF from a list of PIL Image objects.
#     """
#     if not images:
#         raise ValueError("The list of images is empty.")
    
#     # Save as GIF
#     images[0].save(
#         output_path,
#         save_all=True,
#         append_images=images[1:],  # Add the rest of the images
#         duration=duration,         # Duration per frame
#         loop=loop,                 # Number of loops
#         disposal=2                 # Restore to background color (helps with transparency)
#     )
#     print(f"GIF saved at {output_path}")

# def extract_slope(filename):
#     # Extracts the first value (slope) from seq(first, second) using regex
#     match = re.search(r"seq\(([\d\.]+),\s*([\d\.]+)\)", filename)
#     if match:
#         return float(match.group(1))  # First value in the tuple
#     return -1.0

# def process_frame(filepath, slope_val):
#     """
#     Converts SVG to PNG using Inkscape and draws the slope value.
#     """
#     # Create a temporary PNG file
#     with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp_file:
#         tmp_png_path = tmp_file.name
    
#     # Use Inkscape to convert SVG to PNG (full path to CLI)
#     try:
#         subprocess.run([
#             '/Applications/Inkscape.app/Contents/MacOS/inkscape',  # Full path to CLI
#             filepath,
#             '--export-filename', tmp_png_path,
#             '--export-type=png',
#             '--export-dpi', '96'  # Adjust DPI for resolution (96 is standard)
#         ], check=True, capture_output=True)
#     except subprocess.CalledProcessError as e:
#         print(f"Error converting {filepath}: {e}")
#         os.unlink(tmp_png_path)
#         return None
    
#     # Load the PNG with PIL
#     try:
#         img = Image.open(tmp_png_path).convert("RGB")
#     except Exception as e:
#         print(f"Error loading PNG: {e}")
#         os.unlink(tmp_png_path)
#         return None
    
#     # Clean up temp file
#     os.unlink(tmp_png_path)
    
#     # Setup drawing
#     draw = ImageDraw.Draw(img)
    
#     # Try to load a nice font, otherwise use default
#     try:
#         font = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 36)
#     except IOError:
#         font = ImageFont.load_default()
    
#     # Text to draw (show slope value)
#     text = f"slope = {slope_val:.3f}"
    
#     # Calculate text position (bottom right)
#     bbox = draw.textbbox((0, 0), text, font=font)
#     text_width = bbox[2] - bbox[0]
#     text_height = bbox[3] - bbox[1]
    
#     width, height = img.size
#     margin = 20
#     x = width - text_width - margin
#     y = height - text_height - margin
    
#     # Draw text with black color
#     draw.text((x, y), text, font=font, fill=(0, 0, 0))
    
#     return img

# # Example usage
# if __name__ == "__main__":
#     root = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/working_notebooks/PRB_test_results/eigs_vs_mu/N250/svg"
    
#     # 1. Collect all valid files and their slope values
#     files_with_slope = []
    
#     print("Scanning directory...")
#     for filename in os.listdir(root):
#         if filename.endswith(".svg") and "eigenvalues_vs_mu" in filename:  # Filter for your file type
#             slope = extract_slope(filename)
#             if slope != -1.0:
#                 files_with_slope.append((slope, os.path.join(root, filename)))
    
#     # 2. Sort by slope (smallest to largest)
#     files_with_slope.sort(key=lambda x: x[0])
    
#     if not files_with_slope:
#         print("No matching files found.")
#         exit()

#     print(f"Found {len(files_with_slope)} files. Processing frames...")

#     # 3. Process images
#     processed_images = []
#     for i, (slope, filepath) in enumerate(files_with_slope):
#         print(f"[{i+1}/{len(files_with_slope)}] Processing slope={slope:.4f}...")
#         img = process_frame(filepath, slope)
#         if img:
#             processed_images.append(img)
    
#     # 4. Save GIF
#     output_gif = os.path.join(root, "eigenvalues_vs_mu_ordered_by_slope.gif")
#     create_gif(processed_images, output_gif, duration=500, loop=0)





# # ------------------------------------------------------------------------------------------------------------------
# # ------------------------------------------------------------------------------------------------------------------

# ## # For Local wavefunctions (all in dir, or target subdir)
# from PIL import Image, ImageDraw, ImageFont
# import os
# import numpy as np
# import re
# import io
# import tempfile
# import subprocess

# def create_gif(images, output_path, duration=500, loop=0):
#     """
#     Creates a GIF from a list of PIL Image objects.
#     """
#     if not images:
#         raise ValueError("The list of images is empty.")
    
#     # Save as GIF
#     images[0].save(
#         output_path,
#         save_all=True,
#         append_images=images[1:],  # Add the rest of the images
#         duration=duration,         # Duration per frame
#         loop=loop,                 # Restore to background color (helps with transparency)
#         disposal=2
#     )
#     print(f"GIF saved at {output_path}")

# def extract_mu(filename):
#     # Extracts the mu value using regex (handles negatives like mu-0.3)
#     match = re.search(r"mu([-]?[\d\.]+)_", filename)
#     if match:
#         return float(match.group(1))
#     return float('inf')  # Use inf if not found, so it sorts last

# def process_frame(filepath, mu_val):
#     """
#     Converts SVG to PNG using Inkscape and draws the mu value.
#     """
#     # Create a temporary PNG file
#     with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp_file:
#         tmp_png_path = tmp_file.name
    
#     # Use Inkscape to convert SVG to PNG (full path to CLI)
#     try:
#         subprocess.run([
#             '/Applications/Inkscape.app/Contents/MacOS/inkscape',  # Full path to CLI
#             filepath,
#             '--export-filename', tmp_png_path,
#             '--export-type=png',
#             '--export-dpi', '96'  # Adjust DPI for resolution (96 is standard)
#         ], check=True, capture_output=True)
#     except subprocess.CalledProcessError as e:
#         print(f"Error converting {filepath}: {e}")
#         os.unlink(tmp_png_path)
#         return None
    
#     # Load the PNG with PIL
#     try:
#         img = Image.open(tmp_png_path).convert("RGB")
#     except Exception as e:
#         print(f"Error loading PNG: {e}")
#         os.unlink(tmp_png_path)
#         return None
    
#     # Clean up temp file
#     os.unlink(tmp_png_path)
    
#     # Setup drawing
#     draw = ImageDraw.Draw(img)
    
#     # Try to load a nice font, otherwise use default
#     try:
#         font = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 36)
#     except IOError:
#         font = ImageFont.load_default()
    
#     # Text to draw (show mu value)
#     text = f"mu = {mu_val:.3f}"
    
#     # Calculate text position (bottom right)
#     bbox = draw.textbbox((0, 0), text, font=font)
#     text_width = bbox[2] - bbox[0]
#     text_height = bbox[3] - bbox[1]
    
#     width, height = img.size
#     margin = 20
#     x = width - text_width - margin
#     y = height - text_height - margin
    
#     # Draw text with black color
#     draw.text((x, y), text, font=font, fill=(0, 0, 0))
    
#     return img

# # Example usage
# if __name__ == "__main__":
#     parent_root = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/working_notebooks/PRB_test_results/local_wavefunctions/N100_Delta1.0_rho2.0/svg"
    
#     # Get all subdirectories in the parent root
#     subfolders = [os.path.join(parent_root, d) for d in os.listdir(parent_root) if os.path.isdir(os.path.join(parent_root, d))]

#     # Only specific target subfolder (uncomment if needed)
#     subfolder = os.path.join(parent_root, "seq_(0.0, 0.0)")
#     if os.path.isdir(subfolder):
#         subfolders = [subfolder]
#     else:
#         print(f"Warning: Target subfolder {subfolder} not found. Processing all subfolders.")

    
#     for root in subfolders:
#         print(f"\nProcessing subfolder: {root}")
        
#         # 1. Collect all valid files and their mu values
#         files_with_mu = []
        
#         print("Scanning subfolder...")
#         for filename in os.listdir(root):
#             if filename.endswith(".svg") and "local_wavefunction" in filename:
#                 mu = extract_mu(filename)
#                 if mu != float('inf'):
#                     files_with_mu.append((mu, os.path.join(root, filename)))
        
#         # 2. Sort by mu (smallest to largest, handles negatives)
#         files_with_mu.sort(key=lambda x: x[0])
        
#         if not files_with_mu:
#             print("No matching files found in this subfolder. Skipping.")
#             continue
        
#         print(f"Found {len(files_with_mu)} files. Processing frames...")
        
#         # 3. Process images
#         processed_images = []
#         for i, (mu, filepath) in enumerate(files_with_mu):
#             print(f"[{i+1}/{len(files_with_mu)}] Processing mu={mu:.4f}...")
#             img = process_frame(filepath, mu)
#             if img:
#                 processed_images.append(img)
        
#         # 4. Save GIF in the subfolder
#         output_gif = os.path.join(root, "local_wavefunction_ordered_by_mu.gif")
#         create_gif(processed_images, output_gif, duration=500, loop=0)
#         output_gif_fast = os.path.join(root, "local_wavefunction_ordered_by_mu_fast.gif")
#         create_gif(processed_images, output_gif_fast, duration=100, loop=0)






# # ------------------------------------------------------------------------------------------------------------------
# # ------------------------------------------------------------------------------------------------------------------
# ## # For local mp scales (all in dir, or target subdir)

from PIL import Image, ImageDraw, ImageFont
import os
import numpy as np
import re
import io
import tempfile
import subprocess

def create_gif(images, output_path, duration=500, loop=0):
    """
    Creates a GIF from a list of PIL Image objects.
    """
    if not images:
        raise ValueError("The list of images is empty.")
    
    # Save as GIF
    images[0].save(
        output_path,
        save_all=True,
        append_images=images[1:],  # Add the rest of the images
        duration=duration,         # Duration per frame
        loop=loop,                 # Number of loops
        disposal=2                 # Restore to background color (helps with transparency)
    )
    print(f"GIF saved at {output_path}")

def extract_mu(filename):
    # Extracts the mu value using regex (handles negatives like mu_-0.3)
    match = re.search(r"mu_([-]?[\d\.]+)_", filename)
    if match:
        return float(match.group(1))
    return float('inf')  # Use inf if not found, so it sorts last

def process_frame(filepath, mu_val):
    """
    Converts SVG to PNG using Inkscape and draws the mu value.
    """
    # Create a temporary PNG file
    with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp_file:
        tmp_png_path = tmp_file.name
    
    # Use Inkscape to convert SVG to PNG (full path to CLI)
    try:
        subprocess.run([
            '/Applications/Inkscape.app/Contents/MacOS/inkscape',  # Full path to CLI
            filepath,
            '--export-filename', tmp_png_path,
            '--export-type=png',
            '--export-dpi', '96'  # Adjust DPI for resolution (96 is standard)
        ], check=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        print(f"Error converting {filepath}: {e}")
        os.unlink(tmp_png_path)
        return None
    
    # Load the PNG with PIL
    try:
        img = Image.open(tmp_png_path).convert("RGB")
    except Exception as e:
        print(f"Error loading PNG: {e}")
        os.unlink(tmp_png_path)
        return None
    
    # Clean up temp file
    os.unlink(tmp_png_path)
    
    # Setup drawing
    draw = ImageDraw.Draw(img)
    
    # Try to load a nice font, otherwise use default
    try:
        font = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 36)
    except IOError:
        font = ImageFont.load_default()
    
    # Text to draw (show mu value)
    text = f"mu = {mu_val:.3f}"
    
    # Calculate text position (bottom right)
    bbox = draw.textbbox((0, 0), text, font=font)
    text_width = bbox[2] - bbox[0]
    text_height = bbox[3] - bbox[1]
    
    width, height = img.size
    margin = 20
    x = width - text_width - margin
    y = height - text_height - margin
    
    # Draw text with black color
    draw.text((x, y), text, font=font, fill=(0, 0, 0))
    
    return img

# Example usage
if __name__ == "__main__":
    parent_root = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/working_notebooks/PRB_test_results/local_mp_scales/check_2_N(100-100-1)_mu(-3.000-3.000-101)_t11.000-1.000-1_t21.000-1.000-1_Delta0.500-0.500-1_slope(0.618-0.618-1)_phason(0.000-0.000-1)/svg"
    

    # Get all subdirectories in the parent root
    subfolders = [os.path.join(parent_root, d) for d in os.listdir(parent_root) if os.path.isdir(os.path.join(parent_root, d))]

    # # Target specific subfolder
    # subfolder = os.path.join(parent_root, "seq_(0.0, 0.0)")
    # if os.path.isdir(subfolder):
    #     subfolders = [subfolder]
    # else:
    #     print(f"Error: Target subfolder {subfolder} not found.")
    #     exit()
    
    for root in subfolders:
        print(f"\nProcessing subfolder: {root}")
        
        # 1. Collect all valid files and their mu values
        files_with_mu = []
        
        print("Scanning subfolder...")
        for filename in os.listdir(root):
            if filename.endswith(".svg") and "local_mp_scales" in filename:
                mu = extract_mu(filename)
                if mu != float('inf'):
                    files_with_mu.append((mu, os.path.join(root, filename)))
        
        # 2. Sort by mu (smallest to largest, handles negatives)
        files_with_mu.sort(key=lambda x: x[0])
        
        if not files_with_mu:
            print("No matching files found in this subfolder. Skipping.")
            continue
        
        print(f"Found {len(files_with_mu)} files. Processing frames...")
        
        # 3. Process images
        processed_images = []
        for i, (mu, filepath) in enumerate(files_with_mu):
            print(f"[{i+1}/{len(files_with_mu)}] Processing mu={mu:.4f}...")
            img = process_frame(filepath, mu)
            if img:
                processed_images.append(img)
        
        # 4. Save GIF in the subfolder
        output_gif = os.path.join(root, "local_mp_scales_ordered_by_mu.gif")
        create_gif(processed_images, output_gif, duration=500, loop=0)
        output_gif_fast = os.path.join(root, "local_mp_scales_ordered_by_mu_fast.gif")
        create_gif(processed_images, output_gif_fast, duration=100, loop=0)