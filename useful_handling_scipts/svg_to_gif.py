from PIL import Image, ImageDraw, ImageFont
import os
import numpy as np
import re
import io
import tempfile
import subprocess


# try:
#     import cairosvg
# except ImportError:
#     print("Error: 'cairosvg' library is missing.")
#     print("Please install it running: pip install cairosvg")
#     # On macOS you may also need: brew install cairo
#     cairosvg = None

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

def extract_phi(filename):
    # Extracts the phi value using regex
    match = re.search(r"phi([\d\.]+)_", filename)
    if match:
        return float(match.group(1))
    return -1.0


def process_frame(filepath, phi_val):
    """
    Converts SVG to PNG using Inkscape and draws the phi value.
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
    
    # Text to draw
    text = f"phi = {phi_val:.3f}"
    
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

    root = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/results/bp_results/MBData_final/hof_style_slopes_N1000_phason_0.0-1-0.0_nbins1000_npb1_N(500-500-1)_t1(1p0-1p0-1)_t2(2p0-2p0-1)_mu(0p0-5p0-146627)_Delta(0p05-0p05-1)/mp_gap_comp/mptarg-1.0_mptol0.01_idos0.05,0.004,200_N_rangenothing_Delta_rangenothing_tn_rangenothing_phi_rangenothing_phason_range[0.0]/bulk_vs_qc_gaps"
    
    # 1. Collect all valid files and their phi values
    files_with_phi = []
    
    print("Scanning directory...")
    for filename in os.listdir(root):
        if filename.endswith(".svg") and "bulk_vs_qc_gaps_robust_from_eigs_robust_direct_logscalefalse" in filename:
            phi = extract_phi(filename)
            if phi != -1.0:
                files_with_phi.append((phi, os.path.join(root, filename)))
    
    # 2. Sort by phi
    files_with_phi.sort(key=lambda x: x[0])
    
    if not files_with_phi:
        print("No matching files found.")
        exit()

    print(f"Found {len(files_with_phi)} files. Processing frames...")

    # 3. Process images
    processed_images = []
    for i, (phi, filepath) in enumerate(files_with_phi):
        print(f"[{i+1}/{len(files_with_phi)}] Processing phi={phi:.4f}...")
        img = process_frame(filepath, phi)
        if img:
            processed_images.append(img)
    
    # 4. Save GIF
    output_gif = os.path.join(root, "bulk_vs_qc_gaps_ordered_by_phi.gif")
    create_gif(processed_images, output_gif, duration=500, loop=0)