# # Optimized for speed: Adds parallel processing for frame conversion

# from PIL import Image, ImageDraw, ImageFont
# import os
# import re
# from pdf2image import convert_from_path
# from concurrent.futures import ProcessPoolExecutor, as_completed  # For parallel processing

# def create_gif(images, output_path, duration=500, loop=0):
#     """Creates a GIF from a list of PIL Image objects."""
#     if not images:
#         raise ValueError("The list of images is empty.")
    
#     images[0].save(
#         output_path,
#         save_all=True,
#         append_images=images[1:],
#         duration=duration,
#         loop=loop,
#         disposal=2
#     )
#     print(f"GIF saved at {output_path}")

# def extract_mu(filename):
#     match = re.search(r"mu_([-]?[\d\.]+)_", filename)
#     return float(match.group(1)) if match else float('inf')

# def process_single_frame(mu_filepath_tuple):
#     """Helper for parallel processing: (mu, filepath) -> processed image or None."""
#     mu, filepath = mu_filepath_tuple
#     try:
#         images = convert_from_path(filepath, dpi=96)  # Lower DPI for speed
#         if not images:
#             return None
#         img = images[0].convert("RGB")
        
#         draw = ImageDraw.Draw(img)
#         try:
#             font = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 36)
#         except IOError:
#             font = ImageFont.load_default()
        
#         text = f"mu = {mu:.3f}"
#         bbox = draw.textbbox((0, 0), text, font=font)
#         text_width, text_height = bbox[2] - bbox[0], bbox[3] - bbox[1]
#         width, height = img.size
#         x, y = width - text_width - 20, height - text_height - 20
#         draw.text((x, y), text, font=font, fill=(0, 0, 0))
        
#         return img
#     except Exception as e:
#         print(f"Error processing {filepath}: {e}")
#         return None

# if __name__ == "__main__":
#     parent_root = "/path/to/your/pdf/directory"  # Update path
    
#     subfolders = [os.path.join(parent_root, d) for d in os.listdir(parent_root) if os.path.isdir(os.path.join(parent_root, d))]
    
#     for root in subfolders:
#         print(f"\nProcessing subfolder: {root}")
        
#         files_with_mu = []
#         for filename in os.listdir(root):
#             if filename.endswith(".pdf") and "local_mp_scales" in filename:
#                 mu = extract_mu(filename)
#                 if mu != float('inf'):
#                     files_with_mu.append((mu, os.path.join(root, filename)))
        
#         files_with_mu.sort(key=lambda x: x[0])
        
#         if not files_with_mu:
#             print("No matching files. Skipping.")
#             continue
        
#         print(f"Found {len(files_with_mu)} files. Processing in parallel...")
        
#         # Parallel processing
#         processed_images = []
#         with ProcessPoolExecutor() as executor:
#             futures = [executor.submit(process_single_frame, item) for item in files_with_mu]
#             for future in as_completed(futures):
#                 img = future.result()
#                 if img:
#                     processed_images.append(img)
        
#         # Sort images back by mu (since parallel order is not guaranteed)
#         processed_images.sort(key=lambda img: float(img.filename.split('mu_')[1].split('_')[0]) if hasattr(img, 'filename') else 0)
        
#         output_gif = os.path.join(root, "local_mp_scales_ordered_by_mu.gif")
#         create_gif(processed_images, output_gif, duration=500, loop=0)
#         output_gif_fast = os.path.join(root, "local_mp_scales_ordered_by_mu_fast.gif")
#         create_gif(processed_images, output_gif_fast, duration=100, loop=0)




# ## for converting eigenvalues_vs_rho_with_mp plots into gif
# import os
# import re
# from PIL import Image, ImageDraw, ImageFont
# from pdf2image import convert_from_path

# def extract_mu_from_filename(filename):
#     """Extracts the mu value from filenames like eigenvalues_vs_rho_with_mp_N200_Delta0.05_mu0.0_seq(...).pdf"""
#     match = re.search(r"_mu([-\d\.]+)", filename)
#     return float(match.group(1)) if match else float('inf')

# def create_gif_from_pdfs(folder_path, output_gif_path, duration=500, loop=0):
#     pdf_files = [
#         f for f in os.listdir(folder_path)
#         if f.endswith('.pdf') and 'mu' in f
#     ]
#     files_with_mu = [
#         (extract_mu_from_filename(f), os.path.join(folder_path, f))
#         for f in pdf_files
#         if extract_mu_from_filename(f) != float('inf')
#     ]
#     files_with_mu.sort(key=lambda x: x[0])

#     images = []
#     for mu, filepath in files_with_mu:
#         pil_images = convert_from_path(filepath, dpi=96)
#         if not pil_images:
#             continue
#         img = pil_images[0].convert("RGB")
#         draw = ImageDraw.Draw(img)
#         try:
#             font = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 36)
#         except IOError:
#             font = ImageFont.load_default()
#         text = f"mu = {mu:.3f}"
#         bbox = draw.textbbox((0, 0), text, font=font)
#         text_width, text_height = bbox[2] - bbox[0], bbox[3] - bbox[1]
#         width, height = img.size
#         x, y = width - text_width - 20, height - text_height - 20
#         draw.text((x, y), text, font=font, fill=(0, 0, 0))
#         images.append(img)

#     if images:
#         images[0].save(
#             output_gif_path,
#             save_all=True,
#             append_images=images[1:],
#             duration=duration,
#             loop=loop,
#             disposal=2
#         )
#         print(f"GIF saved at {output_gif_path}")
#     else:
#         print("No images to save.")

# # Example usage:
# folder = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/working_notebooks/PRB_test_results/eigs_vs_rho_with_mp/N(200-200-1)_mu(0.000-3.000-101)_t11.000-1.000-1_t20.000-1.000-101_Delta0.050-0.050-1_slope(0.618-0.618-1)_phason(0.000-0.000-1)"
# output_gif = os.path.join(folder, "eigenvalues_vs_rho_with_mp_ordered_by_mu.gif")
# create_gif_from_pdfs(folder, output_gif, duration=500, loop=0)
# output_gif = os.path.join(folder, "eigenvalues_vs_rho_with_mp_ordered_by_mu_fast.gif")
# create_gif_from_pdfs(folder, output_gif, duration=100, loop=0)





# ## for converting eigenvalues_vs_mu_with_mp plots into gif
# import os
# import re
# from PIL import Image, ImageDraw, ImageFont
# from pdf2image import convert_from_path

# def extract_t2_from_filename(filename):
#     """Extracts the t2 value from filenames like eigenvalues_vs_mu_with_mp_N200_Delta0.05_t11.0_t20.0_seq(...).pdf"""
#     match = re.search(r"_t2([-\d\.]+)", filename)
#     return float(match.group(1)) if match else float('inf')

# def create_gif_ordered_by_t2(folder_path, output_gif_path, duration=500, loop=0):
#     pdf_files = [
#         f for f in os.listdir(folder_path)
#         if f.endswith('.pdf') and 't2' in f
#     ]
#     files_with_t2 = [
#         (extract_t2_from_filename(f), os.path.join(folder_path, f))
#         for f in pdf_files
#         if extract_t2_from_filename(f) != float('inf')
#     ]
#     files_with_t2.sort(key=lambda x: x[0])

#     images = []
#     for t2, filepath in files_with_t2:
#         pil_images = convert_from_path(filepath, dpi=96)
#         if not pil_images:
#             continue
#         img = pil_images[0].convert("RGB")
#         draw = ImageDraw.Draw(img)
#         try:
#             font = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 36)
#         except IOError:
#             font = ImageFont.load_default()
#         text = f"t2 = {t2:.3f}"
#         bbox = draw.textbbox((0, 0), text, font=font)
#         text_width, text_height = bbox[2] - bbox[0], bbox[3] - bbox[1]
#         width, height = img.size
#         x, y = width - text_width - 20, height - text_height - 20
#         draw.text((x, y), text, font=font, fill=(0, 0, 0))
#         images.append(img)

#     if images:
#         images[0].save(
#             output_gif_path,
#             save_all=True,
#             append_images=images[1:],
#             duration=duration,
#             loop=loop,
#             disposal=2
#         )
#         print(f"GIF saved at {output_gif_path}")
#     else:
#         print("No images to save.")

# # Example usage:
# folder = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/working_notebooks/PRB_test_results/eigs_vs_mu_with_mp/N(200-200-1)_mu(0.000-3.000-101)_t11.000-1.000-1_t20.000-1.000-101_Delta0.050-0.050-1_slope(0.618-0.618-1)_phason(0.000-0.000-1)"
# output_gif = os.path.join(folder, "eigenvalues_vs_mu_with_mp_ordered_by_t2.gif")
# create_gif_ordered_by_t2(folder, output_gif, duration=500, loop=0)
# output_gif = os.path.join(folder, "eigenvalues_vs_mu_with_mp_ordered_by_t2_fast.gif")
# create_gif_ordered_by_t2(folder, output_gif, duration=100, loop=0)





import os
import re
from PIL import Image, ImageDraw, ImageFont
from pdf2image import convert_from_path

def _decode_compact_float(value):
    """Converts strings like '0p01' or '-3p0' into floats."""
    return float(value.replace("p", "."))


def extract_delta_from_filename(filename):
    """Extracts Delta from filenames like eigenvalues_vs_mu_rho_1p5_Delta_0p01.pdf."""
    match = re.search(r"_Delta_([-\dp]+)\.pdf$", filename)
    if not match:
        return float("inf")
    return _decode_compact_float(match.group(1))


def create_gif_ordered_by_delta(folder_path, output_gif_path, duration=500, loop=0, max_delta=None):
    pdf_files = [
        filename
        for filename in os.listdir(folder_path)
        if filename.endswith(".pdf") and "_Delta_" in filename
    ]

    files_with_delta = []
    for filename in pdf_files:
        delta = extract_delta_from_filename(filename)
        if delta == float("inf"):
            continue
        if max_delta is not None and delta > max_delta:
            continue
        files_with_delta.append((delta, os.path.join(folder_path, filename)))

    files_with_delta.sort(key=lambda item: item[0], reverse=True)

    images = []
    for delta, filepath in files_with_delta:
        pil_images = convert_from_path(filepath, dpi=96)
        if not pil_images:
            continue
        img = pil_images[0].convert("RGB")
        draw = ImageDraw.Draw(img)
        try:
            font = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 36)
        except IOError:
            font = ImageFont.load_default()

        text = f"Delta = {delta:.3f}"
        bbox = draw.textbbox((0, 0), text, font=font)
        text_width = bbox[2] - bbox[0]
        text_height = bbox[3] - bbox[1]
        width, height = img.size
        x = width - text_width - 20
        y = height - text_height - 20
        draw.text((x, y), text, font=font, fill=(0, 0, 0))
        images.append(img)

    if not images:
        print("No images to save.")
        return

    images[0].save(
        output_gif_path,
        save_all=True,
        append_images=images[1:],
        duration=duration,
        loop=loop,
        disposal=2,
    )
    print(f"GIF saved at {output_gif_path}")


# Example usage:
folder = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/results/rho-Delta_vary/hp/GQC_N(500-500-1)_t1(1.0-1.0-1)_t2(1.5-1.5-1)_mu(0.0-3.0-501)_Delta(0.01-1.0-101)/rho_1p5/"
max_delta = 1.0

# output_gif = os.path.join(folder, f"eigenvalues_vs_mu_ordered_by_delta_to_{str(max_delta).replace('.', 'p')}.gif")
# create_gif_ordered_by_delta(folder, output_gif, duration=500, loop=0, max_delta=max_delta)

output_gif = os.path.join(folder, f"eigenvalues_vs_mu_ordered_by_delta_to_{str(max_delta).replace('.', 'p')}_fast.gif")
create_gif_ordered_by_delta(folder, output_gif, duration=150, loop=0, max_delta=max_delta)