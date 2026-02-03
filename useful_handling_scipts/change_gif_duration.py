import sys
from PIL import Image
import os

def change_gif_duration(input_path, output_path, new_duration_ms):
    """
    Change the frame duration of an existing GIF.
    
    Parameters:
    - input_path: Path to the input GIF file.
    - output_path: Path to save the modified GIF.
    - new_duration_ms: New duration for each frame in milliseconds.
    """
    with Image.open(input_path) as im:
        frames = []
        try:
            while True:
                frames.append(im.copy())
                im.seek(im.tell() + 1)
        except EOFError:
            pass
        
        if not frames:
            raise ValueError("No frames found in the GIF.")
        
        # Save with new duration
        frames[0].save(
            output_path,
            save_all=True,
            append_images=frames[1:],
            duration=new_duration_ms,
            loop=0  # Infinite loop; change if needed
        )
        print(f"GIF saved to {output_path} with frame duration {new_duration_ms} ms.")

if __name__ == "__main__":
    # Option 1: Hardcode the paths and duration here (uncomment and set these)
    input_path = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/results/bp_results/MBData_final/hof_style_slopes_N1000_phason_0.0-1-0.0_nbins1000_npb1_N(500-500-1)_t1(1p0-1p0-1)_t2(2p0-2p0-1)_mu(0p0-5p0-146627)_Delta(0p05-0p05-1)/mp_gap_comp/mptarg-1.0_mptol0.01_idos0.05,0.004,200_N_rangenothing_Delta_rangenothing_tn_rangenothing_phi_rangenothing_phason_range[0.0]/bulk_vs_qc_gaps/bulk_vs_qc_gaps_ordered_by_phi.gif"
    output_path = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/results/bp_results/MBData_final/hof_style_slopes_N1000_phason_0.0-1-0.0_nbins1000_npb1_N(500-500-1)_t1(1p0-1p0-1)_t2(2p0-2p0-1)_mu(0p0-5p0-146627)_Delta(0p05-0p05-1)/mp_gap_comp/mptarg-1.0_mptol0.01_idos0.05,0.004,200_N_rangenothing_Delta_rangenothing_tn_rangenothing_phi_rangenothing_phason_range[0.0]/bulk_vs_qc_gaps/bulk_vs_qc_gaps_ordered_by_phi_fast.gif"
    new_duration = 100  # New duration in milliseconds per frame
    
    # os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # Run with hardcoded values
    change_gif_duration(input_path, output_path, new_duration)
    
    # Option 2: Use command-line arguments (comment out the above and uncomment below if preferred)
    # if len(sys.argv) != 4:
    #     print("Usage: python change_gif_duration.py <input.gif> <output.gif> <duration_ms>")
    #     print("Example: python change_gif_duration.py input.gif output.gif 500")
    #     sys.exit(1)
    # 
    # input_path = sys.argv[1]
    # output_path = sys.argv[2]
    # try:
    #     new_duration = int(sys.argv[3])
    # except ValueError:
    #     print("Error: duration_ms must be an integer.")
    #     sys.exit(1)
    # 
    # change_gif_duration(input_path, output_path, new_duration)