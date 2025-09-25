import json
import pandas as pd
import matplotlib.pyplot as plt


def json_to_dataframe(filename):
    """
    Converts the JSON data into a Pandas DataFrame, handling both formats of eigenvalues.

    Parameters:
    - filename (str): Path to the JSON file.

    Returns:
    - df (pd.DataFrame): DataFrame containing extracted data.
    """
    with open(filename, "r") as f:
        data = json.load(f)

    structured_data = []

    for entry in data:
        # Check if eigenvalues are stored as dictionaries (complex format) or as raw numbers (real format)
        if isinstance(entry["eigenvalues"][0], dict):  # Complex format
            eigenvalues = [eig["re"] for eig in entry["eigenvalues"]]  # Extract real part
        else:  # Real format
            eigenvalues = entry["eigenvalues"]

        structured_data.append({
            "N": entry["N"],
            "t_n": tuple(entry["t_n"]),  # Convert list to tuple for easier indexing
            "mu": entry["mu"],
            "Delta": entry["Delta"],
            "sequence_name": entry["sequence_name"],
            "eigenvalues": eigenvalues,
        })

    # Convert to DataFrame
    df = pd.DataFrame(structured_data)

    return df

def plot_eigenvalues(df, fixed_params, save_filename):
    """
    Plots eigenvalues vs. mu for a fixed set of parameters.

    Parameters:
    - df (pd.DataFrame): DataFrame containing eigenvalues and parameters.
    - fixed_params (dict): Dictionary specifying fixed values for parameters except mu.
    """

    mask = (df[list(fixed_params)] == pd.Series(fixed_params)).all(axis=1)
    filtered_df = df[mask]

    if filtered_df.empty:
        print("No matching data found for the given fixed parameters.")
        return

    filtered_df = filtered_df.sort_values("mu")

    mu_values = filtered_df["mu"].tolist()
    eigenvalues = filtered_df["eigenvalues"].tolist()

    # Convert to transposed format (each eigenvalue set becomes a separate line)
    eigenvalues = list(map(list, zip(*eigenvalues)))

    plt.figure(figsize=(8, 6))
    for eig_set in eigenvalues:
        plt.plot(mu_values, eig_set, linestyle="-", alpha=0.7)

    plt.xlabel(r"$\mu$")
    plt.ylabel("Eigenvalues")
    plt.title(rf"Eigenvalues vs. $\mu$ for sequence '{fixed_params['sequence_name']}'")
    
    # plt.xlim(0.0, 0.5)
    # plt.ylim(-1e-5, 1e-5)

    plt.grid(True)
    # plt.show()
    plt.savefig(f"simulations/system_eigensolver/{save_filename}.png", dpi=800)



# Usage
filename = "simulations/system_eigensolver/results.json"

# # Show JSON structure if needed for debugging
# with open(filename, 'r') as f:
#     data = json.load(f)
# print(data[0])

df = json_to_dataframe(filename)

fixed_params = {
    "N": 100,
    "t_n": (1.0,),
    "Delta": 1.0,
    "sequence_name": "NC"
}

save_figname = "arpack_eigenplot"
plot_eigenvalues(df, fixed_params, save_figname)
