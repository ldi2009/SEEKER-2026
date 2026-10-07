import os
import pandas as pd
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
file_path = os.path.join(REPO, "data", "top10000_combine_9mer_2308-Flu-mut-SCD-round123.csv")
OUT_DIR = os.path.join(REPO, "results", "flu_mutation_screen")
os.makedirs(OUT_DIR, exist_ok=True)

# Check that the input file exists
if not os.path.exists(file_path):
    raise FileNotFoundError(f"File {file_path} not found. Please check the path.")

# Read the file and take the top 300 rows with per-round frequency columns
df = pd.read_csv(file_path, sep=",")  # comma-separated
peptide_sequences = df.iloc[:300]  # top 300 rows
peptide_sequences = peptide_sequences[~peptide_sequences['aa'].str.contains(r'\*', na=False)]


# Print column names for inspection
print("Columns in peptide_sequences:", peptide_sequences.columns)

# Clean column names
peptide_sequences.columns = peptide_sequences.columns.str.strip()

# Verify the round columns are present
round_columns = [
    "count_pct_round1",
    "count_pct_round2",
    "count_pct_round3",
    #"count_pct_round4"
]

# Check for missing columns
missing_columns = [col for col in round_columns if col not in peptide_sequences.columns]
if missing_columns:
    raise KeyError(f"The following columns are missing in the file: {missing_columns}")

# Extract the per-round data
round_data = peptide_sequences[round_columns]
print("Extracted round data:")
print(round_data.head())

# Sort by round 3 frequency
sorted_indices = round_data["count_pct_round3"].sort_values(ascending=False).index
sorted_round_data = round_data.loc[sorted_indices]

# Cumulative percentage table
cumulative_data = pd.DataFrame(index=sorted_round_data.index, columns=round_columns)

# Cumulative sum per round (in round-3 sort order)
for round_name in round_columns:
    cumulative_percentages = sorted_round_data[round_name].cumsum()
    cumulative_data[round_name] = cumulative_percentages

# Plot the cumulative percentage curves
plt.figure(figsize=(10, 6))

# One line per peptide
for peptide_idx in cumulative_data.index:
    plt.plot(
        ["Round 1", "Round 2", "Round 3", ],  # x axis: sorting rounds
        #["Round 1", "Round 2", "Round 3", "Round 4"],
        cumulative_data.loc[peptide_idx],  # y axis: cumulative percentage
        # marker="o",
        alpha=0.5,
        label=f"Peptide {peptide_idx}" if peptide_idx < 5 else ""  # legend for the first 5 peptides only
    )

    # Optional: annotate peptide sequences at the last round (disabled)
    # peptide_name = peptide_sequences.iloc[peptide_idx, 1]
    # plt.text(
    #     x=3,
    #     y=cumulative_data.loc[peptide_idx, "count_pct_round4"],
    #     s=peptide_name,
    #     fontsize=8,
    #     alpha=0.7
    # )

# Legend and axis labels
plt.xlabel("Rounds")
plt.ylabel("Cumulative Percentage")
plt.title("Cumulative Percentage of Peptides Across Rounds (Sorted by Round 3)")
plt.tight_layout()

# Save the figure
plt.savefig(os.path.join(OUT_DIR, "cumulative_percentage_peptides_sorted_by_round3_no_labels.pdf"), format="pdf")
plt.close()