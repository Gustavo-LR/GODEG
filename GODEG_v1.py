import pandas as pd
import numpy as np
import argparse
import os
import sys
import re
import subprocess
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import shutil
from concurrent.futures import ProcessPoolExecutor
from functools import partial


# Argument parser setup
parser = argparse.ArgumentParser(
    description="Apply log2FC and adjusted P-value thresholds for a DEG table, classify the DEGs in the complete GO hierarchy and do an enrichment analysis to search overreprensented DEGs."
)
parser.add_argument(
    "-d", "--deg", required=True, help="Path to the DEG output file, provided by, i.e., DESeq2 or edgeR (TSV format), including headers"
)
parser.add_argument(
    "-pvalue", required=False, type=float, default=0.05, help="Adjusted p-value cutoff for differentially expressed gene detection (default: 0.05)"
)
parser.add_argument(
    "-fup", "--log2foldUp", required=False, type=float, default=1.0, help="Log2foldChange or FoldChange cutoff for upregulated gene detection (default: 1)"
)
parser.add_argument(
    "-fdo", "--log2foldDown", required=False, type=float, default=-1.0, help="Log2foldChange or FoldChange cutoff for downregulated gene detection (default: -1)"
)
parser.add_argument(
    "-go", required=True, help="GO annotation file for the specified go_mode, the first column containing the seqIDs and the second the corresponding GO terms (TSV), and must include headers"
)
parser.add_argument(
    "-go_mode", required=True, type=str.upper, choices=["P", "F", "C"], help="GO mode, the options are P (biological process, default), F (molecular function) and C (cellular component)"
)
parser.add_argument(
    "-a", "--annotation", required=False, help="Annotation containing the protein function for the seqIDs (TSV format), including headers for each column"
)
parser.add_argument(
    "-s", "--samples", required=True, help="Samples list from the DEG output file, in the first column the log2FC, in the second colunm the correspondent adjusted p-value and in the third the treatment name. Columns needs to be tab separated and contain no header"
)

parser.add_argument(
    "-bgs","--bgsize", required=False, type=int, help="Gene background size. By default, GODEG uses a minimum background (number of GO annotated genes + DEGs), which is the most conservative setting to identify enriched GO terms. Using the real gene background size will increase the number of enriched terms but sometimes increase false positive discovery. We advise to increase the background size only if there is a significant number of DEGs but a negligible number of enriched terms."
)

parser.add_argument(
    "-epvalue", required=False, type=float, default=0.05, help="Adjusted p-value cutoff for enriched GO terms detection (default: 0.05)."
)

parser.add_argument(
    "-t","--threads", required=True, type=int, help="Number of threads assigned for degs_by_go creation and fisher enrichment test"
)
args = parser.parse_args()


# Load DEG table file
try:
    df = pd.read_csv(args.deg, sep='\t')
except Exception as e:
    print(f"Error reading DEG table file: {e}")
    sys.exit(1)

try:
    df_go = pd.read_csv(args.go, sep='\t')
except Exception as e:
    print(f"Error reading GO annotation file: {e}")
    sys.exit(1)


if args.samples is not None:
    treatments = []
    try:
        # Read samples file (3 columns: log2fc_col, adjustedP_col, treatment)
        samples = pd.read_csv(args.samples, sep='\t', header=None, names=["log2fc", "adjustedP", "treatment"])

        # Build list of columns to keep (first column of df + all log2fc/adjustedP columns)
        cols_to_keep = [df.columns[0]] + [col for pair in samples[["log2fc","adjustedP"]].values for col in pair]

        # Check that all requested columns exist
        missing_cols = [col for col in cols_to_keep[1:] if col not in df.columns]
        if missing_cols:
            print(f"Error: The following columns from {args.samples} were not found in {args.deg}: {', '.join(missing_cols)}")
            sys.exit(1)

        # Rename columns using the treatment names provided
        rename_map = {}
        for log2fc_col, adjp_col, treatment in samples.values:
            treatments.append(treatment)
            rename_map[log2fc_col] = f"log2FC_{treatment}"
            rename_map[adjp_col]   = f"adjustedP_{treatment}"

        # Filter df and rename columns
        df = df[cols_to_keep].rename(columns=rename_map)

        # Space-separated treatments string
        treatments = " ".join(treatments)

    except Exception as e:
        print(f"Error processing samples file: {e}")
        sys.exit(1)


# Define full name for GO mode
go_mode_fullname = {
    'P': 'Biological Process',
    'C': 'Cellular Component',
    'F': 'Molecular Function'
}.get(args.go_mode)


def bgsize_check():
    if args.bgsize is not None:
        return args.bgsize
    else:
        return "Minimun (DEGS+GO Annotated genes)"

# Print parameters
print("==================== PARAMETERS FOR THIS RUN ====================")
print(f"DEG table input file(-d):                                   {args.deg}")
print(f"GO annotation file (-go):                                   {args.go}")
print(f"GO mode (-go_mode):                                         {args.go_mode} ({go_mode_fullname})")
print(f"P-value cutoff for DEG detection (-pvalue):                 {args.pvalue}")
print(f"FoldChange UpRegulated cutoff (-fup):                       > {args.log2foldUp}")
print(f"FoldChange DownRegulated cutoff (-fdo):                     < {args.log2foldDown}")
print(f"Number of threads for degs_by_go and enrichment (-t):       {args.threads}")
print(f"Enrichment background size (-bgs):                          {bgsize_check()}")
print(f"P-value cutoff for enriched GO terms detection (-epvalue):  {args.epvalue}")
print(f"Treatment conditions:                                       {treatments}")
print("=================================================================")


# Ask user for confirmation
proceed = input("Proceed with the run? (yes/no): ").strip().lower()
if proceed.lower() not in ['yes', 'y', 'ye']:
    print("Aborted by user.")
    exit()


# Create a working directory with incremental numbering
base_name = "GODEG"
counter = 1
while True:
    output_dir = f"{base_name}_{counter}"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        break
    counter += 1

print('The intermediary files for this run will be saved in the folder ' + output_dir)

if args.annotation is not None:
    # Load annotation file
    df_anno = pd.read_csv(args.annotation, sep="\t", dtype=str)

    # Key columns
    key_df = df.columns[0]
    key_anno = df_anno.columns[0]

    # Align annotation key name if needed
    if key_df != key_anno:
        df_anno = df_anno.rename(columns={key_anno: key_df})

    # Merge df with annotation
    df_merged = pd.merge(df, df_anno, on=key_df, how="left")

    # Reorder so that annotation columns appear immediately after seqID
    cols = [key_df] + list(df_anno.columns[1:]) + list(df.columns[1:])
    df_merged = df_merged[cols]

    # Fill missing annotation with empty string (optional)
    df_merged[df_anno.columns[1:]] = df_merged[df_anno.columns[1:]].fillna("no_hit")
    df = df_merged


deg_input = args.deg
go_input = args.go

# Copy input files into the new folder
shutil.copy(deg_input, os.path.join(output_dir, os.path.basename(deg_input)))
shutil.copy(go_input, os.path.join(output_dir, os.path.basename(go_input)))


# Save current working directory
original_dir = os.getcwd()
# Change to the new directory
os.chdir(output_dir)


list_of_genes = df.iloc[:, 0].tolist()

# Double columns containing 'log2FC' and 'adjustedP' and add prefix 'value_'
log2_cols = df.columns[df.columns.str.contains('log2FC')]
adjusted_cols = df.columns[df.columns.str.contains('adjustedP')]

for col in log2_cols:
    new_col = f"value_{col}"
    df.insert(df.columns.get_loc(col) + 1, new_col, df[col])

for col in adjusted_cols:
    new_col = f"value_{col}"
    df.insert(df.columns.get_loc(col) + 1, new_col, df[col])

def conditions1(x):
    if x < args.pvalue:
        return 'Significative'
    else:
        return 'null'

func1 = np.vectorize(conditions1)
copy_cols = df.columns[df.columns.str.contains('value_adjusted')]
#DEG_class1 = df[copy_cols].applymap(conditions1)
DEG_class1 = df[copy_cols].apply(lambda col: col.map(conditions1))
df[copy_cols] = DEG_class1

def conditions(x):
    x = float(x)
    if x > args.log2foldUp:
        return "UpRegulated"
    elif x < args.log2foldDown:
        return "DownRegulated"
    else:
        return "NotDEG"

func = np.vectorize(conditions)
copy_cols1 = df.columns[df.columns.str.contains('value_log2')]
#DEG_class2 = df[copy_cols1].applymap(conditions)
DEG_class2 = df[copy_cols1].apply(lambda col: col.map(conditions))
df[copy_cols1] = DEG_class2

null_idx = np.where(df == 'null')

for i, j in zip(null_idx[0], null_idx[1]):
    if j > 0:
        df.iloc[i, j - 2] = 'null'


# Write to output file
output_file_idsupdown = deg_input.replace('.tsv', '_DEG.tsv')
df.to_csv(output_file_idsupdown, sep='\t', index=False, na_rep='null')

print(f"DEGs conditions were applied in {output_file_idsupdown} (Result output 1)")


data = pd.read_csv(output_file_idsupdown, sep='\t')

data = data[data.apply(lambda row: row.astype(str).str.contains('UpRegulated|DownRegulated').any(), axis=1)]

# Function to identify numeric-only columns
def is_numeric_only(col):
    return pd.to_numeric(col, errors='coerce').notna().all()

# Keep the first column always, also keep any column with 'log2FC' or 'adjusted', or any that is not numeric-only
cols_to_keep = [
    col for i, col in enumerate(data.columns)
    if i == 0 or 'log2FC' in col or 'adjusted' in col or not is_numeric_only(data[col])
]

filtered_data = data[cols_to_keep]

filtered_data = filtered_data.fillna('null')

for col in filtered_data.columns:
    if 'value_log2FC' in col:
        col_idx = filtered_data.columns.get_loc(col)
        prev_col = filtered_data.columns[col_idx - 1]

        # Make sure we compare strings safely
        target_series = filtered_data[col].astype(str)

        condition = target_series.isin(['null', 'NotDEG'])
        filtered_data.loc[condition, prev_col] = 'null'


cols_to_drop = [
    col for col in filtered_data.columns
    if col.startswith("value_log2FC") or "adjusted" in col
]

filtered_data = filtered_data.drop(columns=cols_to_drop)

number_of_degs = len(filtered_data)
list_of_degs = filtered_data.iloc[:, 0].tolist()

#Export save
output_file_deseq_filter = 'filter_' + output_file_idsupdown
filtered_data.to_csv(output_file_deseq_filter, sep='\t', index=False, na_rep='null')

print(f"{output_file_deseq_filter} created, containing only DEGs log2FC (Result output 2)")


df = pd.read_csv(output_file_deseq_filter, sep='\t')

# Initialize dictionaries to store the count of positive and negative values
positive_counts = {}
negative_counts = {}

# Iterate through columns with the header "log2"
for column in df.columns:
    if "log2FC" in column:
        # Convert non-empty strings to numeric, treating empty strings as NaN
        numeric_values = pd.to_numeric(df[column], errors='coerce')

        # Calculate counts for positive and negative values, skipping NaN values
        positive_counts[column] = (numeric_values > 0).sum()
        negative_counts[column] = (numeric_values < 0).sum()

# Create a new DataFrame to store the positive and negative counts
count_df = pd.DataFrame({"Sample": positive_counts.keys(),
                         "UpRegulated": positive_counts.values(),
                         "DownRegulated": negative_counts.values()})

count_df["Total"] = count_df["UpRegulated"] + count_df["DownRegulated"]

# Remove the 'log2FC' prefix from the 'Sample' column
count_df["Sample"] = count_df["Sample"].str.replace("^log2FC_", "", regex=True)

# Print the count DataFrame
count_df.to_csv('degs_count.tsv', sep='\t', index=False)


print(f"degs_count.tsv created, containing the number of upregulated and downregulated genes per sample (Result output 3)")

# Read the TSV file into a dataframe
df = pd.read_csv(output_file_deseq_filter, delimiter='\t')

# Identify columns containing "log2FC"
log2FC_columns = [col for col in df.columns if 'log2FC' in col]

# Lists to store the output filenames
positive_filename_list = []
negative_filename_list = []

# Separate positive and negative values for each log2FC column
for col in log2FC_columns:
    col_name = col.replace('log2FC_', '')  # Remove the 'log2FC_' prefix
    positive_ids = df.iloc[:, 0][df[col] > 0]
    negative_ids = df.iloc[:, 0][df[col] < 0]

    # Create filenames for positive and negative IDs
    positive_filename = f'upregulated_{col_name}_seqID.txt'
    negative_filename = f'downregulated_{col_name}_seqID.txt'

    # Append filenames to their respective lists
    positive_filename_list.append(positive_filename)
    negative_filename_list.append(negative_filename)

    # Write positive IDs to the respective file
    positive_ids.to_csv(positive_filename, index=False)

    # Write negative IDs to the respective file
    negative_ids.to_csv(negative_filename, index=False)

print(f"Lists containing downregulated_ and upregulated_ seqIDs created")



data = pd.read_csv(args.go, sep='\t')

#if the GO annotation has more than 2 colunms (i.e. blast2go output table)
if data.shape[1] > 2:
    if 'SeqName' in data.columns and 'GO IDs' in data.columns:
        data = data[['SeqName', 'GO IDs']]
    else:
        raise ValueError("Expected columns 'SeqName' and 'GO IDs' not found in the input file.")

    # Extract only GO terms that match the specified category (e.g., P:GO:XXXXXXX)
    pattern = fr'{args.go_mode}:GO:\d{{7}}'
    data['GO IDs'] = data['GO IDs'].apply(lambda x: ' '.join(re.findall(pattern, str(x))))

    # drop rows with no matches
    data = data[data['GO IDs'] != '']

    # Remove the category prefix (e.g., P:)
    data['GO IDs'] = data['GO IDs'].str.replace(f'{args.go_mode}:', '', regex=False)

    # Drop the first column ('SeqName')
    go_series = data['GO IDs']

    # Replace spaces with newlines and explode into individual GO terms
    go_list = go_series.str.split().explode()

    # Drop duplicates and sort
    go_list = sorted(go_list.drop_duplicates())

    output_name = 'go_terms_' + args.go_mode + '_' + args.go

    # Output to a text file, one GO ID per line
    with open(output_name, "w") as f:
        f.write("GO IDs\n")
        for go_term in go_list:
            f.write(go_term + "\n")

# Only proceed if there are exactly 2 columns
if data.shape[1] == 2:
    data.columns = ['SeqName', 'GO IDs']  # rename for clarity

    # Convert the column to a single string for scanning
    sample_text = ' '.join(data['GO IDs'].astype(str).head(50))

    # Set pattern and post-processing rules
    if f"{args.go_mode}:GO" in sample_text:
        # Case 1: Match category-specific (e.g., P:GO:XXXXXXX)
        pattern = fr'{args.go_mode}:GO:\d{{7}}'
        add_prefix = False
        remove_category_prefix = True
    elif "GO:" in sample_text:
        # Case 2: Match generic GO:XXXXXXX
        pattern = r'GO:\d{7}'
        add_prefix = False
        remove_category_prefix = False
    else:
        # Case 3: Only 7-digit numbers, prefix needed
        pattern = r'\d{7}'
        add_prefix = True
        remove_category_prefix = False

    # Extract matching terms using regex
    data['GO IDs'] = data['GO IDs'].apply(lambda x: ' '.join(re.findall(pattern, str(x))))
    data = data[data['GO IDs'] != '']  # drop empty rows

    # Optionally remove the category prefix (e.g., P:)
    if remove_category_prefix:
        data['GO IDs'] = data['GO IDs'].str.replace(f'{args.go_mode}:', '', regex=False)

    # Explode into one GO term per row and drop duplicates
    go_list = data['GO IDs'].str.split().explode().drop_duplicates()

    # Optionally add "GO:" prefix
    if add_prefix:
        go_list = go_list.apply(lambda x: f'GO:{x}')

    # Sort and save
    go_list = sorted(go_list)

    # Create output filename
    output_basename = os.path.splitext(os.path.basename(args.go))[0]
    output_name = f"go_terms_{args.go_mode}_{output_basename}.txt"

    # Write to file
    with open(output_name, "w") as f:
        f.write("GO ID\n")
        for go_term in go_list:
            f.write(go_term + "\n")
if data.shape[1] < 2:
    raise ValueError("Expected columns 'SeqName' and 'GO IDs' not found in the input file.")



if args.go_mode == 'P':
    go_mode_two = 'BP'
elif args.go_mode == 'C':
    go_mode_two = 'CC'
elif args.go_mode == 'F':
    go_mode_two = 'MF'

print("Starting GO ancestor fetch using GO.DB R package")

# R script starts here as a Python multi-line string
r_script = f"""
library(GO.db)

# Read the GO terms from the text file
go_terms <- readLines("{output_name}")

# Convert the GO{go_mode_two}ANCESTOR object to a list
ancestor_list <- as.list(GO{go_mode_two}ANCESTOR)

# Create an empty data frame to store the results
result_df <- data.frame(GO_Term = character(), Ancestor_Terms = character(), stringsAsFactors = FALSE)

# Iterate over the GO terms and retrieve their ancestor terms
for (term in go_terms) {{
    term <- trimws(term)
    if (term %in% names(ancestor_list)) {{
        ancestor <- ancestor_list[[term]]
        if (length(ancestor) > 0) {{
            ancestor <- paste(ancestor, collapse = ",")
            result_df <- rbind(result_df, data.frame(
                GO_Term = term,
                Ancestor_Terms = paste(term, ancestor, sep = ","),
                stringsAsFactors = FALSE
            ))
        }}
    }}
}}

# Remove ", NA" from the Ancestor_Terms column
result_df$Ancestor_Terms <- gsub(",NA", "", result_df$Ancestor_Terms)
result_df$Ancestor_Terms <- gsub(",all", "", result_df$Ancestor_Terms)
# Write output
write.table(result_df, file = "go_ancestor_{output_name}.tsv", sep = "\\t", row.names = FALSE, quote = FALSE)
"""
go_ancestor_output = 'go_ancestor_'+ output_name + '.tsv'

# Save the R script to a temporary file
r_script_file = "run_ancestor_analysis.R"
with open(r_script_file, "w") as f:
    f.write(r_script)

# Run the R script
with open("run_ancestor_analysis_log.txt", "w") as log:
    subprocess.run(
        ["Rscript", r_script_file],
        check=True,
        stdout=log,
        stderr=log
    )

go_pattern = re.compile(r"GO:\d{7}")

def validate_go_ancestor_output(filepath):
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            cols = line.split("\t")  # change to "," or whitespace if needed
            if len(cols) < 2:
                continue

            go_terms = set(go_pattern.findall(cols[1]))

            if len(go_terms) >= 2:
                return True  # safety condition satisfied

    return False


if not validate_go_ancestor_output(go_ancestor_output):
    print(f"Error: {go_ancestor_output} does not contain ancestor terms, probably because R (v >=4) and GO.DB were not installed properly. Please, check run_ancestor_analysis_log.txt to see if there was any warning",
    file=sys.stderr
    )
    sys.exit(1)



print("Ancestor terms assigned to each genome GO term")

#RScript1 end

#This script will map the seqids to the go ancestor terms

data = pd.read_csv(args.go, sep='\t')

#if the GO annotation has more than 2 colunms (i.e. blast2go output table)
if data.shape[1] > 2:
    if 'SeqName' in data.columns and 'GO IDs' in data.columns:
        data = data[['SeqName', 'GO IDs']]
    else:
        raise ValueError("Expected columns 'SeqName' and 'GO IDs' not found in the input file.")

    # Extract only GO terms that match the specified category (e.g., P:GO:XXXXXXX)
    pattern = fr'{args.go_mode}:GO:\d{{7}}'
    data['GO IDs'] = data['GO IDs'].apply(lambda x: ' '.join(re.findall(pattern, str(x))))

# Only proceed if there are exactly 2 columns
if data.shape[1] == 2:
    data.columns = ['SeqName', 'GO IDs']  # rename for clarity

    # Convert the column to a single string for scanning
    sample_text = ' '.join(data['GO IDs'].astype(str).head(50))

    # Set pattern and post-processing rules
    if f"{args.go_mode}:GO" in sample_text:
        # Case 1: Match category-specific (e.g., P:GO:XXXXXXX)
        pattern = fr'{args.go_mode}:GO:\d{{7}}'
        add_prefix = False
        remove_category_prefix = True
    elif "GO:" in sample_text:
        # Case 2: Match generic GO:XXXXXXX
        pattern = r'GO:\d{7}'
        add_prefix = False
        remove_category_prefix = False
    else:
        # Case 3: Only 7-digit numbers, prefix needed
        pattern = r'\d{7}'
        add_prefix = True
        remove_category_prefix = False

    # Extract matching terms using regex
    data['GO IDs'] = data['GO IDs'].apply(lambda x: ' '.join(re.findall(pattern, str(x))))
    data = data[data['GO IDs'] != '']  # drop empty rows

    # Optionally remove the category prefix (e.g., P:)
    if remove_category_prefix:
        data['GO IDs'] = data['GO IDs'].str.replace(f'{args.go_mode}:', '', regex=False)


data2 = pd.read_csv('go_ancestor_' + output_name + '.tsv', sep='\t')
data2_dict = dict(zip(data2['GO_Term'], data2['Ancestor_Terms']))

# Step 1: Split GO terms (in case they are space-separated strings)
data['GO_ID_List'] = data['GO IDs'].str.split(' ')

# Step 2: Explode to get one row per GO term
data = data.explode('GO_ID_List')

# Step 3: Map each GO term to its ancestors
data['GO_Ancestors'] = data['GO_ID_List'].map(data2_dict)

# Optional: Drop the exploded GO term if you want only the original column names
data = data.drop(columns=['GO_ID_List'])

# Aggregate the mapped values back into a single list for each original row
result = data.groupby('SeqName')['GO_Ancestors'].agg(list).reset_index()

final_result = pd.DataFrame({
    'SeqName': result['SeqName'],
    'GO_Ancestors': result['GO_Ancestors'].apply(
        lambda x: ','.join(sorted(set(i for i in x if pd.notna(i))))
    )
})

# Save the result to a new file
final_result.to_csv(args.go + '_ancestor.tsv', sep='\t', index=None, na_rep='NAN')
go_annotation_ancestor = args.go + '_ancestor.tsv'

print("seqIDs mapped to ancestor GO terms")


positive_and_negative_filename_list = positive_filename_list + negative_filename_list


bp_ancestor_list = []
ancestor_file = f"{args.go}_ancestor.tsv"

# Read the ancestor file into memory
with open(ancestor_file, 'r') as f:
    ancestor_lines = f.readlines()

header = ancestor_lines[0]
data_lines = ancestor_lines[1:]

# Process each filename
for filename in positive_and_negative_filename_list:
    filename = filename.strip()

    # Read patterns from filename
    with open(filename, 'r') as pattern_file:
        patterns = [line.strip() for line in pattern_file if line.strip()]

    # Build the output file name
    output_file = f"{args.go}_ancestor_{filename}"
    bp_ancestor_list.append(output_file)

    # Filter the lines
    with open(output_file, 'w') as out_f:
        out_f.write(header)  # Write the header first
        for line in data_lines:
            if '\tGO:' not in line:
                continue
            if any(pattern in line for pattern in patterns):
                out_f.write(line)


#This script will replace GO terms from the second colunm of a input for the corresponding IDs from the genes

# Load the file
dfago = pd.read_csv(go_ancestor_output, sep='\t')
# Define GO term pattern
pattern = r'GO:\d{7}'
# Extract GO terms using regex
dfago['Ancestor_Terms'] = dfago['Ancestor_Terms'].apply(
    lambda x: ' '.join(re.findall(pattern, str(x)))
)
# Split into separate terms
dfago['Ancestor_Terms'] = dfago['Ancestor_Terms'].str.split(' ')
# Explode to get one term per row
dfago = dfago.explode('Ancestor_Terms')
# Remove duplicates and keep only the 'Ancestor_Terms' column
dfago = dfago[['Ancestor_Terms']].drop_duplicates()

# Optional: Sort the GO terms
dfago = dfago.sort_values(by='Ancestor_Terms')


# Read content from the list file and create a list of input file paths
input_file_paths = bp_ancestor_list
output_ancestor_final_list =[]
# Iterate through each input file path
for file1_path in input_file_paths:
    # Generate the output file name based on the input file name
    # Find the index where '_downregulated' starts
    index = file1_path.find('downregulated')

    # If not found, try '_upregulated'
    if index == -1:
        index = file1_path.find('upregulated')

    # Slice the string from the found index
    output_file_name = file1_path[index:] if index != -1 else file1_path  # fallback if neither found

    # Read content from file1 and create a dictionary mapping GO terms to sequences
    go_dict = {}
    with open(file1_path, 'r') as file1:
        for line in file1:
            seq_id, go_terms = line.strip().split('\t')
            go_terms_list = go_terms.split(',')
            for go_term in go_terms_list:
                if go_term not in go_dict:
                    go_dict[go_term] = []
                go_dict[go_term].append(seq_id)

    # Read content from file2 and create a list of GO terms
    go_terms_to_search = dfago['Ancestor_Terms'].dropna().unique().tolist()
    # Create a dictionary to store the results
    result_dict = {go_term: list(set(go_dict.get(go_term, []))) for go_term in go_terms_to_search}

    # Write the results to the output file with header
    output_file_path = f'ancestor_{output_file_name}'
    output_ancestor_final_list.append(output_file_path)
    lines_to_write = []
    lines_to_write.append("GO Term\t" + output_file_name + "\n")  # header
    for go_term, seq_ids in result_dict.items():
        lines_to_write.append(f"{go_term}\t{','.join(seq_ids)}\n")

    # Remove the second line (index 1) if it exists
    if len(lines_to_write) > 1:
        del lines_to_write[1]

    with open(output_file_path, 'w') as file3:
        file3.writelines(lines_to_write)


# Now to  Initialize with the first file
merged_df = None

for i, file_path in enumerate(output_ancestor_final_list):
    df = pd.read_csv(file_path, sep='\t')  # assuming tab-separated

    if i == 0:
        # Keep both columns for the first file
        merged_df = df
    else:
        # Only keep the second column and rename it to the filename (optional)
        second_col = df.columns[1]
        df = df[[second_col]]
        #df.columns = [file_path]  # rename column to the file name (or extract a label if preferred)

        # Concatenate side-by-side
        merged_df = pd.concat([merged_df, df], axis=1)


# This section is to remove the suffix from the header and order the colunms pairwise for upreg and downreg
original_columns = merged_df.columns.tolist()
first_col = original_columns[0]
other_cols = original_columns[1:]

# --- Step 1: Detect and remove common suffix ---
def common_suffix(strings):
    reversed_strings = [s[::-1] for s in strings]
    min_len = min(map(len, reversed_strings))
    suffix = ''
    for i in range(min_len):
        chars = set(s[i] for s in reversed_strings)
        if len(chars) == 1:
            suffix += chars.pop()
        else:
            break
    return suffix[::-1]

suffix = common_suffix(other_cols)
new_cols = [first_col] + [col.replace(suffix, '_seqID') for col in other_cols]

# Apply new column names
merged_df.columns = new_cols

# --- Step 2: Reorder to pair upregulated/downregulated columns ---
# Separate upregulated and downregulated columns
up_cols = [col for col in new_cols if col.startswith("upregulated_")]
down_cols = [col for col in new_cols if col.startswith("downregulated_")]

# Build new column order: start with first column, then interleave up/down
paired_cols = []
for up in up_cols:
    suffix = up.replace("upregulated_", "")
    down = f"downregulated_{suffix}"
    if down in down_cols:
        paired_cols.extend([up, down])
        down_cols.remove(down)
    else:
        paired_cols.append(up)

# Append any remaining unmatched downregulated columns
paired_cols.extend(down_cols)

# Final column order
final_order = [first_col] + paired_cols
merged_df = merged_df[final_order]

# Now merged_df has the first column (GO Term) once, and all second columns from the files
merged_df.to_csv("ancestor_" + output_file_idsupdown, sep='\t', index=False)
merged_go_IDs = "ancestor_" + output_file_idsupdown + "_GO_terms.tsv"
# Select only the first column by position
merged_df.iloc[:, [0]].to_csv(merged_go_IDs, sep='\t', index=False)

print("DEGs seqIDs assigned for each go term (ancestors and specifics)")

#get the GO names
r_script2 = f"""
# This R script will use a single colunm input with GO ids and output a two colunms file with the GO ids and GO term name
library(GO.db)

# Read the GO terms from the file (excluding the header), this is a list of GO ids
go_terms <- readLines("{merged_go_IDs}")[-1]

keys <- head(keys(GO.db))

result_df = select(GO.db, keys=go_terms, columns=c("TERM"),
      keytype="GOID")

# Write the result data frame to a file
write.table(result_df, file = "{merged_go_IDs}_ids_mapped.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
"""
go_ids_names = merged_go_IDs+'_ids_mapped.tsv'

# Save the R script to a temporary file
r_script_file2 = "get_go_term_by_id.R"
with open(r_script_file2, "w") as f:
    f.write(r_script2)
# Run the R script

with open("get_go_term_by_id_log.txt", "w") as log:
    subprocess.run(
        ["Rscript", r_script_file2],
        check=True,
        stdout=log,
        stderr=log
    )


#now map all the genome seqid to each go category

file1_path = go_annotation_ancestor
cleaned_file_path = file1_path + "_cleaned.tsv"

with open(file1_path, 'r') as infile, open(cleaned_file_path, 'w') as outfile:
    for line in infile:
        parts = line.strip().split('\t')
        if len(parts) == 2 and parts[1].strip():  # keep only lines with non-empty second column
            outfile.write(line)

# Now use the cleaned file
file1_path = cleaned_file_path  # update the path

# Iterate through the cleaned input file
output_file_name = args.go + '_ancestor_genome.tsv'
# Read content from file1 and create a dictionary mapping GO terms to sequences
go_dict = {}
with open(file1_path, 'r') as file1:
    for line in file1:
        seq_id, go_terms = line.strip().split('\t')
        go_terms_list = go_terms.split(',')
        for go_term in go_terms_list:
            if go_term not in go_dict:
                go_dict[go_term] = []
            go_dict[go_term].append(seq_id)

# Read content from file2 and create a list of GO terms
go_terms_to_search = dfago['Ancestor_Terms'].dropna().unique().tolist()
# Create a dictionary to store the results
result_dict = {go_term: list(set(go_dict.get(go_term, []))) for go_term in go_terms_to_search}

# Build all lines in memory first
lines_to_write = []
lines_to_write.append("GO Term\t" + output_file_name + "\n")  # header
for go_term, seq_ids in result_dict.items():
    lines_to_write.append(f"{go_term}\t{','.join(seq_ids)}\n")

# Remove the second line if it exists
if len(lines_to_write) > 1:
    del lines_to_write[1]

# Write final output
output_file_path = f'ancestor_{output_file_name}'
with open(output_file_path, 'w') as file3:
    file3.writelines(lines_to_write)


#this portion is to combine the GO names and genome seqID to the final dataframe
#Load the files
go_ids_names = pd.read_csv(go_ids_names, sep='\t')
output_file = pd.read_csv(output_file_path, sep='\t')

# Reconstruct the dataframe
# Insert second col of go_ids_names after first col of merged_df
# Insert second col of output_file after that
combined_df = pd.concat([
    merged_df.iloc[:, [0]],                  # 1st column of merged_df
    go_ids_names.iloc[:, [1]],               # 2nd column of go_ids_names
    output_file.iloc[:, [1]],                # 2nd column of output_file
    merged_df.iloc[:, 1:]                    # Remaining columns of merged_df
], axis=1)

# Remove duplicates within comma-delimited cells from 4th column onward
for col in combined_df.columns[3:]:
    combined_df[col] = combined_df[col].apply(
        lambda x: ','.join(sorted(set(str(x).split(',')))) if pd.notna(x) and x != 'null' else x
    )

combined_df.rename(columns={combined_df.columns[2]: 'Genome_seqid'}, inplace=True)


def add_count_columns(df, start_index=2):
    """
    For every column from start_index onward, create a duplicate with suffix '_count'
    where each cell contains the number of comma-separated terms from the original column.
    """
    count_columns = {}

    for i, col in enumerate(df.columns[start_index:], start=start_index):
        count_col_name = f"{col}_count"
        count_columns[count_col_name] = df[col].apply(
            lambda x: len(str(x).split(',')) if pd.notna(x) and x != 'null' and str(x).strip() != '' else 0
        )

    # Insert the new columns right after their corresponding originals
    for i, (col_name, count_series) in enumerate(count_columns.items()):
        original_col_index = df.columns.get_loc(col_name.replace('_count', ''))
        df.insert(original_col_index + 1, col_name, count_series)

    return df

combined_df = add_count_columns(combined_df)


def add_deg_sum_columns_adjacent(df):
    df = df.copy()
    # Extract all headers that end with '_seqID_count'
    count_cols = [col for col in df.columns if col.endswith('_seqID_count')]

    # Build a sample map of paired upregulated/downregulated columns
    pattern = r'(upregulated|downregulated)_(.+)_seqID_count'
    sample_map = {}

    for col in count_cols:
        match = re.match(pattern, col)
        if match:
            condition, sample = match.groups()
            if sample not in sample_map:
                sample_map[sample] = {}
            sample_map[sample][condition] = col

    for sample, pair in sample_map.items():
        up_col = pair.get('upregulated')
        down_col = pair.get('downregulated')
        if up_col and down_col:
            deg_col = f'DEG_{sample}_seqID_count'
            # Compute the new column
            df.loc[:, deg_col] = df[up_col].fillna(0).astype(int) + df[down_col].fillna(0).astype(int)

            # Find index of the later of the two columns
            idx_up = df.columns.get_loc(up_col)
            idx_down = df.columns.get_loc(down_col)
            insert_at = max(idx_up, idx_down) + 1

            # Reorder columns to insert DEG col right after the pair
            cols = list(df.columns)
            cols.remove(deg_col)
            cols.insert(insert_at, deg_col)
            df = df[cols]

    return df

combined_df = add_deg_sum_columns_adjacent(combined_df)

def add_deg_total_column(df):
    # Step 1: Select relevant columns starting from the 5th column
    relevant_cols = [col for col in df.columns[4:] if col.endswith('_seqID')]

    # Step 2: Apply per-row processing
    def count_unique_terms(row):
        terms = []
        for col in relevant_cols:
            cell = str(row[col])
            split_terms = re.split(r'[\t,]', cell)
            terms.extend(t.strip() for t in split_terms if t.strip())
        return len(set(terms))

    # Step 3: Apply and assign new column
    df['DEG_total'] = df.apply(count_unique_terms, axis=1)

    return df

combined_df = add_deg_total_column(combined_df)

def load_deg_list(treatment):
    up_file = f"upregulated_{treatment}_seqID.txt"
    down_file = f"downregulated_{treatment}_seqID.txt"

    def read_file(fname):
        with open(fname) as f:
            lines = f.read().strip().splitlines()
        return set(lines[1:])  # skip header

    degs = read_file(up_file) | read_file(down_file)
    return degs


#get deg list from filtered deseq2 file
df = pd.read_csv(output_file_deseq_filter, sep="\t")
degs_full_list = set(df.iloc[:, 0].dropna().astype(str))

max_count = combined_df['Genome_seqid_count'].max()

genome_seqids_set = set()

for seqids in combined_df.loc[
    combined_df['Genome_seqid_count'] == max_count,
    'Genome_seqid'
]:
    genome_seqids_set.update(seqids.split(','))


go_major = len(genome_seqids_set | degs_full_list)

bgsize = args.bgsize

def compute_fisher_chunk(args):
    chunk, col, count_col, seqid_col, degs, non_deg_genes, treatment, bgsize = args
    result = []

    number_of_degs = len(degs)

    for idx, row in chunk.iterrows():
        a = row[col]
        b = row[count_col] - a
        c = number_of_degs - a
        if bgsize is not None:
            d =  bgsize - a -b -c
        else:
            d = go_major - a -b -c

        contingency = [[a, b], [c, d]]

        _, pval = fisher_exact(contingency, alternative='greater')
        if a<5:
            pval = 1
        result.append({
            'index': idx,
            'pval': pval,
            'a': a,
            'b': b,
            'c': c,
            'd': d
        })
    return result


def add_fisher_pvals_parallel(df, threads=args.threads):
    count_col = 'Genome_seqid_count'
    deg_count_cols = [col for col in df.columns if col.startswith('DEG_') and col.endswith('_seqID_count')]

    # Ensure numeric
    df[count_col] = pd.to_numeric(df[count_col], errors='coerce').fillna(0).astype(int)
    for col in deg_count_cols:
        df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0).astype(int)

    genome_go_seqid = [col for col in df.columns if col.endswith("seqid")]
    if len(genome_go_seqid) != 1:
        raise ValueError(f"Expected exactly one column ending in 'seqid', found: {genome_go_seqid}")
    seqid_col = genome_go_seqid[0]

    contingency_rows = []

    for col in deg_count_cols:
        treatment = col.replace("DEG_", "").replace("_seqID_count", "")

        # Load DEGs for this treatment
        degs = load_deg_list(treatment)
        non_deg_genes = set(list_of_genes) - degs

        # Split dataframe for parallel
        #chunks = np.array_split(df, threads)
        indices = np.array_split(df.index, threads)
        chunks = [df.loc[idx] for idx in indices]
        args_list = [
            (chunk, col, count_col, seqid_col, degs, non_deg_genes, treatment, bgsize)
            for chunk in chunks
        ]

        with ProcessPoolExecutor(max_workers=threads) as executor:
            results = list(executor.map(compute_fisher_chunk, args_list))

        all_results = [row for chunk_result in results for row in chunk_result]
        all_results.sort(key=lambda x: x['index'])

        pvals = [r['pval'] for r in all_results]
        a_vals = [r['a'] for r in all_results]
        b_vals = [r['b'] for r in all_results]
        c_vals = [r['c'] for r in all_results]
        d_vals = [r['d'] for r in all_results]

        pval_col = col.replace('_seqID_count', '_f_test_p_value')
        insert_idx = df.columns.get_loc(col) + 1
        df.insert(insert_idx, pval_col, pvals)

        _, pvals_corrected, _, _ = multipletests(pvals, method='fdr_tsbh')
        fdr_col = col.replace('_seqID_count', '_f_test_FDR')
        df.insert(insert_idx + 1, fdr_col, pvals_corrected)
        # _, pvals_corrected, _, _ = multipletests(pvals, method='bonferroni')
        # bonferroni_col = col.replace('_seqID_count', '_f_test_bonferroni')
        # df.insert(insert_idx + 2, bonferroni_col, pvals_corrected)

        contingency_rows.extend([
            {
                'row_index': i,
                'treatment': treatment,
                'a': a_vals[i],
                'b': b_vals[i],
                'c': c_vals[i],
                'd': d_vals[i]
            }
            for i in range(len(a_vals))
        ])

    pd.DataFrame(contingency_rows).to_csv("contigency_variables.tsv", sep='\t', index=False)
    return df

print("GO Enrichment analyses using fisher exact test started")
combined_df = add_fisher_pvals_parallel(combined_df, threads=args.threads)
print("GO Enrichment analyses using fisher exact test complete")

combined_output = f"Complete_GO_{go_mode_two}_classification_{os.path.splitext(deg_input)[0]}"
combined_output_tsv = combined_output + '.tsv'
combined_output_xlsx = combined_output + '.xlsx'

combined_df = combined_df.sort_values(by='DEG_total', ascending=False) #sort dataframe by the last colunm before send to output

# Save to new TSV
combined_df.to_csv(combined_output_tsv, sep='\t', index=False)
# Save as Excel file
#combined_df.to_excel("combined_output.xlsx", index=False)
# Write with freeze pane at E2 (i.e., scroll starts at row 2, column E)
with pd.ExcelWriter(combined_output_xlsx, engine='xlsxwriter') as writer:
    combined_df.to_excel(writer, sheet_name='Sheet1', index=False)

    # Get the workbook and worksheet objects
    workbook  = writer.book
    worksheet = writer.sheets['Sheet1']

    # Freeze pane at E2 (row=1, col=4 because 0-indexed)
    worksheet.freeze_panes(1, 4)

    # Adjust column width
    worksheet.set_column(0, 0, 11)
    worksheet.set_column(1, 1, 18)
    worksheet.set_column(2, 2, 11)
    #worksheet.set_column(0, len(combined_df.columns) - 1, 13)



# Example: load your dataframe (replace with your actual file)
df = combined_df

# Select only columns ending with "_FDR"
fdr_cols = [col for col in df.columns if col.endswith('_FDR')]

# Filter rows where at least one _FDR column > 0.05
#filtered_df = df[df[fdr_cols].lt(0.05).any(axis=1)]
filtered_df = df[df[fdr_cols].lt(args.epvalue).any(axis=1)]

# Define which columns to keep:
#  - the first 4 columns
#  - columns ending with "_FDR"
#  - columns starting with "DEG_" and ending with "_count"
cols_to_keep = [
    col for i, col in enumerate(df.columns)
    if i < 4
    or col.endswith('_FDR')
    or (col.startswith('DEG_') and col.endswith('_count'))
]

# Keep only those columns (and remove duplicates, in case of overlap)
filtered_df = filtered_df[cols_to_keep]

def add_count_columns(df):
    """
    For every column ending with '_FDR':
      - Create a new column immediately after it named 'count_<FDR column name>'.
      - Then, for each row:
          - If the column before the new one (position -1) < 0.05,
            copy the value from the column position -2.
          - Else, put NaN.
    """
    df = df.copy()
    fdr_cols = [col for col in df.columns if col.endswith('_FDR')]
    new_columns = list(df.columns)

    # Step 1: Create new empty columns immediately after each _FDR column
    for fdr_col in fdr_cols:
        fdr_index = new_columns.index(fdr_col)
        new_col_name = f"count_{fdr_col}"
        df.insert(fdr_index + 1, new_col_name, np.nan)
        new_columns = list(df.columns)  # update column list after insert

    # Step 2: Fill new columns based on the logic
    count_cols = [col for col in df.columns if col.startswith('count_') and col.endswith('_FDR')]
    for count_col in count_cols:
        new_index = df.columns.get_loc(count_col)
        if new_index < 2:
            continue  # need at least two columns before
        prev_col = df.columns[new_index - 1]
        prev2_col = df.columns[new_index - 2]

        df[count_col] = np.where(
            #df[prev_col] < 0.05,
            df[prev_col] < args.epvalue,
            df[prev2_col],
            np.nan
        )

    return df

filtered_df = add_count_columns(filtered_df)


cols_to_keep = [
    col for i, col in enumerate(filtered_df.columns)
    if i < 4
    or (col.startswith('count_') and col.endswith('_FDR'))
]

# Keep only those columns (and remove duplicates, in case of overlap)
filtered_df = filtered_df[cols_to_keep]


with pd.ExcelWriter("enriched_only.xlsx", engine='xlsxwriter') as writer:
    filtered_df.to_excel(writer, sheet_name='Sheet1', index=False)

    # Get the workbook and worksheet objects
    workbook  = writer.book
    worksheet = writer.sheets['Sheet1']

    # Freeze pane at E2 (row=1, col=4 because 0-indexed)
    worksheet.freeze_panes(1, 4)

    # Adjust column width
    worksheet.set_column(0, 0, 11)
    worksheet.set_column(1, 1, 18)
    worksheet.set_column(2, 2, 11)
    worksheet.set_column(2, 2, 11)



print("GO complete classification complete and saved to " + combined_output_tsv + " , " + combined_output_xlsx + " and enriched_only.xlsx " + "(Result output 4, 5 and 6)")



#The next portion will extract the

def sanitize_filename(filename):
    # Sanitize filename by removing forbidden characters for filenames, which are abundant in GO TERM
    return re.sub(r'[\\/*?:"<>|]', '', filename)

def process_row(row_data, seq_columns, filter_file_lines, output_dir):
    index, row = row_data
    #raw_output_file = row[1]  # second column
    raw_output_file = row.iloc[1]  # second column
    output_filename = sanitize_filename(raw_output_file) + '.tsv'
    output_file = os.path.join(output_dir, output_filename)

    seq_data = row[seq_columns].astype(str)
    seq_data = seq_data[seq_data.str.strip() != ""]
    all_ids = ",".join(seq_data)
    unique_seq_ids = sorted(set(all_ids.split(",")))

    with open(output_file, 'w') as f:
        # Write header
        f.write(filter_file_lines[0] + '\n')
        # Write matching lines
        for line in filter_file_lines[1:]:
            if any(seq_id in line for seq_id in unique_seq_ids):
                f.write(line.strip() + '\n')

    return 1  # Return 1 for progress counting


def process_table_parallel(input_file, output_file_deseq_filter, max_workers=args.threads):
    # Output dir
    output_dir = 'degs_by_go_term'
    os.makedirs(output_dir, exist_ok=True)

    # Read input
    df = pd.read_csv(input_file, delimiter='\t', low_memory=False)
    total_rows = len(df)

    # Get seqID columns
    seq_columns = [col for col in df.columns[3:] if col.endswith('seqID')]

    # Read filter file once
    with open(output_file_deseq_filter, 'r') as f:
        filter_file_lines = [line.strip() for line in f.readlines()]

    # Prepare the function with constant arguments
    worker = partial(process_row, seq_columns=seq_columns, filter_file_lines=filter_file_lines, output_dir=output_dir)

    # Progress bar
    processed = 0
    bar_length = 40

    # Process in parallel
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        for _ in executor.map(worker, df.iterrows()):
            processed += 1
            progress = processed / total_rows
            filled_length = int(bar_length * progress)
            bar = '#' * filled_length + '-' * (bar_length - filled_length)
            sys.stdout.write(f'\rProgress: [{bar}] {progress * 100:.1f}% ({processed}/{total_rows})')
            sys.stdout.flush()



print("Creating GO categories files containing DEGs")
input_file = combined_output_tsv
filter_file = output_file_deseq_filter
process_table_parallel(input_file, filter_file)
print("\nAll GO categories files saved to degs_by_go_term (Result output 6)")

# Return to the original directory
os.chdir(original_dir)

# Create result directory
result_dir = output_dir + "_result"
os.makedirs(result_dir, exist_ok=True)


# Move the main output files
shutil.copy(os.path.join(output_dir, combined_output_xlsx), os.path.join(result_dir, combined_output_xlsx))
shutil.copy(os.path.join(output_dir, combined_output_tsv), os.path.join(result_dir, combined_output_tsv))
shutil.copy(os.path.join(output_dir, output_file_idsupdown), os.path.join(result_dir, output_file_idsupdown))
shutil.copy(os.path.join(output_dir, output_file_deseq_filter), os.path.join(result_dir, output_file_deseq_filter))
shutil.copy(os.path.join(output_dir, "degs_count.tsv"), os.path.join(result_dir, "degs_count.tsv"))
shutil.copy(os.path.join(output_dir, "enriched_only.xlsx"), os.path.join(result_dir, "enriched_only.xlsx"))

src = os.path.join(output_dir, "degs_by_go_term")
dst = os.path.join(result_dir, "degs_by_go_term")

subprocess.run(["mv", src, dst], check=True)
print(f"Result files saved in the folder {result_dir}")
