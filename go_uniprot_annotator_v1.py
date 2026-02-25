import aiohttp
import asyncio
import pandas as pd
from tqdm import tqdm
from io import StringIO
import argparse
import os
import re
import shutil

parser = argparse.ArgumentParser(
    description="Annotate preteomes go terms and provides the protein name for the best hits based on a blast output file against Uniprot and SwissProt databases."
)
parser.add_argument(
    "-bu", "--blastuni", required=True, help="Path to the Blast output file against Uniprot database (TSV format)"
)
parser.add_argument(
    "-bs", "--blastswp", required=False, help="Path to the Blast output file against SwissProt database (TSV format)"
)

parser.add_argument(
    "-t", "--top", required=False, default=200, type=int, help="Number of top blast hits to search GO terms for each gene. Higher numbers can increase the number of annotated terms but increase the uniprot server search significantly"
)

args = parser.parse_args()

# Create a working directory with incremental numbering
base_name = "go_uniprot_annotator"
counter = 1
while True:
    output_dir = f"{base_name}_{counter}"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        break
    counter += 1

print('The intermediary files for this run will be saved in the folder ' + output_dir)

# Copy input files into the new folder
shutil.copy(args.blastuni, os.path.join(output_dir, os.path.basename(args.blastuni)))
if args.blastswp is not None:
    shutil.copy(args.blastswp, os.path.join(output_dir, os.path.basename(args.blastswp)))

# Save current working directory
original_dir = os.getcwd()
# Change to the new directory
os.chdir(output_dir)
def las():
    # Read the Diamond output file
    input_file = args.blastuni
    output_top20 = input_file.replace('.tsv', f'_best_{args.top}.tsv')
    output_top1 = input_file.replace('.tsv', '_best_1.tsv')
    output_collapsed = output_top20.replace(f'_best_{args.top}.tsv', f'_collapsed_{args.top}.tsv')
    output_collapsed1 = output_top1.replace('_best_1.tsv', '_collapsed_top1.tsv')


    # Load input file (no header, tab-separated)
    df = pd.read_csv(input_file, sep='\t', header=None)

    # Sort and select top hits per query based on column 12 (index 11)
    df_sorted = df.sort_values(by=[0, 11], ascending=[True, False])
    top20 = df_sorted.groupby(0, group_keys=False).head(args.top)
    top1 = df_sorted.groupby(0, group_keys=False).head(1)

    # Save top results (raw) to file
    top20.to_csv(output_top20, sep='\t', index=False, header=False)
    top1.to_csv(output_top1, sep='\t', index=False, header=False)

    # Extract columns 0 and 1 (query ID and target)
    df_top2 = top20[[0, 1]].copy()
    df_top1 = top1[[0, 1]].copy()

    # Clean column 1: remove 'tr|' prefix and everything after second '|'
    df_top2[1] = df_top2[1].str.replace(r'^.*?\|([^|]+)\|.*$', r'\1', regex=True)
    df_top1[1] = df_top1[1].str.replace(r'^.*?\|([^|]+)\|.*$', r'\1', regex=True)

    # Group by query ID (column 0), aggregate unique targets in column 1
    collapsed = df_top2.groupby(0)[1].apply(lambda x: ','.join(sorted(set(x)))).reset_index()
    collapsed1 = df_top1.groupby(0)[1].apply(lambda x: ','.join(sorted(set(x)))).reset_index()

    # Save final collapsed output
    collapsed.to_csv(output_collapsed, sep='\t', index=False, header=False)
    collapsed1.to_csv(output_collapsed1, sep='\t', index=False, header=False)

    print(f" Top {args.top} results written to: {output_top20}")
    print(f" Top 1 results written to: {output_top1}")
    print(f" Collapsed top {args.top} written to: {output_collapsed}")
    print(f" Collapsed top 1 written to: {output_collapsed1}")

    # Load the collapsed file (no header)
    collapsed_file = output_collapsed
    collapsed_file1 = output_collapsed1
    output_ids_file = f"uniprot_ids_best_{args.top}.txt"
    output_ids_file1 = "uniprot_ids_best_1.txt"

    # Step 1: Read the second column (UniProt IDs)
    df = pd.read_csv(collapsed_file, sep='\t', header=None, usecols=[1])
    df1 = pd.read_csv(collapsed_file1, sep='\t', header=None, usecols=[1])

    # Step 2: Split by comma and flatten the list
    all_ids = df[1].str.split(',').explode()
    all_ids1 = df1[1].str.split(',').explode()

    # Step 3: Remove duplicates and sort (optional)
    unique_ids = sorted(all_ids.drop_duplicates())
    unique_ids1 = sorted(all_ids1.drop_duplicates())

    # Step 4: Write to output file (one ID per line)
    with open(output_ids_file, "w") as f:
        for uid in unique_ids:
            f.write(uid + "\n")

    # Step 4: Write to output file (one ID per line)
    with open(output_ids_file1, "w") as f:
        for uid in unique_ids1:
            f.write(uid + "\n")

    print(f" Unique UniProt IDs from best {args.top} hits written to: {output_ids_file}")
    print(f" Unique UniProt IDs from best hit written to: {output_ids_file1}")
    return input_file, df_top1, collapsed

#keep df_top1 and collapsed to use later in the script, and clean all other variables
input_file, df_top1, collapsed = las()
#Now to map the uniprot ids to go terms and protein names using uniprot API

INPUT_FILE = f"uniprot_ids_best_{args.top}.txt"
OUTPUT_FILE = "uniprot_annotation_go_terms.tsv"
BATCH_SIZE = 100000
FIELDS = "accession%2Cgo_p%2Cgo_c%2Cgo_f&format=tsv"


def chunk_list(lst, size):
    for i in range(0, len(lst), size):
        yield lst[i:i + size]

async def safe_json(resp):
    try:
        return await resp.json()
    except Exception:
        text = await resp.text()
        print(f"Failed to parse JSON. Response text:\n{text}")
        return {}

async def submit_id_mapping(session, ids):
    url = f"https://rest.uniprot.org/idmapping/run"
    data = {
        "from": "UniProtKB_AC-ID",
        "to": "UniProtKB",
        "ids": ",".join(ids)
    }
    async with session.post(url, data=data) as resp:
        if resp.status != 200:
            raise Exception(f"Submit failed with status {resp.status}")
        json_data = await safe_json(resp)
        job_id = json_data.get("jobId")
        if not job_id:
            raise Exception(f"Missing jobId in response: {json_data}")
        return job_id

async def wait_for_job(session, job_id):
    await asyncio.sleep(5)  # wait before first status check
    url = f"https://rest.uniprot.org/idmapping/status/{job_id}"
    while True:
        async with session.get(url) as resp:
            if resp.status != 200:
                raise Exception(f"Status check failed with {resp.status}")
            json_data = await safe_json(resp)
            status = json_data.get("jobStatus")
            if status in [None, "FINISHED"]:
                return
            elif status in ("RUNNING", "NEW"):
                await asyncio.sleep(5)
            else:
                raise Exception(f"Unexpected job status: {json_data}")

async def fetch_results_tsv(session, job_id):
    url = f"https://rest.uniprot.org/idmapping/uniprotkb/results/stream/{job_id}?fields={FIELDS}"
    async with session.get(url) as resp:
        if resp.status != 200:
            raise Exception(f"Failed to download results with status {resp.status}")
        text_data = await resp.text()
        return pd.read_csv(StringIO(text_data), sep='\t')

def should_skip_query():
    if os.path.isfile(OUTPUT_FILE):
        print(f"{OUTPUT_FILE} already present in the folder.\n "
              f"Please, exclude or rename {OUTPUT_FILE} to process this query again.")
        return True
    return False

async def process_batch(session, batch, semaphore, batch_idx):
    async with semaphore:  # UniProt concurrency gate
        while True:
            retries = 10
            for attempt in range(retries):
                try:
                    #print(f"Batch {batch_idx}: attempt {attempt + 1}/{retries}")
                    job_id = await submit_id_mapping(session, batch)
                    await wait_for_job(session, job_id)
                    df = await fetch_results_tsv(session, job_id)
                    return df  # success

                except (aiohttp.ClientError, asyncio.TimeoutError) as e:
                    print(f"Batch {batch_idx}: network error → {e}")
                    await asyncio.sleep(5)

                except Exception as e:
                    print(f"Batch {batch_idx}: error → {e}")
                    await asyncio.sleep(5)

            print(f"Batch {batch_idx}: resubmitting after cooldown")
            await asyncio.sleep(30)


async def main():
    if should_skip_query():
        return

    with open(INPUT_FILE) as f:
        ids = [line.strip() for line in f if line.strip()]

    batches = list(chunk_list(ids, BATCH_SIZE))

    MAX_CONCURRENT_JOBS = 10  #  number of batcher to send at the same time
    semaphore = asyncio.Semaphore(MAX_CONCURRENT_JOBS)

    async with aiohttp.ClientSession() as session:
        tasks = [
            process_batch(session, batch, semaphore, i)
            for i, batch in enumerate(batches, start=1)
        ]

        all_dfs = []
        for coro in tqdm(
            asyncio.as_completed(tasks),
            total=len(tasks),
            desc="Processing UniProt batches"
        ):
            df = await coro
            all_dfs.append(df)


    if all_dfs:
        final_df = pd.concat(all_dfs)
        final_df.to_csv(OUTPUT_FILE, sep='\t', index=False)
        print(f" Done. Results written to {OUTPUT_FILE}")
    else:
        print("No results to write.")

if __name__ == "__main__":
    asyncio.run(main())


INPUT_FILE = "uniprot_ids_best_1.txt"
OUTPUT_FILE = "uniprot_annotation_protein_name_best_hit.tsv"
BATCH_SIZE = 100000
FIELDS = "accession%2Cprotein_name&format=tsv"

if __name__ == "__main__":
    asyncio.run(main())


# Determine long GO term mapping
go_mode_map = {
    "P": "biological process",
    "F": "molecular function",
    "C": "cellular component"
}

final_df = pd.read_csv("uniprot_annotation_go_terms.tsv", sep='\t')

new_collapses_output_names = []  # create empty list BEFORE the loop

for short_mode, go_mode_long in go_mode_map.items():
    # Find the correct GO column name using a case-insensitive match
    pattern = re.compile(re.escape(go_mode_long), flags=re.IGNORECASE)
    matches = [col for col in final_df.columns if pattern.search(col)]
    if not matches:
        print(f"No GO column found matching '{go_mode_long}', skipping.")
        continue
    go_column = matches[0]

    # Create new dataframe with first column and extracted GO terms
    df_go_mode = pd.DataFrame({
        "db_id": final_df.iloc[:, 0],
        "go_term_" + go_mode_long.replace(" ", "_"): final_df[go_column].astype(str).str.extractall(r"(GO:\d{7})")[0]
                      .groupby(level=0).agg(",".join)
    })

    # Create dictionary from df_go_mode
    go_dict = dict(zip(df_go_mode['db_id'], df_go_mode.iloc[:, 1]))

    # Function to replace UniProt IDs with GO terms and remove duplicates
    def replace_ids_with_go(uniprot_ids_str):
        uniprot_ids = uniprot_ids_str.split(",")
        go_terms = []
        for uid in uniprot_ids:
            terms = go_dict.get(uid)
            if isinstance(terms, str) and terms:
                go_terms.extend(terms.split(","))
        unique_terms = sorted(set(go_terms))
        return ",".join(unique_terms)

    # Apply the transformation to second column of collapsed
    new_collapsed = collapsed.copy()
    new_collapsed.iloc[:, 1] = collapsed.iloc[:, 1].apply(replace_ids_with_go)
    new_collapsed.columns = ['seq_id', f"go_term_{go_mode_long.replace(' ', '_')}"]

    # Write output per GO type
    new_collapses_output = input_file.replace('.tsv', f"_GO_{go_mode_long.replace(' ', '_')}.tsv")
    new_collapses_output_names.append(new_collapses_output)
    new_collapsed.to_csv(new_collapses_output, sep='\t', index=False)

    print(f"Saved {go_mode_long} results to {new_collapses_output}")

# read protein annotation
df_protname = pd.read_csv("uniprot_annotation_protein_name_best_hit.tsv", sep='\t')

# make dictionary col1 -> col3
df_protname_dic = dict(zip(df_protname.iloc[:, 0], df_protname.iloc[:, 2]))

# duplicate second column into a new third column
df_top1.insert(2, "col3", df_top1.iloc[:, 1])

# replace values in the new column using the dictionary
df_top1.iloc[:, 2] = df_top1.iloc[:, 2].map(df_protname_dic).fillna(df_top1.iloc[:, 2])

#rename the columns
df_top1.columns = ["seqID", "UniProtID", "UniProtFunction"]

df_top1.to_csv("uniprot_protein_annotation.tsv", sep='\t', index=False)

if args.blastswp is not None:
    input_file_swp = args.blastswp
    output_top1_swp = input_file_swp.replace('.tsv', '_best_1.tsv')
    output_collapsed2 = output_top1_swp.replace('_best_1.tsv', '_collapsed_top1.tsv')

    df_swp = pd.read_csv(input_file_swp, sep='\t', header=None)
    df_sorted_swp = df_swp.sort_values(by=[0, 11], ascending=[True, False])
    top1_swp = df_sorted_swp.groupby(0, group_keys=False).head(1)

    top1_swp.to_csv(output_top1_swp, sep='\t', index=False, header=False)
    df_top1_swp = top1_swp[[0, 1]].copy()
    df_top1_swp[1] = df_top1_swp[1].str.replace(r'^.*?\|([^|]+)\|.*$', r'\1', regex=True)
    collapsed2 = df_top1_swp.groupby(0)[1].apply(lambda x: ','.join(sorted(set(x)))).reset_index()
    collapsed2.to_csv(output_collapsed2, sep='\t', index=False, header=False)

    collapsed_file2 = output_collapsed2
    output_ids_file1_swp = "swissprot_ids_best_1.txt"
    df1_swp = pd.read_csv(collapsed_file2, sep='\t', header=None, usecols=[1])
    all_ids1_swp = df1_swp[1].str.split(',').explode()
    unique_ids1_swp = sorted(all_ids1_swp.drop_duplicates())

    with open(output_ids_file1_swp, "w") as f:
        for uid in unique_ids1_swp:
            f.write(uid + "\n")


    INPUT_FILE = "swissprot_ids_best_1.txt"
    OUTPUT_FILE = "swissprot_annotation_protein_name_best_hit.tsv"
    BATCH_SIZE = 100000
    FIELDS = "accession%2Cprotein_name&format=tsv"

    if __name__ == "__main__":
        asyncio.run(main())

    df_protname_swp = pd.read_csv("swissprot_annotation_protein_name_best_hit.tsv", sep='\t')
    df_protname_swp_dic = dict(zip(df_protname_swp.iloc[:, 0], df_protname_swp.iloc[:, 2]))
    df_top1_swp.insert(2, "col3", df_top1_swp.iloc[:, 1])
    df_top1_swp.iloc[:, 2] = df_top1_swp.iloc[:, 2].map(df_protname_swp_dic).fillna(df_top1_swp.iloc[:, 2])
    df_top1_swp.columns = ["seqID", "SwissProtID", "SwissProtFunction"]
    df_top1_swp.to_csv("swissprot_protein_annotation.tsv", sep='\t', index=False)

    # Load the two tables
    df_uniprot = pd.read_csv("uniprot_protein_annotation.tsv", sep="\t", dtype=str)
    df_swissprot = pd.read_csv("swissprot_protein_annotation.tsv", sep="\t", dtype=str)

    # Get the first column name (key) from each
    key_uni = df_uniprot.columns[0]
    key_swi = df_swissprot.columns[0]

    # Align the key column names if needed
    if key_uni != key_swi:
        df_swissprot = df_swissprot.rename(columns={key_swi: key_uni})

    # Merge tables (outer join keeps all IDs from both)
    df_merged = pd.merge(df_uniprot, df_swissprot, on=key_uni, how="outer", suffixes=("_uniprot", "_swissprot"))

    # Save result
    df_merged.to_csv("uniprot_swissprot_protein_annotation.tsv", sep="\t", index=False, na_rep='no_hit')

# Return to the original directory
os.chdir(original_dir)

# Create result directory
result_dir = output_dir + "_result"
os.makedirs(result_dir, exist_ok=True)

# Move the main output files
for filename in new_collapses_output_names:
    shutil.copy(os.path.join(output_dir, filename), os.path.join(result_dir, filename))

if args.blastswp is not None:
    shutil.copy(os.path.join(output_dir, "uniprot_swissprot_protein_annotation.tsv"), os.path.join(result_dir, "uniprot_swissprot_protein_annotation.tsv"))
else:
    shutil.copy(os.path.join(output_dir, "uniprot_protein_annotation.tsv"), os.path.join(result_dir, "uniprot_protein_annotation.tsv"))

print(f"GO Annotation and best hit function files saved to {result_dir}.")
