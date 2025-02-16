#!/usr/bin/env python3
"""
Combined Workflow: Data Processing, BEAST XML Update, BEAST Run, and TreeAnnotator

This script does the following:
1. Processes variant folders to extract sequences and metadata.
2. Combines sequences into a single FASTA (with metadata saved as TSV).
3. Runs MAFFT to align the combined FASTA.
4. Removes the top (reference) sequence from the alignment.
5. Reads the BEAST XML template from the data folder and replaces the taxa and alignment blocks
   with new blocks based on the aligned FASTA.
6. Also updates output file names in the XML (log, trees, checkpoint) to use the reference variant name.
7. Saves all outputs in a project folder:
       projects/[reference_variant]/
           (combined FASTA, TSV, aligned FASTA files, etc.)
           beast/[reference_variant].xml  (modified BEAST XML file)
8. Runs BEAST with its working directory set to the beast folder so that BEAST output files are generated there.
9. Finally, applies TreeAnnotator with "-heights ca" to the trees file and saves the result as mcc.tree.
"""

import os
import csv
import subprocess
from Bio import SeqIO

def run_data_processing():
    # ---------------------------
    # Configuration and Setup
    # ---------------------------
    variants_dir = 'variants'             # Folder containing variant subfolders
    reference_variant = 'SARS-CoV-1'        # Folder name used as reference

    # Use a project folder for all outputs:
    project_dir = os.path.join("/home/conor/covid_dendrogram_comparison/projects", reference_variant)
    os.makedirs(project_dir, exist_ok=True)

    # Header for combined TSV file.
    new_header = ["AccessionID", "CollectionDate", "Location", "Lineage"]

    combined_fasta_records = []  # To store combined sequence records
    combined_tsv_rows = []       # To store combined metadata rows

    # ---------------------------
    # 1. Extract the reference sequence
    # ---------------------------
    ref_folder_path = os.path.join(variants_dir, reference_variant)
    if not os.path.isdir(ref_folder_path):
        raise FileNotFoundError(f"Reference folder not found: {ref_folder_path}")

    # Look for the first FASTA file in the reference folder.
    ref_fasta = None
    for fname in os.listdir(ref_folder_path):
        if fname.lower().endswith((".fasta", ".fa")):
            ref_fasta = os.path.join(ref_folder_path, fname)
            break

    if not ref_fasta or not os.path.isfile(ref_fasta):
        raise FileNotFoundError(f"No FASTA found in reference folder {reference_variant} ({ref_folder_path}).")

    # Parse the first sequence from the reference FASTA.
    ref_records = list(SeqIO.parse(ref_fasta, "fasta"))
    if not ref_records:
        raise ValueError(f"Reference FASTA {ref_fasta} is empty.")
    reference_record = ref_records[0]
    print(f"Reference record is {reference_record.id} from {ref_fasta}")

    # ---------------------------
    # 2. Process each variant folder (including the reference folder)
    # ---------------------------
    for subfolder in os.listdir(variants_dir):
        subfolder_path = os.path.join(variants_dir, subfolder)
        if not os.path.isdir(subfolder_path):
            continue

        # Example condition: if reference is Wuhan-Hu-1, skip SARS-CoV-1.
        if reference_variant == "Wuhan-Hu-1" and subfolder == "SARS-CoV-1":
            print("Skipping SARS-CoV-1 because reference is Wuhan-Hu-1.")
            continue

        # Identify FASTA and TSV (or CSV) files in the subfolder.
        fasta_file = None
        tsv_file = None
        csv_file = None
        for fname in os.listdir(subfolder_path):
            if fname.lower().endswith((".fasta", ".fa")):
                fasta_file = os.path.join(subfolder_path, fname)
            elif fname.endswith('.csv'):
                csv_file = os.path.join(subfolder_path, fname)
            elif fname.endswith('.tsv'):
                tsv_file = os.path.join(subfolder_path, fname)

        # If only CSV exists, convert it to TSV.
        if tsv_file is None and csv_file is not None:
            base = os.path.splitext(csv_file)[0]
            tsv_file = base + ".tsv"
            with open(csv_file, 'r', newline='') as infile, open(tsv_file, 'w', newline='') as outfile:
                reader = csv.reader(infile, delimiter=',')
                writer = csv.writer(outfile, delimiter='\t')
                for row in reader:
                    writer.writerow(row)
            print(f"Converted {csv_file} to {tsv_file}")

        if not fasta_file or not tsv_file:
            print(f"Skipping {subfolder_path}: missing FASTA or TSV file.")
            continue

        # Read metadata from TSV.
        with open(tsv_file, 'r', newline='') as tsv_handle:
            reader = csv.DictReader(tsv_handle, delimiter='\t')
            tsv_rows = list(reader)

        if not tsv_rows:
            print(f"No data in {tsv_file}; skipping.")
            continue

        # Filter rows by "Lineage" from the first row, then take the top 5.
        top5_rows = tsv_rows[:5]
        if not top5_rows:
            print(f"No rows found in {tsv_file}; skipping.")
            continue

        # Map original accession IDs to new IDs including the subfolder name.
        mapping_dict = {}
        for row in top5_rows:
            orig_acc = row["Accession ID"]
            new_acc = f"{subfolder}_{orig_acc}"
            reduced = {
                "AccessionID": new_acc,
                "CollectionDate": row["Collection date"],
                "Location": row["Location"],
                "Lineage": subfolder
            }
            combined_tsv_rows.append(reduced)
            mapping_dict[orig_acc] = new_acc

        # Process the FASTA: rename records that match the original accession.
        filtered_records = []
        for record in SeqIO.parse(fasta_file, "fasta"):
            for orig_acc, new_acc in mapping_dict.items():
                if orig_acc in record.id:
                    record.id = new_acc
                    record.name = new_acc
                    record.description = new_acc
                    filtered_records.append(record)
                    break

        combined_fasta_records.extend(filtered_records)

    # ---------------------------
    # 3. Insert the reference at the top
    # ---------------------------
    combined_fasta_records.insert(0, reference_record)

    # ---------------------------
    # 4. Write combined FASTA and TSV outputs in the project folder
    # ---------------------------
    prefix = reference_variant
    combined_fasta_filename = os.path.join(project_dir, f"{prefix}_combined.fasta")
    combined_tsv_filename = os.path.join(project_dir, f"{prefix}_combined.tsv")

    # Write the combined FASTA file.
    SeqIO.write(combined_fasta_records, combined_fasta_filename, "fasta")
    print(f"Combined FASTA saved to {combined_fasta_filename}")

    # Write the combined TSV file.
    with open(combined_tsv_filename, "w", newline='') as out_handle:
        writer = csv.DictWriter(out_handle, fieldnames=new_header, delimiter='\t')
        writer.writeheader()
        for row in combined_tsv_rows:
            writer.writerow(row)
    print(f"Combined TSV saved to {combined_tsv_filename}")

    # ---------------------------
    # 5. Align sequences using MAFFT
    # ---------------------------
    aligned_fasta_filename = os.path.join(project_dir, f"{prefix}_combined_aligned.fasta")
    mafft_command = ["mafft", "--auto", combined_fasta_filename]
    try:
        with open(aligned_fasta_filename, "w") as aligned_handle:
            subprocess.run(mafft_command, stdout=aligned_handle, check=True)
        print(f"Aligned FASTA saved to {aligned_fasta_filename}")
    except Exception as e:
        print("Error running MAFFT:", e)
        raise

    # ---------------------------
    # 6. Remove the reference sequence from the alignment
    # ---------------------------
    final_aligned_no_ref_filename = os.path.join(project_dir, f"{prefix}_combined_aligned_no_ref.fasta")
    aligned_records = list(SeqIO.parse(aligned_fasta_filename, "fasta"))
    if aligned_records:
        # Remove the first record (reference).
        trimmed_aligned_records = aligned_records[1:]
        SeqIO.write(trimmed_aligned_records, final_aligned_no_ref_filename, "fasta")
        print(f"Final aligned FASTA (without the top reference) saved to {final_aligned_no_ref_filename}")
    else:
        print(f"No records found in {aligned_fasta_filename} to remove reference from.")

    # Return the final aligned FASTA path, prefix, and project directory.
    return final_aligned_no_ref_filename, prefix, project_dir

def run_beast_xml_modification(aligned_fasta_path, prefix, project_dir):
    """
    Update the BEAST XML template by replacing the taxa and alignment blocks
    with new blocks based on the aligned FASTA. Also update the output file names
    in the XML (log, trees, checkpoint) to use the reference variant name.
    The modified BEAST XML is saved in:
         projects/[reference_variant]/beast/[reference_variant].xml
    """
    # ---------------------------
    # Configuration for BEAST XML
    # ---------------------------
    # Read the BEAST template from the data folder.
    TEMPLATE_PATH = os.path.join("data", "beast_template.xml")
    # Create the output folder for BEAST XML inside the project folder.
    beast_dir = os.path.join(project_dir, "beast")
    os.makedirs(beast_dir, exist_ok=True)
    OUTPUT_XML = os.path.join(beast_dir, f"{prefix}.xml")

    # ---------------------------
    # 1. Read the BEAST XML template
    # ---------------------------
    if not os.path.exists(TEMPLATE_PATH):
        raise FileNotFoundError(f"Template not found: {TEMPLATE_PATH}")
    with open(TEMPLATE_PATH, "r") as f:
        template_str = f.read()

    # ---------------------------
    # 2. Remove old taxa/alignment blocks from the template
    # ---------------------------
    start_marker = "<!-- The list of taxa to be analysed"
    end_marker = "</alignment>"
    start_idx = template_str.find(start_marker)
    if start_idx == -1:
        raise ValueError(f"Could not find start marker in template:\n{start_marker}")
    end_idx = template_str.find(end_marker)
    if end_idx == -1:
        raise ValueError(f"Could not find end marker ({end_marker}) in template.")
    end_idx += len(end_marker)  # Include the end marker

    # Split the template into parts outside the taxa/alignment block.
    head_part = template_str[:start_idx]
    tail_part = template_str[end_idx:]

    # ---------------------------
    # 3. Read the aligned FASTA and group sequences by lineage
    # ---------------------------
    if not os.path.exists(aligned_fasta_path):
        raise FileNotFoundError(f"Aligned FASTA not found: {aligned_fasta_path}")
    
    lineage_dict = {}  # Group sequences by lineage.
    for record in SeqIO.parse(aligned_fasta_path, "fasta"):
        # Assume IDs like "Alpha_EPI_12345"; extract "Alpha" as the lineage.
        lineage = record.id.split("_", 1)[0]
        if lineage not in lineage_dict:
            lineage_dict[lineage] = []
        lineage_dict[lineage].append(record)

    # ---------------------------
    # 4. Build new <taxa> blocks
    # ---------------------------
    taxa_master_lines = ['\t<taxa id="taxa">']
    lineage_blocks = []
    for lineage, records in lineage_dict.items():
        block_lines = [f'\t<taxa id="{lineage}">']
        for rec in records:
            taxa_master_lines.append(f'\t\t<taxon id="{rec.id}"/>')
            block_lines.append(f'\t\t<taxon idref="{rec.id}"/>')
        block_lines.append('\t</taxa>')
        lineage_blocks.append("\n".join(block_lines))
    taxa_master_lines.append('\t</taxa>')
    taxa_block = "\n".join(taxa_master_lines) + "\n\n" + "\n\n".join(lineage_blocks) + "\n"

    # ---------------------------
    # 5. Build new <alignment> block
    # ---------------------------
    alignment_lines = ['\t<alignment id="alignment" dataType="nucleotide">']
    for lineage, records in lineage_dict.items():
        for rec in records:
            alignment_lines.append('\t\t<sequence>')
            alignment_lines.append(f'\t\t\t<taxon idref="{rec.id}"/>')
            alignment_lines.append(f'\t\t\t{str(rec.seq)}')
            alignment_lines.append('\t\t</sequence>')
    alignment_lines.append('\t</alignment>')
    alignment_block = "\n".join(alignment_lines)

    # ---------------------------
    # 6. Combine new blocks with the unchanged parts of the template
    # ---------------------------
    new_mid_section = (
        "\n\t<!-- The list of taxa to be analysed (replaced by script) -->\n"
        + taxa_block +
        "\n" +
        alignment_block +
        "\n\n"
    )
    modified_xml = head_part + new_mid_section + tail_part

    # ---------------------------
    # 7. Update output file names in the XML to use the reference variant name.
    # ---------------------------
    modified_xml = modified_xml.replace("sars_cov_2_cpass.log", f"{prefix}.log")
    modified_xml = modified_xml.replace("sars_cov_2_cpass.trees", f"{prefix}.trees")
    modified_xml = modified_xml.replace("sars_cov_2_cpass.chkpt", f"{prefix}.chkpt")

    # ---------------------------
    # 8. Write out the modified BEAST XML file
    # ---------------------------
    with open(OUTPUT_XML, "w") as f:
        f.write(modified_xml)
    print(f"Created updated BEAST XML with new <taxa>, <alignment>, and output file names: {OUTPUT_XML}")

    # Return the path to the modified XML.
    return OUTPUT_XML

def run_beast_analysis(xml_path):
    """
    Run BEAST on the provided XML file using the specified BEAST executable.
    The working directory is set to the XML file's folder (i.e. the beast folder)
    so that output files (log, trees, chkpt) are generated there.
    """
    BEAST_EXECUTABLE = "/home/conor/tools/BEASTv10.5.0/bin/beast"
    beast_command = [BEAST_EXECUTABLE, "-overwrite", xml_path]
    work_dir = os.path.dirname(xml_path)
    print(f"Running BEAST analysis on {xml_path} in folder {work_dir} using {BEAST_EXECUTABLE}...")
    try:
        subprocess.run(beast_command, check=True, cwd=work_dir)
        print("BEAST analysis completed successfully.")
    except Exception as e:
        print("Error running BEAST analysis:", e)
        raise

def run_treeannotator(prefix, project_dir):
    """
    After BEAST finishes, run TreeAnnotator with the command:
         treeannotator -heights ca [prefix].trees mcc.tree
    The input trees file ([prefix].trees) is expected in the beast folder,
    and the output mcc.tree is saved there.
    """
    beast_dir = os.path.join(project_dir, "beast")
    input_trees = os.path.join(beast_dir, f"{prefix}.trees")
    output_mcc = os.path.join(beast_dir, "mcc.tree")
    treeannotator_command = ["treeannotator", "-heights", "ca", input_trees, output_mcc]
    print(f"Running TreeAnnotator on {input_trees} to produce {output_mcc}...")
    try:
        subprocess.run(treeannotator_command, check=True, cwd=beast_dir)
        print("TreeAnnotator completed successfully.")
    except Exception as e:
        print("Error running TreeAnnotator:", e)
        raise

def main():
    # Step 1: Process data, combine FASTAs/TSVs, align with MAFFT, and remove reference.
    aligned_fasta_path, prefix, project_dir = run_data_processing()
    
    # Step 2: Update BEAST XML using the new alignment.
    modified_xml_path = run_beast_xml_modification(aligned_fasta_path, prefix, project_dir)
    
    # Step 3: Run BEAST on the modified XML.
    run_beast_analysis(modified_xml_path)
    
    # Step 4: Run TreeAnnotator on the generated trees file.
    run_treeannotator(prefix, project_dir)
    
    print("Script completed successfully.")

if __name__ == "__main__":
    main()