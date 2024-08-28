import os
import subprocess

def run_bowtie2(forward, reverse, base_name, bowtie2_index, output_dir, threads):
    unmapped_r1 = os.path.join(output_dir, f"{base_name}_unmapped_1.fastq.gz")
    unmapped_r2 = os.path.join(output_dir, f"{base_name}_unmapped_2.fastq.gz") if reverse else None

    bowtie2_cmd = [
        "bowtie2", "--threads", str(threads),
        "-x", bowtie2_index,
        "-1", forward, "-2", reverse,
        "--un-conc-gz", os.path.join(output_dir, f"{base_name}_unmapped_%.fastq.gz"),
        "-S", "/dev/null"
    ] if reverse else [
        "bowtie2", "--threads", str(threads),
        "-x", bowtie2_index,
        "-U", forward,
        "--un-gz", unmapped_r1,
        "-S", "/dev/null"
    ]

    print("Running Bowtie2 command:", " ".join(bowtie2_cmd))  # Debug
    subprocess.run(bowtie2_cmd, check=True)
    
    return unmapped_r1, unmapped_r2
