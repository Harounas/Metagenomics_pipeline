import os
import subprocess

def run_trimmomatic(forward, reverse, base_name, output_dir, threads):
    trimmed_forward = os.path.join(output_dir, f"{base_name}_trimmed_R1.fastq.gz")
    trimmed_reverse = os.path.join(output_dir, f"{base_name}_trimmed_R2.fastq.gz")
    unpaired_forward = os.path.join(output_dir, f"{base_name}_unpaired_R1.fastq.gz")
    unpaired_reverse = os.path.join(output_dir, f"{base_name}_unpaired_R2.fastq.gz")

    if reverse:
        trimmomatic_cmd = [
            "trimmomatic", "PE", "-threads", str(threads),
            forward, reverse,
            trimmed_forward, unpaired_forward,
            trimmed_reverse, unpaired_reverse,
            "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:15", "MINLEN:36"
        ]
    else:
        trimmed_forward = os.path.join(output_dir, f"{base_name}_trimmed_R1.fastq.gz")
        trimmomatic_cmd = [
            "trimmomatic", "SE", "-threads", str(threads),
            forward, trimmed_forward,
            "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:15", "MINLEN:36"
        ]

    print("Running Trimmomatic command:", " ".join(trimmomatic_cmd))  # Debug
    subprocess.run(trimmomatic_cmd, check=True)
    
    return trimmed_forward, trimmed_reverse if reverse else None
