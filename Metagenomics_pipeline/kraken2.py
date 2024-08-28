import os
import subprocess

def run_kraken2(forward, reverse, base_name, kraken_db, output_dir, threads):
    kraken_report = os.path.join(output_dir, f"{base_name}_report.txt")
    kraken_output = os.path.join(output_dir, f"{base_name}_kraken.txt")

    kraken_cmd = [
        "kraken2", "--db", kraken_db,
        "--threads", str(threads),
        "--report", kraken_report,
        "--output", kraken_output,
    ]

    if reverse:
        kraken_cmd.extend(["--paired", "--gzip-compressed", forward, reverse])
    else:
        kraken_cmd.extend(["--gzip-compressed", forward])

    print("Running Kraken2 command:", " ".join(kraken_cmd))  # Debug
    subprocess.run(kraken_cmd, check=True)

    return kraken_report
