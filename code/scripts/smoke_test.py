#!/usr/bin/env python3
"""Smoke test: build a tiny DB from the release and verify one query hits."""

import os
import subprocess
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FASTA = os.path.join(ROOT, "release", "sequences", "maize_pathogens_all.fasta")
MAKEBLASTDB = os.environ.get("MAKEBLASTDB", "makeblastdb")
BLASTN = os.environ.get("BLASTN", "blastn")


def main():
    records = []
    header = None
    seq = []
    with open(FASTA, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                if header:
                    records.append((header, "".join(seq)))
                header = line[1:].strip()
                seq = []
            else:
                seq.append(line.strip())
        if header:
            records.append((header, "".join(seq)))

    with tempfile.TemporaryDirectory() as tmp:
        tiny = os.path.join(tmp, "tiny.fasta")
        with open(tiny, "w", encoding="utf-8") as fh:
            for header, seq in records[:10]:
                fh.write(">" + header + "\n")
                for i in range(0, len(seq), 80):
                    fh.write(seq[i:i + 80] + "\n")
        db = os.path.join(tmp, "tiny_db")
        subprocess.run([MAKEBLASTDB,
                        "-in", tiny, "-dbtype", "nucl", "-out", db],
                       check=True, capture_output=True)
        query = os.path.join(tmp, "query.fasta")
        with open(query, "w", encoding="utf-8") as fh:
            fh.write(">query\n" + records[0][1] + "\n")
        result = subprocess.run(
            [BLASTN, "-query", query, "-db", db,
             "-outfmt", "6 sseqid pident", "-max_target_seqs", "1"],
            check=True, capture_output=True, text=True)
        line = result.stdout.strip().splitlines()[0]
        assert line.split("\t")[1] == "100.000", line
        print("SMOKE TEST PASS:", line)


if __name__ == "__main__":
    main()
