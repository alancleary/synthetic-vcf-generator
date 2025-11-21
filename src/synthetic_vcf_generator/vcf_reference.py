import json
import mmap
from pathlib import Path

import synthetic_vcf_generator

METADATA_FILE_NAME = "sequence_metadata.json"


class ReferenceData:
    def __init__(self, reference_file):
        self.reference_file = reference_file

    def __enter__(self):
        self.open()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def open(self):
        self.file_obj = open(self.reference_file, mode="rb")
        self.mmap_obj = mmap.mmap(
            self.file_obj.fileno(), length=0, access=mmap.ACCESS_READ
        )

    def close(self):
        self.mmap_obj.close()
        self.file_obj.close()

    def ref_length(self):
        return self.mmap_obj.size()

    def get_ref_at_pos(self, position):
        return self.mmap_obj[position : position + 1].decode(encoding="utf-8")


def load_reference_data(reference_file):
    return ReferenceData(reference_file)


def parse_fasta(file_path, include_sequences):
    include_ids = set(include_sequences) if include_sequences else None
    parsed_ids = set()
    current_id = ""
    current_sequence = []
    with open(file_path, encoding="utf-8") as fasta_file:
        for line in fasta_file:
            line = line.strip()
            # New sequence header
            if line.startswith(">"):
                if current_id and (include_ids is None or current_id in include_ids):
                    parsed_ids.add(current_id)
                    yield {"id": current_id, "sequence": current_sequence}
                    if include_ids is not None and parsed_ids == include_ids:
                        return
                current_id = line[1:].split(" ")[0]
                current_sequence = []
            elif current_id and (include_ids is None or current_id in include_ids):
                current_sequence.extend(line.upper())
        # Add the last sequence in the file
        if current_id and (include_ids is None or current_id in include_ids):
            parsed_ids.add(current_id)
            yield {"id": current_id, "sequence": current_sequence}


def import_reference(file_path, output_dir, include_sequences=None):
    output_dir = Path(output_dir)

    if not output_dir.exists():
        output_dir.mkdir(parents=True)

    sequence_metadata_path = output_dir / METADATA_FILE_NAME

    sequence_metadata = {
        "reference_file": file_path.name,
        "synthetic-vcf-generator-version": synthetic_vcf_generator.version,
        "reference_files": {},
    }

    for parsed_sequence in parse_fasta(file_path, include_sequences=include_sequences):
        seq_file = output_dir / f"reference_{parsed_sequence['id']}.seq"
        sequence_metadata["reference_files"][parsed_sequence["id"]] = seq_file.name
        with seq_file.open("w", encoding="utf-8") as file:
            file.writelines(parsed_sequence["sequence"])

    with open(sequence_metadata_path, "w") as metadata_file:
        json.dump(sequence_metadata, metadata_file, ensure_ascii=False, indent=4)
