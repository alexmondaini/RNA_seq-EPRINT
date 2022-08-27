from pathlib import Path
import json

files = sorted(Path('/groups/cgsd/alexandre/EPRINT_workflow/inputs/10774_ep').glob('*R2*'))
p = Path('inputs','adapters.txt')

adapters = []
with p.open() as f:
    for adapter in f:
        adapters.append(adapter.strip())


def create_data():
    data = []
    for file in files:
        data.append(
            {
                "Eprint.hg19_tar": "/groups/cgsd/alexandre/EPRINT_workflow/inputs/HG19.tar",
                "Eprint.hg19_dup_tar": "/groups/cgsd/alexandre/EPRINT_workflow/inputs/REP.tar",
                "Eprint.samples": {
                "fastq_r2": f"{file}",
                "adapters": adapters,
                "bc_pattern": "NNNNNNNNNN"
            }
            }
        )
    
    return data


if __name__=='__main__':
    output = Path('eprint_ALS.json')
    with output.open('w') as f:
        json.dump(create_data(),f,indent=4)


