from pathlib import Path
import json


r2_files = sorted(Path('/groups/cgsd/alexandre/cromwell-executions/Eprint/').glob('*/call-STAR_genome_map/execution/*bam'))
r2_files = [str(x) for x in r2_files]

data = {
    "Call_Peaks.samples": r2_files
}

output = Path('clipper_ALS.json')
with output.open('w') as f:
    json.dump(data,f,indent=4)

