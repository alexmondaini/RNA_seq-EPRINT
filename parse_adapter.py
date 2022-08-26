from pathlib import Path

p = Path('inputs','adapters.txt')

data = []
with p.open() as f:
    for line in f:
        data.append(line.strip())

