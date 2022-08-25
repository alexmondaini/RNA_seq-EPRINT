from pathlib import Path
import pandas as pd

p = Path('inputs','C9-eprint-sample-details.xlsx')

df = pd.read_excel(p)
print(df['filename_R2'])
print(df['I7_barcode'])