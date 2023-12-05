import matplotlib.pyplot as plt
import pandas as pd

fig = plt.figure(figsize=(16,8),
                 dpi = 200)
in_file = snakemake.input[0]
out_file = snakemake.output[0]

df_file = pd.read_csv(in_file, sep='\t', header=None, names=['col1','col2','col3'])  

plt.xlabel('Index')  
plt.ylabel('Depth')
plt.title('Mapping chromosome depth')


plt.plot("col2",
         "col3",
         data=df_file)

plt.savefig(out_file)
