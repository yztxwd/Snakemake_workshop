#!python3

# snakemake varaibles
input = snakemake.input[0]
output = snakemake.output[0]

import pandas as pd
import matplotlib.pyplot as plt

def main():
    df = pd.read_csv(input)

    plt.figure()
    plt.scatter('x', 'y', data=df)
    plt.savefig(output)

if __name__ == "__main__":
    main()