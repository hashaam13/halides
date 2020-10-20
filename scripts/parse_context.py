import re
import pandas as pd
import argparse
from glob import glob
import os

def main(args):
    
    dir_path = args.input.rstrip('/')
    g = glob(f'{dir_path}')

    df = pd.DataFrame({0:[]})
    for i in g:
        if os.stat(i).st_size != 0:
              df = pd.concat([df, pd.read_csv(i, sep="\t", header=None)], axis=0)
    df.index = [i for i in range(0, len(df))]

    model = df[0].str.split(":", expand=True)
    model_name = model[0]
    asa = model[1]
    resolution = model[2]
    halide = model[3]
    atom_id = model[4]
    n_chain = model[5]

    datka = pd.DataFrame({'model_name': model_name,
                          'asa': asa,
                          'resolution': resolution,
                          'halide': halide,
                          'id':atom_id,
							'n_chain':n_chain})
    datka = pd.concat([datka, df[1].str.split(',')], axis=1).rename(index=str,columns={1:'samples'})
    spread_datka = datka.apply(lambda x: pd.Series(x['samples']),axis=1).stack().reset_index(level=1, drop=True)
    spread_datka.name = 'sample'
    datka = datka.drop('samples', axis=1).join(spread_datka)
    datka.index = [i for i in range(0, len(datka))]
    datka = pd.concat([datka, datka['sample'].str.split(':', expand=True).\
             rename(columns={0:'atom_name', 1:'residue_aa', 2:'distance', 3:'angle',4:'atom_type',5: 'chain_index', 6: 'secondary_structure'})], axis=1).\
             drop('sample', axis=1)

    datka.to_csv(args.output, sep='\t', index=None)
    
if __name__=='__main__':

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)  
    parser.add_argument('-f', nargs='+', type=str)

    parser.add_argument('-input', type=str,
                        help='Path to context file or dir with context files (glob pattern should be usefull).')
    parser.add_argument('-output', type=str, help='Output file path.')

    args = parser.parse_args()
    main(args)
