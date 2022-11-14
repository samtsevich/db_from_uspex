#!/usr/bin/env python3

from ase.io.espresso import read_espresso_out
from pathlib import Path

import argparse
import pandas as pd
import re
import ase.db


def getResFolderName(path: Path = Path.cwd()) -> Path:
    '''
    :param path: path to the folder where USPEX has been started
    :return: path of results folder
    '''

    folderNum = 1
    tmpFold = path / f'results{folderNum}'
    resFolder = tmpFold

    while tmpFold.is_dir():
        resFolder = tmpFold
        folderNum += 1
        tmpFold = path / f'results{folderNum}'
    return resFolder


def parse_ascii_table(ascii_table):
    table_re_list = re.split(r'\+[+-]+', ascii_table)
    table_list = [l.strip().replace('\n', '') for l in table_re_list if l.strip()]
    columns = [ch.strip() for ch in table_list[0].split('|') if ch.strip()]
    table = pd.DataFrame(columns=columns)
    for row in table_list[1:]:
        for lines in [line for line in row.split('||')]:
            line_part = [i.strip() for i in lines.split('|') if i]
            table = pd.concat([table, pd.DataFrame([line_part], columns=columns)], ignore_index=True)
    return table


def read_structures(calcfold: Path):
    assert calcfold.exists()
    output_file = calcfold/'output'
    structures = []
    if not output_file.exists():
        return
    with open(output_file, 'r') as fp:
        structures.extend(list(read_espresso_out(fp, index=slice(None))))
    return structures


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='parse results of USPEX calculation')
    parser.add_argument("-i", dest="input_path", default=Path.cwd(),
                        help="Path to the USPEX calculation folder. Should contain CalcFold* and results*")
    parser.add_argument("-db", dest="db_path", default='results.db',
                        help='path to database (if exist) that will store all data')
    args = parser.parse_args()

    uspex_fold = Path(args.input_path)
    assert uspex_fold.is_dir()

    db_path = Path(args.db_path)
    db = ase.db.connect(db_path)

    res_folder = getResFolderName(uspex_fold)
    ext_ch_file = res_folder/'extended_convex_hull'
    good_structs_file = res_folder/'goodStructures'
    datafile = good_structs_file if good_structs_file.exists() else ext_ch_file
    assert datafile.exists()

    with open(datafile) as fp:
        x = fp.read()
        data = parse_ascii_table(x)
        # print(data)

    print(f'TOTAL: {len(data)} points on the extended CH')

    # Reading structures
    for j, id in data['ID'].items():
        folders = list(uspex_fold.glob(f'CalcFold{id}_*'))
        for f in folders:
            structures = read_structures(f)
            for i, s in enumerate(structures):
                spec_data = {'uid': int(id), 'opt_step': i, 'var_cell': 1}
                db.write(s, key_value_pairs=spec_data)
            print(f'{j+1}\tfrom folder {f.name} has been added {len(structures)} structures')

