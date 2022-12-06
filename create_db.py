#!/usr/bin/env python3

from ase.io.espresso import read_espresso_out
from copy import copy
from pathlib import Path
from raw_parser import parse

import ase.db
import argparse
import numpy as np
import pandas as pd
import re

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

def read_params(params_file: Path):
    assert params_file.exists()
    with open(params_file, 'rt') as f:
        sections = f.read().split('#define ')

    def _process(input, definitions: dict):
        if isinstance(input, str) and input in definitions:
            input = copy(definitions[input])
        if isinstance(input, list):
            items = enumerate(input)
        elif isinstance(input, dict):
            items = input.items()
        else:
            items = []
        for i, element in items:
            if i != 'name':
                input[i] = _process(element, definitions)
        return input

    definitions = {}
    for section in sections[1:]:
        name, definition = section.split('\n', 1)
        name = name.strip()
        definitions[name] = parse(definition)
        definitions[name]['name'] = name
    return _process(parse(sections[0]), definitions)


def read_structures(calcfold: Path):
    assert calcfold.exists(), f'Please, check whether calcfold {calcfold} exists.'
    output_file = calcfold/'output'
    structures = []
    if not output_file.exists():
        return
    with open(output_file, 'r') as fp:
        structures.extend(list(read_espresso_out(fp, index=slice(None))))
    return structures

def get_metadata(params_path: Path) -> dict:
    assert params_path.exists()
    params = read_params(params_path)
    assert 'optimizer' in params, 'Seems like there is no optimizer set in input.uspex'
    assert 'stages' in params, 'Seems like there is no optimization stages set in input.uspex'

    metadata = {}
    target = params['optimizer']['target']
    if 'environmentUtility' in target:
        assert target['environmentUtility']
        system = {'system': 'surface', **target['environmentUtility']}
    else:
        system = {'system': 'bulk'}
    metadata.update(system)

    compositionSpace = target['compositionSpace']
    if 'minAt' in compositionSpace and 'maxAt' in compositionSpace:
        var_comp = int(compositionSpace['minAt'] != compositionSpace['maxAt'])
    else:
        var_comp = int(len(compositionSpace['blocks']) > 1)
    metadata.update({'var_comp': var_comp})

    opt_stages = params['stages'].copy()
    other_keys = ['commandExecutable', 'keepFolders']
    for key in other_keys:
        for stage in opt_stages:
            del stage[key]
    metadata.update({'opt_stages': opt_stages})
    return metadata

# ==========================================

# number of intermediate optimization steps that will be written to DB
num_selected_steps = 5

MODES = ['all', 'selected', 'generated']

modes_explanation = '''
Mode of reading optimization steps. Can be:\n
 * generated – only initial (generated) structures\n
 * selected – initial (generated) structure, final (optimized) structure and 4 random optimization steps.
 * all – reading all optimizations steps
'''


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='parse results of USPEX calculation')
    parser.add_argument("-i", dest="input_path", default=Path.cwd(),
                        help="Path to the USPEX calculation folder. Should contain CalcFold* and results*")
    parser.add_argument("-db", dest="db_path", default='results.db',
                        help='path to database (if exist) that will store all data')
    parser.add_argument("-m", dest="mode", default='all',
                        help=modes_explanation)
    args = parser.parse_args()

    print('-------------- START --------------')

    uspex_fold = Path(args.input_path)
    assert uspex_fold.is_dir()
    print(f'Calculation folder: \t\t{uspex_fold}')

    input_file = uspex_fold/'input.uspex'
    print(f'Input file for USPEX: \t\t{input_file}')

    res_folder = getResFolderName(uspex_fold)
    print(f'Results folder for USPEX: \t{res_folder}')
    ext_ch_file = res_folder/'extended_convex_hull'
    good_structs_file = res_folder/'goodStructures'
    datafile = good_structs_file if good_structs_file.is_file() else ext_ch_file
    assert datafile.exists()
    print(f'Data file for USPEX: \t\t{datafile}')

    assert args.mode in MODES, 'Please, check the mode you are trying to read structures'
    mode = args.mode
    print(f'Mode of writing to the DB: \t{mode}')

    db_path = Path(args.db_path)
    db = ase.db.connect(db_path, append=False)
    print(f'DB will be written to \t\t{db_path}')
    print('----------------------------------')

    metadata = get_metadata(input_file)
    metadata.update({'source': str(uspex_fold), 'mode': mode})

    num_stages = len(metadata['opt_stages'])
    is_surface = metadata['system'] == 'surface'
    is_bulk = metadata['system'] == 'bulk'
    assert is_surface or is_bulk, 'Seems like we have new type of system'
    db.metadata = metadata
    print('METADATA from USPEX has been read and appended to the database.')

    with open(datafile) as fp:
        x = fp.read()
        data = parse_ascii_table(x)

    print('----------------------------------')
    print(f'TOTAL: {len(data)} points on the extended CH\n')

    # Reading structures
    for j, id in data['ID'].items():
        folders = list(uspex_fold.glob(f'CalcFold{id}_*'))
        assert len(folders)
        for stage_id, _ in enumerate(metadata['opt_stages']):
            f = uspex_fold/f'CalcFold{id}_{stage_id+1}'

            # Last stage of surface calculation is just substrate – skipping this step
            if is_surface and stage_id+1 == num_stages:
                continue

            structures = read_structures(f)
            var_cell = np.allclose(structures[0].get_cell(), structures[-1].get_cell())
            spec_data = {'uid': int(id), 'opt_step': 'generated', 'var_cell': var_cell}
            db.write(structures[0], key_value_pairs=spec_data)
            if mode == 'generated':
                print(f'{j+1}\t {f.name} \thas been added only generated structure')
                continue
            elif mode == 'selected':
                if len(structures[1:-1]) < num_selected_steps:  # number of opt steps are less than threshold
                    selected_ids = range(1, len(structures)-1)
                else:
                    selected_ids = sorted(np.random.choice(range(1,len(structures)-1), size=num_selected_steps, replace=False))

                for i in selected_ids:
                    spec_data.update({'opt_step': i+1})
                    db.write(structures[i], key_value_pairs=spec_data)
                spec_data.update({'opt_step': 'optimized'})
                db.write(structures[-1], key_value_pairs=spec_data)
                print(f'{j+1}\t {f.name} \thave been added initial, final and some selected structures')
            elif mode == 'all':
                for i, s in enumerate(structures[1:-1]):
                    spec_data.update({'opt_step': i+1})
                    db.write(s, key_value_pairs=spec_data)
                spec_data.update({'opt_step': 'optimized'})
                db.write(structures[-1], key_value_pairs=spec_data)
                print(f'{j+1}\t {f.name} \thas been added all {len(structures)} structures')

    print('-------------- FINISH --------------')
