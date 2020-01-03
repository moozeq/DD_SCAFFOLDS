#!/usr/bin/env python3

import os
import sys
import subprocess
import argparse

modes = [
    'RINGS_WITH_LINKERS_1', 'RINGS_WITH_LINKERS_2', 'MURCKO_1', 'MURCKO_2', 'OPREA_1', 'OPREA_2', 'OPREA_3',
    'SCHUFFENHAUER_1', 'SCHUFFENHAUER_2', 'SCHUFFENHAUER_3', 'SCHUFFENHAUER_4', 'SCHUFFENHAUER_5']


def main():
    parser = argparse.ArgumentParser(description='Searching for same scaffolds in Homo sapiens enzymes inhibitors')
    parser.add_argument('target', type=str, help='enzyme target')
    parser.add_argument('mode', type=str, help='inhibitors compering method', choices=modes)
    parser.add_argument('-o', '--output', type=str, help='output file')
    args = parser.parse_args()

    download_inhibitors(args.target)
    strip(args.target)
    scaffolds = merge(args.target, args.mode)
    show_results(scaffolds)


def download_inhibitors(target: str):
    # download inhibitors from chembl
    from chembl_webresource_client.new_client import new_client
    if not os.path.isfile(target):
        print('[*] Downloading inhibitors from ChEMBL...')
        chembl_target = new_client.target
        targets = chembl_target.filter(target_synonym__icontains=target, organism='Homo sapiens').only(['target_chembl_id'])
        if len(targets) > 1:
            print('[-] Found multiple targets')
            sys.exit(1)
        elif len(targets) == 0:
            print('[-] Not found any targets')
            sys.exit(1)

        enzyme_target_id = targets[0]['target_chembl_id']
        activities = new_client.activity
        activities.query.limit = 100
        inhibitors = activities.filter(target_chembl_id=enzyme_target_id, standard_type='IC50').only(
            ['canonical_smiles', 'molecule_chembl_id'])
        print('[+] Done!')

        print('[*] Writing inhibitors to file...')
        with open(target, 'w') as f:
            for i in inhibitors:
                f.write(i['canonical_smiles'] + ' ' + i['molecule_chembl_id'] + '\n')
        print('[+] Done!')


def strip(target: str):
    # strip-it on downloaded inhibitors
    scaffolds_filename = f'{target}_scaffolds'
    if not os.path.isfile(scaffolds_filename):
        print(f'[*] Perfoming strip-it on inhibitors to file: {scaffolds_filename}')
        ret = subprocess.call(
            ['strip-it', '--inputFormat', 'smiles', '--input', target, '--output', scaffolds_filename])
        if ret:
            print('[-] strip-it error occurred')
            sys.exit(1)
        print('[+] Done!')

    return scaffolds_filename


def merge(target: str, mode: str):
    name_index = 0
    molecule_index = 1
    method_index = modes.index(mode) + 2  # because id and molecule go first

    scaffolds = {}
    print('[*] Merging molecules with same scaffold...')
    with open(target + '_scaffolds') as s:
        s.readline()
        for line in s:
            inhib = line.split()
            name = inhib[name_index]
            molecule = inhib[molecule_index]
            scaff = inhib[method_index]

            if scaff not in scaffolds:
                scaffolds[scaff] = []
            scaffolds[scaff].append(f'{name}\t{molecule}')
    print('[+] Done!\n')
    return scaffolds


def show_results(scaffolds, output: str = None):
    result = []
    for scaffold in scaffolds:
        result.append(scaffold + ':\n')
        for name in scaffolds[scaffold]:
            result.append('\t' + name + '\n')
        result.append('\n')

    if output:
        print('[*] Writing results to file...')
        with open(output, 'w') as o:
            o.writelines(result)
        print('[+] Done!\n')

    print('Scaffold:\n\t<ChEMBL id>\t<SMILE>\n')
    print(''.join(result))


if __name__ == '__main__':
    main()
