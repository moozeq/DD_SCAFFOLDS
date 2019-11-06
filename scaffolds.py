#!/usr/bin/env python3

import os
import sys

from chembl_webresource_client.new_client import new_client
import subprocess
import argparse
import argcomplete

if __name__ == '__main__':
    modes = [
        'RINGS_WITH_LINKERS_1', 'RINGS_WITH_LINKERS_2', 'MURCKO_1', 'MURCKO_2', 'OPREA_1', 'OPREA_2', 'OPREA_3',
        'SCHUFFENHAUER_1', 'SCHUFFENHAUER_2', 'SCHUFFENHAUER_3', 'SCHUFFENHAUER_4', 'SCHUFFENHAUER_5']

    parser = argparse.ArgumentParser(description='Searching for same scaffolds in Homo sapiens enzymes inhibitors')
    parser.add_argument('target', type=str, help='enzyme target')
    parser.add_argument('mode', type=str, help='inhibitors compering method', choices=modes)
    parser.add_argument('-o', '--output', type=str, help='output file')
    argcomplete.autocomplete(parser)
    args = parser.parse_args()

    # download inhibitors from chembl
    if not os.path.isfile(args.target):
        print('[*] Downloading inhibitors from ChEMBL...')
        target = new_client.target
        targets = target.filter(target_synonym__icontains=args.target, organism='Homo sapiens').only(
            ['target_chembl_id'])
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
        with open(args.target, 'w') as f:
            for i in inhibitors:
                f.write(i['canonical_smiles'] + ' ' + i['molecule_chembl_id'] + '\n')
        print('[+] Done!')

    # strip-it on downloaded inhibitors
    if not os.path.isfile(args.target + '_scaffolds'):
        print('[*] Perfoming strip-it on inhibitors...')
        ret = subprocess.call(
            ['strip-it', '--inputFormat', 'smiles', '--input', args.target, '--output', args.target + '_scaffolds'])
        if ret:
            print('[-] strip-it error occurred')
            sys.exit(1)
        print('[+] Done!')

    name_index = 0
    molecule_index = 1
    method_index = modes.index(args.mode) + 2  # because id and molecule go first
    scaffolds = {}

    print('[*] Merging molecules with same scaffold...')
    with open(args.target + '_scaffolds') as s:
        s.readline()
        for line in s:
            inhib = line.split()
            name = inhib[name_index]
            molecule = inhib[molecule_index]
            scaff = inhib[method_index]

            if scaff not in scaffolds:
                scaffolds[scaff] = []
            scaffolds[scaff].append(name)
    print('[+] Done!\n')

    result = []
    for scaffold in scaffolds:
        result.append(scaffold + ':\n')
        for name in scaffolds[scaffold]:
            result.append('\t' + name + '\n')
        result.append('\n')

    if args.output:
        print('[*] Writing results to file...')
        with open(args.output, 'w') as o:
            o.writelines(result)
        print('[+] Done!\n')

    print('Scaffold:\n\t<list of inhibitors ChEMBL ids>\n')
    print(''.join(result))

