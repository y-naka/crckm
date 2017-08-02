#!/usr/bin/env python
# -*- coding: utf-8 -*-
# name: download.py
# author: ynakamura

import argparse
import os
import re

from urllib.request import urlopen


def parser_settings():
    parser = argparse.ArgumentParser(
        prog='download.py',
        description='Download KEGG MODULE DEFINITION.\n')
    # general argument
    parser.add_argument(
        '-d', '--definition_file', action='store', type=str,
        help='Downloaded KEGG MODULE DEFINITION File Name.')
    args = parser.parse_args()
    return (parser, args)


def main(definition_file):
    with urlopen('http://rest.kegg.jp/list/module') as responce:
        module_list = \
            [r.decode('utf-8')[3:9] for r in responce.readlines()]

    data = dict()
    m_in_m = list()
    for module in module_list:
        url = \
            'http://togows.org/entry/kegg-module/{}/definition'.format(module)
        with urlopen(url) as responce:
            # data += '{0}\t{1}\n'.format(
            #     module, responce.read().decode('utf-8').strip('\n'))
            data[module] = responce.read().decode('utf-8').strip('\n')
            if 'M' in data[module]:
                m_in_m.append(module)
    for module in m_in_m:
        iterator = re.finditer('M\d\d\d\d\d', data[module])
        for match in iterator:
            data[module] = re.sub(
                match.group(), data[match.group()], data[module], 1)

    if definition_file is None:
        base = os.path.dirname(os.path.abspath(__file__))
        definition_file = os.path.normpath(
            os.path.join(base, '../data/module_definition.tsv'))

    with open(definition_file, 'w') as outfile:
        outfile.write(
            ''.join('{0}\t{1}\n'.format(k, v) for k, v in data.items()))


if __name__ == '__main__':
    parser, args = parser_settings()
    definition_file = args.definition_file

    main(definition_file)
