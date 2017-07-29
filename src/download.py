#!/usr/bin/env python
# -*- coding: utf-8 -*-
# name: download.py
# author: ynakamura

import argparse
from urllib.request import urlopen


def parser_settings():
    parser = argparse.ArgumentParser(
        prog='download.py',
        description='Download KEGG MODULE DEFINITION.\n')
    # general argument
    parser.add_argument(
        '-d', '--definition_file', action='store', type=str,
        default='module_definition.tsv',
        help='Downloaded KEGG MODULE DEFINITION File Name.')
    args = parser.parse_args()
    return (parser, args)


def main(definition_file):
    with urlopen('http://rest.kegg.jp/list/module') as responce:
        module_list = \
            [r.decode('utf-8')[3:9] for r in responce.readlines()]

    data = str()
    for module in module_list:
        url = \
            'http://togows.org/entry/kegg-module/{}/definition'.format(module)
        with urlopen(url) as responce:
            data += '{0}\t{1}\n'.format(
                module, responce.read().decode('utf-8').strip('\n'))

    with open(definition_file, 'w') as outfile:
        outfile.write(data)


if __name__ == '__main__':
    parser, args = parser_settings()
    definition_file = args.definition_file

    main(definition_file)
