#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# name: format_and_calculation.py

import argparse
import os
import sys

import pandas as pd

from submodule import Module
from submodule import bracket_parser
from submodule import graph_convert


def parser_settings():
    parser = argparse.ArgumentParser(
        prog='format_and_calculation.py',
        description=(
            'Format KEGG MODULE DEFINITION to graph style data.'
            'And, calculate Reaction Coverage (RC) in KEGG MODULE.'
        ))
    # general argument
    parser.add_argument(
        'ko_file', action='store', type=str,
        help='Input File Name(KEGG ORTHOLOGY list).')
    parser.add_argument(
        '-d', '--definition_file', action='store', type=str,
        help='KEGG MODULE DEFINITION File Name.')
    parser.add_argument(
        '-r', '--result_file', action='store', type=str,
        help='Result File Name.')
    # optional argument
    parser.add_argument(
        '-m', '--method', action='store', type=str,
        choices=['minimum', 'mean', 'sum', 'bool'], default='mean',
        help='KO Abundance Mapping Method.')
    parser.add_argument(
        '-t', '--threshold', action='store', type=float,
        default=0.0, help='Threshold for KO exist.')
    args = parser.parse_args()
    return (parser, args)


def main(ko_file, result_file, definition_file, method, threshold):
    # FORMAT KEGG MODULE DEFINITION
    if definition_file is None:
        base = os.path.dirname(os.path.abspath(__file__))
        name = os.path.normpath(
            os.path.join(base, '../data/module_definition.tsv'))
        if os.path.exists(name):
            definition_file = name
        else:
            message = (
                'Module definition file does not exist.\n'
                'Please run src/download.py')
            raise ValueError(message)

    module_definitions = dict()
    with open(definition_file, 'r') as file:
        for line in file:
            line_list = line.strip('\n').split('\t')
            module_definitions[line_list[0]] = line_list[1]

    module_graphs = dict()
    for module in module_definitions:
        parsed, _ = bracket_parser(module_definitions[module], dict(), 0)
        module_graphs[module], _ = graph_convert(parsed, '0', Module())

    # MAPPING KEGG ORTHOLOGY TO FORMATTED KEGG MODULE
    # CALCULATION REACTION COVERAGE IN KEGG MODULE AND OUTPUT FILE
    matrix = pd.read_csv(ko_file, sep='\t', header=0, index_col=0, comment='#')
    columns = matrix.columns
    out = pd.DataFrame(index=module_graphs.keys(), columns=columns)
    out = out.fillna(0)
    for col in columns:
        ko_scores = matrix.ix[:, col].to_dict()
        for module in module_graphs:
            module_graphs[module].map(ko_scores, method)
            out.ix[module, col] = \
                module_graphs[module].reactionCoverage(threshold)
    out = out.sort_index()

    if result_file is not None:
        out.to_csv(result_file, sep='\t')
    else:
        out.to_csv(sys.stdout, sep='\t')


if __name__ == '__main__':
    parser, args = parser_settings()
    ko_file = args.ko_file
    result_file = args.result_file
    definition_file = args.definition_file
    method = args.method
    threshold = args.threshold

    main(ko_file, result_file, definition_file, method, threshold)
