#!/usr/bin/env python
# -*- coding: utf-8 -*-
# name: format_and_calculation.py
# author: ynakamura

import argparse
import copy


class Module():
    '''class for module graph data parsed from KEGG MODULE DEFINITION.'''
    def __init__(self):
        self._graph = list()
        self._count = -1
        self._last_steps = list()
        self._reaction_coverage = 0

    def addReaction(self, name, pre_reaction):
        self._count += 1
        node = {
            'name': name,
            'score': 0,
            'source': pre_reaction
        }
        self._graph.append(node)

    def addScore(self, node_id, score):
        self._graph[node_id]['score'] = score

    def setLastSteps(self, last_steps):
        self._last_steps = last_steps

    def getGraph(self):
        return self._graph

    def getCount(self):
        return self._count

    def getLastSteps(self):
        return self._last_steps

    def calculateReactionCoverage(
            self, keys=None, score_sum=0, num=0, calculated_score=0):
        if keys is None:
            keys = self._last_steps
        if len(keys) == 0:
            result = float(score_sum) / num
            if calculated_score < result:
                calculated_score = result
        else:
            for k in keys:
                calculated_score = self.calculateReactionCoverage(
                    self._graph[k]['source'],
                    score_sum + self._graph[k]['score'],
                    num + 1, calculated_score)
        return calculated_score


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
        'result_file', action='store', type=str,
        help='Result File Name.')
    parser.add_argument(
        'definition_file', action='store', type=str,
        help='KEGG MODULE DEFINITION File Name.')
    args = parser.parse_args()
    return (parser, args)


def bracket_parser(definition='', parsed=dict(), count=0):
    current = str(count)
    parsed[current] = ''
    inner_bracket = ''
    bracket_depth = 0
    for d in definition:
        if d == ')':
            bracket_depth -= 1
            if bracket_depth == 0:
                count += 1
                parsed[current] += str(count)
                parsed, count = bracket_parser(inner_bracket, parsed, count)
                inner_bracket = ''
        if bracket_depth >= 1:
            inner_bracket += d
        else:
            if d not in '()':
                parsed[current] += d
        if d == '(':
            bracket_depth += 1
    return parsed, count


def graph_convert(parsed, key, module, ab_pre=[], re_pre=[]):
    definition = parsed[key] + ','
    ortholog = ''
    products = []
    current = ['']
    complex_flag = False

    for d in definition:
        if d in '+-':
            complex_flag = True

        if d not in ', +-':
            ortholog += d
            continue

        if ortholog in parsed:
            if complex_flag:
                module, inheritance = graph_convert_complex(
                    parsed, ortholog, module, re_pre, re_pre, current)
                current = copy.deepcopy(inheritance)
                if d == ' ':
                    pre_steps = list()
                    for c in current:
                        module.addReaction(c, copy.deepcopy(re_pre))
                        pre_steps.append(module.getCount())
                    re_pre = copy.deepcopy(pre_steps)
                    current = ['']
                elif d == ',':
                    for c in current:
                        module.addReaction(c, copy.deepcopy(re_pre))
                        products.append(module.getCount())
                    re_pre = copy.deepcopy(ab_pre)
                    current = ['']
            else:
                module, inheritance = graph_convert(
                    parsed, ortholog, module, re_pre, re_pre)
                if d == ' ':
                    re_pre = copy.deepcopy(inheritance)
                elif d == ',':
                    re_pre = ab_pre
                    for i in inheritance:
                        products.append(i)
        else:
            current = [c + ortholog for c in current]
            if d == ' ':
                pre_steps = list()
                for c in current:
                    module.addReaction(c, copy.deepcopy(re_pre))
                    pre_steps.append(module.getCount())
                re_pre = copy.deepcopy(pre_steps)
                current = ['']
            elif d == ',':
                for c in current:
                    module.addReaction(c, copy.deepcopy(re_pre))
                    products.append(module.getCount())
                re_pre = copy.deepcopy(ab_pre)
                current = ['']
            elif d in '+-':
                current = [c + d for c in current]
        ortholog = ''

        if d in ', ':
            complex_flag = False

    module.setLastSteps(copy.deepcopy(products))
    return module, products


def graph_convert_complex(parsed, key, module,
                          ab_pre=[], re_pre=[], memory=['']):
    definition = parsed[key] + ','
    ortholog = ''
    products = []
    current = ['']
    complex_flag = False

    for d in definition:
        if d in '+-':
            complex_flag = True

        if d not in ', +-':
            ortholog += d
            continue
        elif d == ' ':
            message = (
                'This module definition is biologically unnatural.\n'
                'Please check if there are some errors in this.\n'
                '{}\n').format(parsed)
            raise ValueError(message)

        if ortholog in parsed:
            if complex_flag:
                module, inheritance = graph_convert_complex(
                    parsed, ortholog, module, re_pre, re_pre, current)
                current = copy.deepcopy(inheritance)
                if d == ',':
                    for i in inheritance:
                        products.append(i)
            else:
                message = (
                    'This module definition is biologically unnatural.\n'
                    'Please check if there are some errors in this.\n'
                    '{}\n').format(parsed)
                raise ValueError(message)
        else:
            current = [c + ortholog for c in current]
            if d == ',':
                for c in current:
                    products.append(c)
                current = ['']
            elif d in '+-':
                current = [c + d for c in current]
        ortholog = ''

        if d == ',':
            complex_flag = False

    items = list()
    for m in memory:
        for p in products:
            items.append('{0}{1}'.format(m, p))
    return module, items


def mapping(kos, module_graphs):
    for module in module_graphs:
        graph = module_graphs[module].getGraph()
        for n in range(len(graph)):
            sub_units = graph[n]['name'].split('+')
            ownership = list()
            for s in sub_units:
                essencial = s.split('-')[0]
                if essencial == '':
                    ownership.append(True)
                else:
                    ownership.append(essencial in kos)
            if all(ownership):
                module_graphs[module].addScore(n, 1)
    return module_graphs


def reaction_coverage(module, key, score_sum=0, num=0, calculated_score=0):
    if len(module[key]['previous']) == 0:
        result = float(score_sum) / num
        if calculated_score < result:
            calculated_score = result
    else:
        for pre_key in module[key]['previous']:
            calculated_score = reaction_coverage(
                module, pre_key, score_sum + module[pre_key]['score'],
                num + 1, calculated_score)
    return calculated_score


def main(ko_file, result_file, definition_file):
    # FORMAT KEGG MODULE DEFINITION
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
    with open(ko_file, 'r') as file:
        ko_list = file.read().split('\n')
    module_graphs = mapping(ko_list, module_graphs)

    # CALCULATION REACTION COVERAGE IN KEGG MODULE AND OUTPUT FILE
    out = ''
    for module in module_graphs:
        out += '{0}\t{1}\n'.format(
            module, module_graphs[module].calculateReactionCoverage())

    with open(result_file, 'w') as result:
        result.write(out)


if __name__ == '__main__':
    parser, args = parser_settings()
    ko_file = args.ko_file
    result_file = args.result_file
    definition_file = args.definition_file

    main(ko_file, result_file, definition_file)
