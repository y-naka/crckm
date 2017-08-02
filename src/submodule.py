#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# name: submodule.py

import copy

import numpy as np


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

    def map(self, data, method='mean'):
        for n in range(len(self._graph)):
            sub_units = self._graph[n]['name'].split('+')
            scores = list()
            for s in sub_units:
                essencial = s.split('-')[0]
                if essencial == '':
                    scores.append(1)
                else:
                    if essencial in data:
                        scores.append(data[essencial])
                    else:
                        scores.append(0)
            if method == 'minimum':
                self.addScore(n, np.amin(scores))
            elif method == 'mean':
                self.addScore(n, np.mean(scores))
            elif method == 'sum':
                self.addScore(n, np.sum(scores))
            elif method == 'bool':
                if np.amin(scores) > 0:
                    self.addScore(n, 1)
            else:
                message = 'Undefined mapping method: {}'.format(method)
                raise ValueError(message)

    def reactionCoverage(self, threshold=0.0, keys=None,
                         score_sum=0, num=0, calculated_score=0):
        if keys is None:
            keys = self._last_steps
        if len(keys) == 0:
            result = float(score_sum) / num
            if calculated_score < result:
                calculated_score = result
        else:
            for k in keys:
                add = 0
                if self._graph[k]['score'] > threshold:
                    add = 1
                calculated_score = self.reactionCoverage(
                    threshold, self._graph[k]['source'], score_sum + add,
                    num + 1, calculated_score)
        return calculated_score


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
                elif d in '+-':
                    current = [c + d for c in current]
            else:
                module, inheritance = graph_convert(
                    parsed, ortholog, module, re_pre, re_pre)
                if d == ' ':
                    re_pre = copy.deepcopy(inheritance)
                elif d == ',':
                    re_pre = ab_pre
                    for i in inheritance:
                        products.append(i)
                elif d in '+-':
                    current = [c + d for c in current]
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
