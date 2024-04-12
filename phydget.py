#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fileencoding=utf-8

"""
PhyDGET: Phylogenetic Differential Gene Expression Tool
Author: James B. Pease
"""

from time import time
import os
import sys
import argparse
import subprocess
from multiprocessing import Manager, Pool
import numpy as np
from scipy.stats import hmean


_LICENSE = """
http://www.github.org/jbpease/phydget
PhyDGET: Phylogenetic Differential Gene Expression Tool

This file is part of PhyDGET

PhyDGET is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

PhyDGET is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with PhyDGET.  If not, see <http://www.gnu.org/licenses/>.
"""


def transform_data(datafilepath, mode='log2cpm'):
    with open(datafilepath) as infile:
        headers = infile.readline().rstrip().split('\t')
        libsizes = dict.fromkeys(headers[1:], 0)
        for line in infile:
            line = line.rstrip().split('\t')
            for i in range(1, len(line)):
                libsizes[headers[i]] += float(line[i])
    with open(datafilepath) as infile, open(
            datafilepath + ".transform", 'w') as transform_file:
        transform_file.write(infile.readline())
        for line in infile:
            line = line.rstrip().split('\t')
            gene = line[0]
            if mode == 'log2':
                transform_data = [np.log2(float(x) + 0.5) for x in line[1:]]
            elif mode == 'log2cpm':
                transform_data = [np.log2((float(x) + 0.5) /
                                  libsizes[headers[i + 1]] * 1e6)
                                  for i, x in enumerate(line[1:])]
            elif mode == 'cpm':
                transform_data = [
                    (float(x) + 0.5) / libsizes[headers[i + 1]] * 1e6
                    for i, x in enumerate(line[1:])]
            transform_file.write(
                "\t".join([gene] + [str(x) for x in transform_data]) +
                "\n")
    return libsizes


def read_config(configfilepath):
    """Reads the config file and returns command line arguments
    """
    argdict = {}
    args = []
    argorder = []
    with open(configfilepath, 'r') as cfile:
        for line in [l.strip() for l in cfile.readlines()]:
            if not line or line[0] == '#':
                continue
            if '=' in line:
                line = line.replace("=")
            elems = [d.strip() for d in line.split()]
            if elems[0] not in argdict:
                argdict[elems[0]] = None if len(elems) == 1 else elems[1:]
            else:
                argdict[elems[0]].extend(elems[1:])
            if elems[0] not in argorder:
                argorder.append(elems[0])
    for arg in argorder:
        if not argdict[arg]:
            args.append(arg)
        else:
            args.extend([arg] + argdict[arg])
    return args


def parse_models(models, species):
    """Parse the models. Syntax should be NAME:R1A+R1B,R1C:R2D
    """
    # models = ['NAME:R1A+R1B,R1C:R2D']
    # species = ["R1A", "Z"]
    ret = {'models': {},
           'model_order': [],
           'all_species': [],
           'sample_groups': {}}
    for entry in species:
        entry = entry.split(':')
        ret['all_species'].append(entry[0])
        ret['sample_groups'][entry[0]] = entry[1].split(",")
    ret['all_species'] = list(set(ret['all_species']))
    ret['models']['null'] = {'groups': [],
                             'fg': ([],),
                             'bg': (ret['all_species'],)}
    ret['model_order'].append("null")
    for modelentry in models:
        modelentry = modelentry.split(':')
        modelname = modelentry[0]
        ret['model_order'].append(modelname)
        ret['models'][modelname] = {}
        ret['models'][modelname]['groups'] = [tuple(
            tuple(y.split("+")) if "+" in y else (y,) for y in x.split(","))
            for x in modelentry[1:]]
        ret['models'][modelname]['fg'] = []
        ret['models'][modelname]['bg'] = []
        for groups in ret['models'][modelname]['groups']:
            fg = []
            for subgroup in groups:
                fg.extend(list(subgroup))
            ret['models'][modelname]['fg'].append(tuple(set(fg)))
            ret['models'][modelname]['bg'].append(tuple([
                    x for x in ret['all_species'] if x not in fg]))
    return ret


def generate_btfiles(modelname, datafile, groups,
                     btprefix, btparams):
    """Generates the control files for BayesTrait3
    """
    if modelname == 'null':
        with open("{}".format(btprefix), 'w') as nullfile:
            nullfile.write((
                "7\n"
                "2\n"
                "Burnin {}\n"
                "Iterations {}\n"
                "PriorAll uniform {} {}\n"
                "DistData {}\n"
                "Stones {} {}\n"
                "Run").format(
                    btparams['bt_burnin'],
                    btparams['bt_iter'],
                    btparams['bt_priors'][0],
                    btparams['bt_priors'][1],
                    datafile,
                    btparams['bt_stones'],
                    btparams['bt_stoneiter']
                ))
    else:
        with open("{}".format(btprefix), 'w') as altfile:
            altfile.write((
                "7\n"
                "2\n"
                "Burnin {}\n"
                "Iterations {}\n"
                "PriorAll uniform {} {}\n"
                "DistData {}\n").format(
                    btparams['bt_burnin'],
                    btparams['bt_iter'],
                    btparams['bt_priors'][0],
                    btparams['bt_priors'][1],
                    datafile
                ))
            for i, grp in enumerate(groups):
                subtags = []
                for j, subgroup in enumerate(grp):
                    altfile.write(
                        "AddTag Target{}_{} {}\n".format(
                            i+1, j+1, " ".join(list(subgroup))))
                    subtags.append("Target{}_{}".format(i+1, j+1))
                altfile.write(
                    "LocalTransform TransBranch{} {} Branch\n".format(
                        i+1, " ".join(subtags)))
            altfile.write(("Stones {} {}\nRun").format(
                btparams['bt_stones'],
                btparams['bt_stoneiter']
            ))
    return ''


def preprocess_entries(params, queue, lock):
    """Preprocessing for all entries from the output of BayesTrait
    """
    all_samples = []
    rawdata = {}
    with open(params['datafile']) as dfile:
        dfile.readline()
        for line in dfile:
            entry = line.rstrip().split()
            rawdata[entry[0]] = dict([(x, entry[i]) for i, x in
                                      enumerate(params['headers'])
                                      if i > 0])
    k = 0
    with open(params['transform_file']) as transform_file:
        transform_file.readline()
        for entry in transform_file:
            sample = {'data': {}}
            entry = entry.rstrip().split()
            sample['gene'] = entry[0]
            if params['args.test_gene'] is not None:
                if sample['gene'] not in params['args.test_gene']:
                    continue
            sample['rawdata'] = rawdata[sample['gene']]
            sample['transform_data'] = dict([(x, entry[i]) for i, x in
                                             enumerate(params['headers'])
                                             if i > 0])
            sample['transbranchsigns'] = {'null': []}
            for model in params['models']:
                # print(model, params['models'][model])
                if model == 'null':
                    continue
                group_signs = []
                for igrp, group in enumerate(params['models'][model]['fg']):
                    # print(igrp, group)
                    fg_data = []
                    bg_data = []
                    for specname in params['models'][model]['fg'][igrp]:
                        # print(specname)
                        # print(params['col_data'])
                        sample['data'][specname] = (
                            [float(entry[x])
                             for x in params['col_data'][specname]])
                        fg_data.extend(sample['data'][specname])
                    for specname in params['models'][model]['bg'][igrp]:
                        sample['data'][specname] = (
                            [float(entry[x])
                             for x in params['col_data'][specname]])
                        bg_data.extend(sample['data'][specname])
                    group_signs.append(
                        "+" if np.mean(fg_data) > np.mean(bg_data) else "-")
                sample['transbranchsigns'][model] = group_signs[:]
            sample['btparams'] = {
                'bt_iter': params['args.bt_iter'],
                'bt_burnin': params['args.bt_burnin'],
                'bt_stones': params['args.bt_stones'],
                'bt_stoneiter': params['args.bt_stoneiter'],
                'bt_priors': params['args.bt_priors']
                }
            sample['args.bt_exec'] = params['args.bt_exec']
            sample['args.tree'] = params['args.tree']
            sample['args.temp_dir'] = params['args.temp_dir']
            sample['args.keep_files'] = params['args.keep_files']
            sample['prefix'] = params['args.temp_prefix']
            sample['lock'] = lock
            sample['queue'] = queue
            sample['models'] = params['models']
            sample['model_order'] = params['model_order']
            sample['headers'] = params['headers']
            sample['species'] = params['all_species']
            # sample['transbranchsigns'] = transbranchsigns
            sample['k'] = k + 0
            all_samples.append(sample)
            k += 1

    return all_samples


def test_gene(data):
    """Main analytical function per gene
    """
    data['outputfilepath'] = (
        os.path.join(data['args.temp_dir'],
                     '{}.{}.out'.format(data['prefix'], data['k'])))
    outfile = open(data['outputfilepath'], 'w')
    result = {}
    for model in data['models']:
        result[model] = {}
    for model in data['models']:
        testid = "{}.{}".format(data['k'], model)
        data['datafilepath.{}'.format(model)] = os.path.join(
            data['args.temp_dir'],
            '{}.{}.data'.format(data['prefix'], testid))
        data['linkfilepath.{}'.format(model)] = os.path.join(
            data['args.temp_dir'],
            '{}.{}.link'.format(data['prefix'], testid))
        data['stonesfilepath.{}'.format(model)] = os.path.join(
            data['args.temp_dir'],
            '{}.{}.data.Stones.txt'.format(data['prefix'], testid))
        data['logfilepath.{}'.format(model)] = os.path.join(
            data['args.temp_dir'],
            '{}.{}.data.Log.txt'.format(data['prefix'], testid))
        data['schedulefilepath.{}'.format(model)] = os.path.join(
            data['args.temp_dir'],
            '{}.{}.data.Schedule.txt'.format(data['prefix'], testid))
        data['btfilepath.{}'.format(model)] = os.path.join(
            data['args.temp_dir'],
            '{}.{}.bt'.format(data['prefix'], testid))
        data['treesfilepath.{}'.format(model)] = os.path.join(
            data['args.temp_dir'],
            '{}.{}.data.Output.trees'.format(data['prefix'], testid))
        with open(data['linkfilepath.{}'.format(model)], 'w') as linkfile:
            for xspec in data['data']:
                linkfile.write("{}\t{}\t{}\n".format(
                    xspec,
                    "Unlinked",
                    ",".join([str(x) for x in data['data'][xspec]])))
        with open(data['datafilepath.{}'.format(model)], 'w') as tmpfile:
            for xspec in data['data']:
                tmpfile.write("{}\t{}\n".format(
                    xspec,
                    np.log2(np.mean([2**x for x in data['data'][xspec]]))))
        print()
        generate_btfiles(model,
                         data['linkfilepath.{}'.format(model)],
                         data['models'][model]['groups'],
                         data['btfilepath.{}'.format(model)],
                         data['btparams'])
        cmd = [data['args.bt_exec'], data['args.tree'],
               data['datafilepath.{}'.format(model)]]
        print('Calling: ', ' '.join(cmd))
        proc = subprocess.Popen(
            cmd,
            stdin=open(data['btfilepath.{}'.format(model)]),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        proc.communicate()
        likelihood_entry = open(
            data['stonesfilepath.{}'.format(model)]).readlines()[-1]
        result[model]['L'] = likelihood_entry.split()[-1].strip()
        print(data['k'], model, result[model]['L'])
        # BEGIN LOG FILE PROCESSING
        result[model]['alpha'] = []
        alphas = []
        result[model]['sigma'] = []
        sigmas = []
        if model != 'null':
            result[model]['branchmeans'] = []
            transbranches = [[] for _ in range(len(data['models'][model]))]
        with open(data['logfilepath.{}'.format(model)]) as lfile:
            collection = False
            for line in lfile:
                if "Sigma^2" in line:
                    collection = True
                    continue
                if collection is False:
                    continue
                entry = line.rstrip().split("\t")
                alphas.append(float(entry[3]))
                sigmas.append(float(entry[4]))
                if model != 'null':
                    for i in range(len(data['models'][model]['groups'])):
                        transbranches[i].append(float(entry[5+i]))
        result[model]['alpha'] = (
            np.log2(hmean([2**x for x in alphas])))
        result[model]['sigma'] = (
            np.log2(hmean([2**x for x in sigmas])))
        if model != 'null':
            result[model]['branchmeans'] = [
                np.log2(hmean([2**x for x in xdata]))
                for xdata in transbranches]
    # Remove files
    if data['args.keep_files'] is False:
        for model in data['models']:
            for fileprefix in ('stonesfilepath', 'logfilepath', 'btfilepath',
                               'datafilepath', 'linkfilepath',
                               'schedulefilepath', 'treesfilepath'):
                if os.path.exists(data['{}.{}'.format(fileprefix, model)]):
                    os.remove(data['{}.{}'.format(fileprefix, model)])
    # Process DeltaL values
    for model in data['models']:
        if model == 'null':
            continue
        result[model]['BF'] = (
            2 * (float(result[model]['L']) - float(result['null']['L'])))
    sorted_likelihoods = list(sorted((float(result[x]['L']), x)
                              for x in data['models']))
    print(sorted_likelihoods)
    # print(transbranchsigns)
    best_model = sorted_likelihoods[-1][1]
    best_delta = (float(sorted_likelihoods[-1][0]) -
                  float(sorted_likelihoods[-2][0]))
    output_entry = (
        [data['gene']] +
        [data['rawdata'][x] for x in data['headers'][1:]] +
        [data['transform_data'][x] for x in data['headers'][1:]] +
        [result[model]['L'] for model in data['model_order']] +
        [result[model]['BF'] for model in data['model_order'][1:]] +
        [result[model]['alpha'] for model in data['model_order']] +
        [result[model]['sigma'] for model in data['model_order']])

    for model in data['model_order'][1:]:
        # print(data['transbranchsigns'][model])
        # print(data['models'][model])
        for i, _ in enumerate(data['models'][model]['groups']):
            output_entry.append(result[model]['branchmeans'][i])
            output_entry.append(data['transbranchsigns'][model][i])
    output_entry.append(best_model)
    output_entry.append(str(best_delta))
    outfile.write("\t".join([str(x) for x in output_entry]) + "\n")
    outfile.close()
    data['queue'].put(data['outputfilepath'])
    data['lock'].acquire()
    print(data['k'], "complete")
    sys.stdout.flush()
    data['lock'].release()
    return ''


def generate_argparser():
    """argparse generation
    """
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog=_LICENSE)
    parser.add_argument("--data",
                        type=os.path.abspath,
                        required=True,
                        help="input expression data filepath (csv format)")
    parser.add_argument("--tree",
                        type=os.path.abspath,
                        required=True,
                        help="input tree file path (Nexus format)")
    parser.add_argument('--out',
                        type=os.path.abspath,
                        required=True,
                        help="output file path (csv format)")
    parser.add_argument("--transform",
                        default="log2cpm",
                        choices=("none",
                                 "log2",
                                 "cpm",
                                 "log2cpm"),
                        help=("Data transformation type "
                              "(see manual for details)."))
    parser.add_argument("--model",
                        required=True,
                        nargs="+",
                        help=("Enter one model or specify a file path"
                              "of a text file containing multiple models."
                              "Models syntax described in the manual."))
    parser.add_argument("--sample",
                        required=True,
                        nargs="+",
                        help=("SA:A1,A1,A3 SB:B1,B2,B3, ... "
                              "Species id (must match the phylogeny tip "
                              "labels) followed by comma-separated "
                              "individual sample ids "
                              "(must match the headers on the data file)."))
    parser.add_argument("--threads",
                        type=int,
                        default=2,
                        help="Number of threads for parallelization")
    parser.add_argument("--keep-files", "--keepfiles",
                        action='store_true',
                        help="Keep all temporary files")
    parser.add_argument("--temp-dir", "--tempdir",
                        default='PhyDGETtmp',
                        help="temporary folder for files")
    parser.add_argument("--temp-prefix", "--temp-prefix",
                        default="PhyDGETRun",
                        help="Temporary Directory Prefix")
    parser.add_argument("--bt-exec", "--btexec",
                        default="BayesTraitV3",
                        help="BayesTrait executable path")
    parser.add_argument("--bt-burnin", "--btburnin",
                        type=int,
                        default=1000000,
                        help=("BayesTrait number of burn-in steps."))
    parser.add_argument("--bt-iter", "--btiter",
                        type=int,
                        default=10000000,
                        help=("BayesTrait number of iterations used"
                              "per stone in the stepping stone sampling."))
    parser.add_argument("--bt-priors", "--btpriors",
                        type=float,
                        nargs=2,
                        default=(-10, 30),
                        help=("BayesTrait uniform prior range."))
    parser.add_argument("--bt-stones", "--bt-stones",
                        type=int,
                        default=200,
                        help=("BayesTrait number of stones used"
                              "in the stepping stone sampling."))
    parser.add_argument("--bt-stoneiter", "--btstoneiter",
                        type=int,
                        default=20000,
                        help=("BayesTrait number of iterations used"
                              "per stone in the stepping stone sampling."))
    parser.add_argument("--test-gene", "--testgene",
                        help=("Enter exacty gene name from first column of"
                              "input csv file to do a test run on a "
                              "single gene."))
    parser.add_argument('--verbose',
                        action="store_true",
                        help="extra screen output")
    return parser


def main(arguments=None):
    """Main PhyDGET Function
    """
    if arguments is None:
        if (len(sys.argv) == 2 and
                sys.argv[1] not in ('-h', '--help', '--version')):
            arguments = read_config(sys.argv[1])
            print("Config file used.")
            print("Executing with arguments: ", " ".join(arguments))
        else:
            arguments = sys.argv[1:]
    time0 = time()
    parser = generate_argparser()
    args = parser.parse_args(args=arguments)
    logfile = open(args.out + ".log", 'w')
    logfile.write(" ".join(arguments) + "\n")
    params = {}
    params.update(parse_models(args.model, args.sample))
    params['col_data'] = {}
    params['args.bt_exec'] = args.bt_exec
    params['args.bt_iter'] = args.bt_iter
    params['args.bt_burnin'] = args.bt_burnin
    params['args.bt_stones'] = args.bt_stones
    params['args.bt_stoneiter'] = args.bt_stoneiter
    params['args.bt_priors'] = args.bt_priors
    params['args.tree'] = args.tree
    params['args.temp_prefix'] = args.temp_prefix
    params['args.keep_files'] = args.keep_files
    params['args.temp_dir'] = os.path.abspath(
        os.path.join(os.curdir, args.temp_dir))
    params['datafile'] = args.data
    if args.transform 1= 'none':
        params['transform_file'] = args.data + ".transform"
        libsizes = transform_data(params['datafile'], mode=args.transform)
        for name, size in libsizes.items():
            logfile.write("{}\t{}\n".format(name, size))
    else:
        params['transform_file'] = args.data
    params['args.test_gene'] = (tuple(args.test_gene.split(','))
                                if args.test_gene is not None else None)
    if not os.path.exists(params['args.temp_dir']):
        os.mkdir(params['args.temp_dir'])
    # Process Data File Headers and Prepare Voom Design File
    with open(args.data) as datafile:
        params['headers'] = datafile.readline().rstrip().split("\t")
    # Establish Output File Headers
    output_headers = ['gene']
    output_headers.extend(
        ["{}.raw".format(x) for x in params['headers'][1:]] +
        ["{}.nrm".format(x) for x in params['headers'][1:]])
    output_headers.extend(
        ["L.{}".format(x) for x in params['model_order']] +
        ["BF.{}".format(x) for x in params['model_order'][1:]] +
        ["alpha.{}".format(x) for x in params['model_order']] +
        ["sigmasq.{}".format(x) for x in params['model_order']])
    for model in params['model_order']:
        for i, _ in enumerate(params['models'][model]['groups']):
            output_headers.extend(['transB.{}.{}'.format(model, i),
                                   'signB.{}.{}'.format(model, i)])
    output_headers.append("bestmodel")
    output_headers.append("bestmargin")
    with open(args.out, 'w') as outfile:
        outfile.write("\t".join(output_headers) + "\n")
    # Prepare Voom Design File
    for specname in set(params['sample_groups']):
        params['col_data'][specname] = [
            i for i, x in enumerate(params['headers'])
            if x in params['sample_groups'][specname]]
    # Prepare Shared Objects for Multithreading
    manager = Manager()
    lock = manager.RLock()
    results_queue = manager.Queue()
    pool = Pool(processes=args.threads, maxtasksperchild=1)
    # Pre-process the data file into a list of dictionary object to prepare for
    # multi-processing
    data_entries = preprocess_entries(params, results_queue, lock)
    # Initial Pool
    pool.map(test_gene, data_entries)
    # Post-process unified output
    with open(args.out, 'w') as outfile:
        outfile.write("\t".join(output_headers) + "\n")
        while not results_queue.empty():
            tempoutfilepath = results_queue.get()
            outfile.write(open(tempoutfilepath).readline())
            if args.keep_files is False:
                os.remove(tempoutfilepath)
    print("PhyDGET completed in {} seconds".format(time() - time0))
    return ''


if __name__ == '__main__':
    main()
