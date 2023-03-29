#!/usr/bin/env python3.9
import argparse
from pprint import pprint
from os import listdir
from os.path import isfile, join, basename
import re


class Settings:
    allowed_extensions = [".tre", ".tree", ".nex", ".nxs", ".treefile"]
    verbose = False
    default_output_file = 'bootstraps.csv'
    relative_tolerance_rounding = 3
    args = None


def main():
    args = parse_args()
    Settings.args = args
    Settings.verbose = args.v

    if Settings.verbose:
        print("DEBUG-ARGS:", end=' ')
        pprint(args)

    if len(args.criteria) < 1:
        exit("No sorting criteria specified.")

    files = get_file_list(args)
    crit_tree = build_criteria_tree(args.criteria)

    # CSV TESTING - START HERE

    if args.output:
        output_file = args.output[0]
    else:
        output_file = Settings.default_output_file

    csv = CSVOutput(output_file)
    csv.write_headers(args.criteria)

    # exit("TESTING - STOP HERE")

    i_t = 0
    for file in files:
        i_t += 1
        # print(file)
        tree = PTree()
        tree.parse_file(read_tree_file(file[0]))
        seed_taxons = []

        if args.seedtaxon:
            seed_taxon_re = regize(args.seedtaxon)
            for taxon in tree.taxons:
                if re.match(seed_taxon_re, taxon.name):
                    seed_taxons.append(taxon.name)
        else:
            seed_taxons = [file[1]]
        for seed_taxon in seed_taxons:
            res = sort_one_tree_file(tree, seed_taxon, crit_tree, file[0])
            # print(i_t)
            # pprint(res)
            row = csv.csv_row_from_list(csv.list_from_result_dict(res))
            # print(row)
            csv.write_row(row)
    csv.close_file()

    # io = IOHandler(directory=args.directory, files=args.files, flist=args.list, output=args.output)
    # while f := io.read_file():
    #     print(f)


class CSVOutput:
    def __init__(self, out_path):
        self.out_path = out_path
        self.rows_written = 0
        self.ordered_crits = []
        try:
            self.file_object = open(out_path, 'w')
        except OSError as e:
            exit(f"Wrong output path: {e}")

    def close_file(self):
        self.file_object.close()

    def list_from_result_dict(self, result: dict):
        build_list = [result['file'], result['taxon'], result['size']]
        for key in self.ordered_crits:
            if result['crits'][key]:
                build_list += [result['crits'][key][0],
                               result['crits'][key][1],
                               result['crits'][key][2]]
            else:
                build_list += ['', '', '']
            pass
        return build_list

    @staticmethod
    def csv_row_from_list(items: list):
        row = ""
        items_n = len(items)
        for i, item in enumerate(items):
            # print(f"{i+1=}/{items_n} {item=}")
            row += f"\"{item}\"" + ("," if i + 1 < items_n else "")

        row += "\n"
        return row

    def write_row(self, s: str):
        self.file_object.write(s)

    def write_headers(self, crit_list: list):
        if self.rows_written == 0:
            items = ["filename", "seed_taxon", "total_taxons"]
            crit_items = []
            for item in crit_list:
                crit_items += [item, "tu_R", "tu_A"]
                self.ordered_crits += [item[:item.find("=")]]
            items += crit_items
            self.file_object.write(self.csv_row_from_list(items))
        else:
            exit("Headers already written")
            return False


def test_crit_checker():
    criteria = [
        "crit_1=10+(Aa*,Ab*),Ca*",
        "crit_2=0.125+Da*,Db*",
        "crit_3=3+Ga*,2+Ha,0.2+Ia*",
        "crit_4=Ha,Ma*,0.22+(Fa,Fe,Fo*),0.2+(Ia,Ib*,Ic),12+(Ko*,Lo,To)"]
    tolerance = "0.125"
    minimum = 2
    build_criteria_checker(criteria, tolerance, minimum)
    pass


def criteria_checker(taxons: list, criterium: list, tolerance: float, minimum: int, seed=None):
    taxons_n = len(taxons)

    tolerance_is_relative = tolerance < 1.0
    tolerance_is_absolute = not tolerance_is_relative

    used_relative_tolerance = 0.0
    used_absolute_tolerance = 0

    if taxons_n < minimum:
        # too small subtree
        return False, 0, 0

    # check for seed taxon
    if seed:
        if seed not in [_.name for _ in taxons]:
            # seed taxon not in this subtree, not a valid tree
            return False, 0, 0

    passing_taxon_set = set()

    for quantified in criterium:
        passing_taxon_count = 0
        for tested_taxon in taxons:
            for re_pattern in quantified[1]:
                if re.match(re_pattern, tested_taxon.name):
                    passing_taxon_count += 1
                    passing_taxon_set.add(tested_taxon)
                    break

        # check if passes set quantification
        if quantified[0] < 1:
            # relative quantification
            if (passing_taxon_count / taxons_n) >= quantified[0]:
                # enough occurences of specified set
                pass
            else:
                # not a valid bootstrap
                return False, 0, 0
        else:
            # absolute quantification
            if passing_taxon_count > quantified[0]:
                # enough occurences of specified set
                pass
            else:
                # not a valid bootstrap
                return False, 0, 0

    # check if passes total tolerance after all criteria are tested
    passing_taxon_count = len(passing_taxon_set)
    if tolerance_is_relative:
        # relative quantification
        used_relative_tolerance = ((taxons_n - passing_taxon_count) / taxons_n)
        if used_relative_tolerance <= tolerance:
            # enough occurences of specified set
            pass
        else:
            # not a valid bootstrap
            return False, 0, 0
    elif tolerance_is_absolute:
        # absolute quantification
        used_absolute_tolerance = (taxons_n - passing_taxon_count)
        if used_absolute_tolerance <= tolerance:
            # enough occurences of specified set
            pass
        else:
            # not a valid bootstrap
            return False, 0, 0

    # if nothing else fails, it is a valid bootstrap

    if tolerance_is_relative:
        used_absolute_tolerance = round(used_relative_tolerance * taxons_n)
        # print(f"This should be quite WHOLE  : {used_absolute_tolerance}")
    elif tolerance_is_absolute:
        used_relative_tolerance = used_absolute_tolerance / taxons_n
        # print(f"This should be quite DECIMAL: {used_relative_tolerance:1.3f}")

    return True, round(used_relative_tolerance, Settings.relative_tolerance_rounding), used_absolute_tolerance

def regize(s):
    s = s.replace("*", ".*")
    s = f"^{s}$"
    return s

def build_criteria_tree(criteria: str):


    def parse_crit(s):
        crit_list = []
        q = 0
        while s:
            # print(f"{j:04d} string is: {s}")
            if s[0].isdigit():
                # print("is digit")
                i = s.find("+")
                q = float(s[:i])
                # print(f"{q=}")
                s = s[i + 1:].lstrip(',')
            elif s[0] == "(":
                # print("is (")
                i = s.find(")")
                tup = (q, [regize(_) for _ in s[1:i].split(',')])
                crit_list.append(tup)
                q = 0
                # print(f"{tup=}")
                s = s[i + 1:].lstrip(',')
            else:
                # print("else")
                i = len(s) + 1
                i = s.find(',') if s.find(',') >= 0 else i
                tax = s[0:i]
                tup = (q, [regize(tax)])
                crit_list.append(tup)
                q = 0
                s = s[i:].lstrip(',')

        return crit_list

    crits_tree = {}
    for crit in criteria:
        # go through all string criteria and prepare filter
        crit = crit.split('=')
        crits_tree[crit[0]] = parse_crit(crit[1])

        # print(crits_tree)

    return crits_tree


def sort_one_tree_file(tree, seed_taxon, crit_tree, tree_path):
    highest_bootstrap = -1
    lowest_rel_tolerance_used = 1.0
    lowest_abs_tolerance_used = 10e6
    subtree_size = 0

    valid_count = 0

    results = dict()

    results['file'] = basename(tree_path)
    results['taxon'] = seed_taxon
    results['crits'] = {}

    for column, criterium in crit_tree.items():
        for edge in tree.edges:
            for d in range(2):
                taxons = edge.get_subtree_taxons(d)
                is_valid, r_t, a_t = criteria_checker(taxons, criterium, float(Settings.args.tolerance),
                                                      int(Settings.args.mintaxons), seed=seed_taxon)
                if is_valid:
                    valid_count += 1
                    # print(f"Bootstrap {edge.bs} is {'VALID' if is_valid else 'INVALID'} with tolerance {r_t}/{a_t}")

                    if is_valid and edge.bs > highest_bootstrap:
                        # new best bootstrap
                        highest_bootstrap = edge.bs
                        subtree_size = len(taxons)
                        lowest_rel_tolerance_used = r_t
                        lowest_abs_tolerance_used = a_t
                    elif edge.bs == highest_bootstrap and r_t < lowest_rel_tolerance_used:
                        lowest_rel_tolerance_used = r_t
                        lowest_abs_tolerance_used = a_t

        if valid_count > 0:
            results['size'] = subtree_size
            results['crits'][column] = (highest_bootstrap, lowest_rel_tolerance_used, lowest_abs_tolerance_used)
            pass
        else:
            results['size'] = ''
            results['crits'][column] = None
            pass

    return results


def read_tree_file(path):
    with open(path, 'r') as f:
        data = ''.join(f.readlines()).rstrip()
    return data


def strip_name_to_taxon(path):
    bn = basename(path)
    return bn[:bn.find(".")]


def get_file_list(args):
    files = []
    # filename and root taxon
    if args.directory:
        files = get_files_in_dir(args.directory[0])
        files = [[f, strip_name_to_taxon(f)] for f in files]
    elif args.list:
        files = parse_csv_input(args.list[0])
        if not files[0][1]:
            files = [[f[0], strip_name_to_taxon(f[0])] for f in files]
    elif args.files:
        files = [[f, strip_name_to_taxon(f)] for f in args.files]
    else:
        exit("No source of tree files specified.")

    # pprint(files)

    missing_files = 0
    for file in files:
        if not isfile(file[0]):
            missing_files += 1
            print(f"File '{file[0]}' does not exist")
    if missing_files > 0:
        exit(f"EXIT: There were {missing_files} non-existent files.")

    return files


def parse_csv_input(path):
    # print(f"{isfile(path)=}")
    if not isfile(path):
        exit(f"EXIT: {path} does not exist")

    files = []
    with open(path) as csv_file:
        i = 0
        while line := csv_file.readline().rstrip():
            if line[0] == '#':
                continue
            i += 1
            line_items = line.split(',')
            files.append([line_items[0].strip("\""), line_items[1].strip("\"") if len(line_items) > 1 else None])
            # if line_items[1]:

    nones = 0
    for file in files:
        nones += 1 if not file[1] else 0

    if nones == 0 or nones == len(files):
        # good
        return files
        # print(f"{nones}/{len(files)}")
    else:
        # print(f"{nones}/{len(files)}")
        exit(f"Seed taxon definitions in {path} are inconsistent ({nones} undefined, {len(files) - nones} defined)")


def get_files_in_dir(tree_dir):
    if tree_dir[-1] != '/':
        tree_dir += '/'
    files_unfiltered = listdir(tree_dir)
    files = []
    for name in files_unfiltered:
        for ext in Settings.allowed_extensions:
            if name[-len(ext):] == ext:
                files.append(name)
                continue

    files = [tree_dir + _ for _ in files]
    # print(f"total of {len(files)} files")
    return files


class PTree:
    def __init__(self):
        self.root = None
        self.edges = []
        self.taxons = []
        self.nodes = []

        # self.subtree_parse(self.tree_date)

    def parse_file(self, data):
        # start building tree from data
        # get root taxon name from the start
        data = data.rstrip()[1:-2]
        data = re.sub(r':[0-9]+\.[0-9]*', '', data)

        i = data.find(',')
        self.root = self.Taxon(data[:i])
        self.taxons.append(self.root)
        # run recursive tree builder
        self.node_recurse(data[i + 1:], self.root)

        for edge in self.edges:
            for obj in edge.nodes:
                if isinstance(obj, self.Taxon):
                    obj.edge = edge

        for node in self.nodes:
            for obj in node.edges:
                if isinstance(obj, self.Taxon):
                    obj.edge = node

    @staticmethod
    def opposite_bracket(data, i):
        subsegment = data[i:]
        i_o = 0
        depth = 0
        for c in subsegment:
            i_o += 1
            if c == '(':
                depth += 1
            if c == ')':
                depth -= 1
                if depth == 0:
                    return i + i_o - 1

    def node_recurse(self, data, conn):
        # create a new node to build around
        new_node = self.Node()
        self.nodes.append(new_node)
        # link "parent" element to this node
        new_node.edges[0] = conn

        # go through data string and find corresponding two elements
        i_left_bracket = 0
        if data[i_left_bracket] == '(':
            # starts with subtree, find its end bracked and recursively process
            i_right_bracket = self.opposite_bracket(data, i_left_bracket)
            bso = data[i_right_bracket + 1:].find(',')
            if bso > 0:
                bootstrap = int(data[i_right_bracket + 1:i_right_bracket + 1 + bso])
            else:
                # identical sequences not bearing any information, let bootstrap be 1
                bootstrap = 1
            # create connecting edge
            edge = self.Edge()
            self.edges.append(edge)
            # connect the edge to current node
            edge.nodes[0] = new_node
            # connect the edge to new subtree
            edge.nodes[1] = self.node_recurse(data[i_left_bracket + 1:i_right_bracket], edge)
            add_to_node = edge
            edge.bs = bootstrap
            # increase position in data to process next element
            i_right_bracket += bso + 2
        else:
            # starts with taxon (leaf), create and append the taxon
            i_right_bracket = data[i_left_bracket:].find(',')
            add_to_node = self.Taxon(data[i_left_bracket:i_right_bracket])
            self.taxons.append(add_to_node)
            i_right_bracket += 1

        # append new link to this node
        new_node.edges[1] = add_to_node

        # continue to the other node of this fork
        i_left_bracket = i_right_bracket
        if data[i_left_bracket] == '(':
            # element is a subtree, process it
            i_right_bracket = self.opposite_bracket(data, i_left_bracket)
            edge = self.Edge()
            self.edges.append(edge)

            edge.nodes[0] = new_node
            edge.nodes[1] = self.node_recurse(data[i_left_bracket + 1:i_right_bracket], edge)
            add_to_node = edge
            if i_right_bracket + 1 < len(data):
                bootstrap = int(data[i_right_bracket + 1:])
            else:
                # if bootstrap is missing
                bootstrap = 1
            edge.bs = bootstrap
        else:
            # element is a taxon, create and connect it
            add_to_node = self.Taxon(data[i_left_bracket:])
            self.taxons.append(add_to_node)

        new_node.edges[2] = add_to_node
        # new node was created and will be returned to be connected to parent node
        return new_node

    class Taxon:
        def __init__(self, name):
            self.name = name
            self.edge = None

    class Edge:
        def __init__(self):
            self.bs = None
            self.nodes = [None, None]

        def get_other_node(self, node):
            if self.nodes[0] == node:
                return self.nodes[1]
            elif self.nodes[1] == node:
                return self.nodes[0]
            else:
                raise Exception(f"Node {node} is not connected to Edge {self}")

        def get_subtree_taxons(self, d):
            st_taxons = []
            d %= 2

            def rec_gather(node, bad_edge):
                for obj in node.edges:
                    if isinstance(obj, PTree.Taxon):
                        st_taxons.append(obj)
                    elif isinstance(obj, PTree.Edge) and not obj == bad_edge:
                        rec_gather(obj.get_other_node(node), obj)
                # print(st_taxons)

            rec_gather(self.nodes[d], self)

            return st_taxons

    class Node:
        def __init__(self):
            self.edges = [None, None, None]


def parse_args():
    parser = argparse.ArgumentParser(
        prog='TreeSorter',
        description="Analyses your phylogenetic tree files to determine highest bootstrap\
                        for specified subtree criteria.")
    files_source = parser.add_mutually_exclusive_group()
    files_source.add_argument('-d', '--directory', nargs=1,
                              help="Directory of tree files to analyze")
    files_source.add_argument('-l', '--list', nargs=1,
                              help="Files containing list of files to analyze")
    files_source.add_argument('-f', '--files', nargs='*',
                              help=".tre files to analyze")

    root_taxon = parser.add_mutually_exclusive_group(required=True)
    root_taxon.add_argument('-s', '--seedtaxon', nargs='?', default='',
                            help="Use seed taxon defined in list file or from filename (in this order of preference)")
    root_taxon.add_argument('-n', action='store_true',
                            help="Don't use root taxon")

    parser.add_argument('-m', '--mintaxons', default=0,
                        help="Minimal number of taxons that pass criteria for a bootstrap to be acceptable.")

    parser.add_argument('-t', '--tolerance', default=0,
                        help="Absolute (1, 2) or relative (0.05, 0.1) tolerance \
                            of taxons outside of a set. Default is no tolerance.")

    parser.add_argument('-o', '--output', nargs=1,
                        help="Output CSV file")
    parser.add_argument('-v', action='store_true',
                        help="Run verbose")

    parser.add_argument('criteria', nargs='*',
                        help="String list of criteria to apply on subtrees")

    return parser.parse_args()


if __name__ == "__main__":
    main()
    # test_crit_checker()
    # test_build_from_file()