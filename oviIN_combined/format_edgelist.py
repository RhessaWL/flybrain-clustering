"""Format an edgelist for use with the generalized modularity density code.
This will take the given edgelist and renumber the nodes 1...N,
where N is the number of unique nodes in the graph.

Output is two files, one has suffix _formatted and one has suffix _key.
[filename]_[suffix] has the renumbered edgelist
[filename]_key has the mapping from original node index to new node index

See https://github.com/prameshsingh/generalized-modularity-density"""
import argparse
import numpy as np
import pathlib
import logging

def main(args):
    add_name = args.name
    if not add_name.startswith("_"):
        add_name = "_" + add_name
    if add_name == "_key":
        logging.critical("Cannot use '_key' as filename as it is reserved for the re-indexing mapping")
        exit(1)

    sep = {"space": None,
        "comma": ",",
        "semicolon": ";"}.get(args.sep, None)

    source_path = pathlib.Path(args.file)
    output_file = source_path.with_stem(source_path.stem + add_name)
    key_file = source_path.with_stem(source_path.stem + "_key")

    # read the file once to store all the nodes
    with open(source_path, "r") as f:
        logging.info(f"Reading file {source_path}, skipping {args.skip} line(s), splitting at {sep=}")
        # logging.debug(f"Skipping first {args.skip} line(s)")
        nodes = []
        n_edges = 0
        for l in f.readlines()[args.skip:]:
            if l.strip():
                n_edges += 1
                try:
                    u,v,wt = l.split(sep)
                except ValueError as ve:
                    logging.critical(f"Failed to split line using separator {sep=}:\nError: {ve}\nOffending line: {l}")
                    exit(1)
                nodes.append(u.strip())
                nodes.append(v.strip())

    nodes = np.unique(nodes)  # nodes are now sorted
    logging.info(f"Read in {len(nodes)} nodes, {n_edges} edges")
    mapper = dict(zip(nodes, range(1, len(nodes)+1)))

    # read the file a second time to write the new file
    with open(args.file, "r") as source:
        logging.debug(f"Reading file {source_path}")
        with open(output_file, "w") as dest:
            logging.info(f"Writing formatted edgelist to {output_file}")
            for l in source.readlines()[args.skip:]:
                if l.strip():
                    u,v,wt = l.split(sep)
                    print(f"{mapper[u.strip()]} {mapper[v.strip()]} {wt.strip()}", file=dest)

    with open(key_file, "w") as f:
        logging.debug(f"Writing mapping to {key_file}")
        for v in sorted(mapper.keys()):
            print(f"{v} {mapper[v]}", file=f)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("file", help="Edge list file")
    ap.add_argument("--name", default="_formatted",
                    help="Output file will be named [file]_[suffix].[ext]. Default is 'formatted'. Underscore is optional.")
    ap.add_argument("--skip",
                    type=int, default=0,
                    help="Number of rows in the edgelist file to skip (useful if the file has a header row)")
    ap.add_argument("--sep",
                    choices=["space", "comma", "semicolon"], default="space",
                    help="How to split rows of the input file. By default, `l.split()` is used which splits at whitespace (see python docs for more info); using --s comma will use `l.split(',')`. This can be used if reading, e.g. a csv file.")
    ap.add_argument("-v", "--verbose",
                    choices=["debug", "info", "warn", "error", "critical"],
                    default="warn",
                    help="How verbose the output should be, from most verbose to least verbose. Default is 'warn'")
    args = ap.parse_args()

    logging.basicConfig(format="%(asctime)s %(levelname)s\t%(message)s",
                        datefmt='%m/%d/%Y %I:%M:%S %p',
                        level=getattr(logging, args.verbose.upper()))
    
    main(args)
