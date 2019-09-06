import argparse

tools = [
    'features_to_partition', 'info', 'obs2idx', 'idx2obs', 'match_partitions',
    'simulate', 'demote', 'compare',
    'merge_simulations',
]

def main():
    parser = argparse.ArgumentParser(prog='tickles')
    subparsers = parser.add_subparsers(title = 'subcommands')

    # Set up the command parsers for each subcommand
    for mod_name in tools:
        mod = getattr(__import__('biggles.' + mod_name, globals(), locals(), [], -1), mod_name)
        mod_parser = subparsers.add_parser(mod_name, help=mod.__doc__)
        mod_parser.formatter_class = argparse.ArgumentDefaultsHelpFormatter
        mod.setup_parser(mod_parser)
        mod_parser.set_defaults(func=mod.main)

    # Parse the command line and invoke the subcommand if any
    args = parser.parse_args()
    if hasattr(args, 'func') and args.func:
        args.func(args)

if __name__ == '__main__':
    main()
