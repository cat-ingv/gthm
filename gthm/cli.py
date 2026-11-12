import argparse
import sys

import gthm.workflow as gthm

def get_parser_args():
    parser = argparse.ArgumentParser(description='Script and input.')
    parser.add_argument('--subd', type=str, help='please, put subduction name', required=True)
    parser.add_argument('--poi', type=str, help='please, put the poi id', required=True)
    parser.add_argument('--subd_id', type=str, help='please, put the subduction id', required=True)
    parser.add_argument('--rupt_file', type=str, help='please, put the path of your ruptures input file', required=True)
    parser.add_argument('--output_path', type=str, help='provide dir for output_folder', required=True)
    parser.add_argument('--plot', action='store_true', help='activate the creation of automatic figures')
    parser.add_argument('--gf_path', type=str, help='please, put the path of your greens function local folder', required=True)

    return parser.parse_args()


def main():
    args = get_parser_args()
    try:
        
        gthm.run(args)
        return 0
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1

