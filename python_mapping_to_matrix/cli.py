import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        description="a tool to compare mappings from multiple files"
    )
    parser.add_argument("--files", type=str, default="",
                        help="files to parse, comma separated")
    parser.add_argument("-i", "--input-directory", type=str, default="",
                        help="input directory")
    parser.add_argument("-m", "--max-reads", type=int, default=500000,
                        help="maximum reads to keep")
    parser.add_argument("-f", "--frequency-threshold", type=float, default=0.1,
                        help="minimum frequency threshold across files")
    parser.add_argument("-t", "--threshold", type=float, default=0.6,
                        help="private reads proportion threshold")
    parser.add_argument("-o", "--output-directory", type=str, default="output",
                        help="output directory")
    parser.add_argument("--cluster-analysis", action=argparse.BooleanOptionalAction,
                        default=True,
                        help="enable cluster analysis (default: enabled)")
    return parser.parse_args()
