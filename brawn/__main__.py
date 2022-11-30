""" The main entry point for the command line invocation of brawn """
import sys

from brawn._cmd import entrypoint

if __name__ == "__main__":
    sys.argv[0] = "brawn"  # 'python -m brawn' would otherwise report as "__main__: error..."
    sys.exit(entrypoint(sys.argv))
