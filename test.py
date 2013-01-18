import doctest
import importlib

module_names = (
    "synbio",
    "synbio.clump",
    "synbio.sequence",
    "synbio.shape",
    "synbio.shapes",
)

def test(name):
    return doctest.testmod(importlib.import_module(name))

def main():
    results = [(test(name), name) for name in module_names]
    for result, name in results:
        print("{0} tests in {1}".format(result.attempted, name))

main()
