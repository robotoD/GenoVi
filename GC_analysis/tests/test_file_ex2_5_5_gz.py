"""
Test_1
"""

import filecmp
import subprocess


def test_1():
    """Test_1"""
    subprocess.run(["python3", "./GC_analysis/GC_analysis.py",
                    "-i", "./tests/ex2.fasta",
                    "-o", "./tests/ex2_5_5_gz_test",
                    "-w", "5",
                    "-s", "5",
                    "-f", "gzip"])
    subprocess.run(["gzip", "-d", "-f", "./tests/ex2_5_5_gz_test.wig.gz"])
    assert filecmp.cmp("./tests/ex2_5_5_gz_test.wig", "./tests/ex2_5_5.wig")


if __name__ == "__main__":
    test_1()
