"""
.. See the NOTICE file distributed with this work for additional information
   regarding copyright ownership.
   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at
       http://www.apache.org/licenses/LICENSE-2.0
   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
"""

"""
Test_1
"""

import filecmp
import subprocess


def test_1():
    """Test_1"""
    subprocess.run(["python3", "./GC_analysis/GC_analysis.py",
                    "-i", "./tests/ex1.fasta",
                    "-o", "./tests/ex1_5_5_ot_bw_test",
                    "-w", "5",
                    "-s", "5",
                    "-ot",
                    "-f", "bigwig"])
    subprocess.run(["./bigWigToWig", "./tests/ex1_5_5_ot_bw_test.bw", "./tests/ex1_5_5_ot_bw_test.wig"])
    assert filecmp.cmp("./tests/ex1_5_5_ot_bw_test.wig", "./tests/ex1_5_5_ot_bw.wig")


if __name__ == "__main__":
    test_1()
