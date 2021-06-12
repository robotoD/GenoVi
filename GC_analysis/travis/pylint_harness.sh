#!/bin/bash

# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

disabled="--disable=similarities,invalid-name,too-many-statements,too-many-arguments,too-many-locals,too-few-public-methods,relative-import,no-self-use"

# pylint ${disabled} --rcfile pylintrc process*.py > output.err
# pylint ${disabled} --rcfile pylintrc tool >> output.err
pylint ${disabled} --rcfile pylintrc GC_analysis > output.err

grep -v "\-\-\-\-\-\-\-\-\-" output.err | grep -v "Your code has been rated" | grep -v "\n\n" | sed '/^$/d' > pylint.err


if [ -s pylint.err ]
then
    cat pylint.err
    rm pylint.err
    exit 1
fi
