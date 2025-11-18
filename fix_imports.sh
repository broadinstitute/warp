#!/bin/bash

# Fix import paths in all WDL files to use the new flattened structure

find pipelines/wdl -name "*.wdl" -exec sed -i 's|tasks/broad/|tasks/wdl/|g' {} \;
find pipelines/wdl -name "*.wdl" -exec sed -i 's|tasks/skylab/|tasks/wdl/|g' {} \;

echo "Import paths updated in all WDL files"