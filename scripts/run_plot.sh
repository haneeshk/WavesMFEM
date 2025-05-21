#!/bin/bash

# read testName from JSON file
TEST_NAME=$(jq -r '.testName' scripts/SimulationParameters.json)
echo "testName = $TEST_NAME"

# run the Wolfram script
wolframscript -file scripts/plot_point_displacement.wls "$TEST_NAME"

# open the generated PNG
PNG_PATH="results/${TEST_NAME}/pointDisplacement.png"
if [ -f "$PNG_PATH" ]; then
    open "$PNG_PATH"
else
    echo "Could not find $PNG_PATH"
fi