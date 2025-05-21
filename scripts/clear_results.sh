#!/bin/bash

# Get the test name from the JSON file
TEST_NAME=$(jq -r '.testName' scripts/SimulationParameters.json)

# Construct the path to the results folder
RESULTS_DIR="results/${TEST_NAME}"

# Make sure the folder exists
if [ -d "$RESULTS_DIR" ]; then
    echo "Clearing contents of $RESULTS_DIR..."
    rm -rf "${RESULTS_DIR:?}/"*  # the ?: guard prevents accidental deletion of root if var is empty
else
    echo "Directory $RESULTS_DIR does not exist. Nothing to clear."
fi