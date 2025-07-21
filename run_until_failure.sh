#!/bin/bash

# Number of iterations
MAX_RUNS=100

# Path to your executable
EXEC="./build/debug/WaveDynamics2D"

for ((i=1; i<=MAX_RUNS; i++))
do
    echo "Run #$i..."
    
    # Run the program and capture both stdout and stderr
    OUTPUT=$($EXEC 2>&1)
    
    # Check for failure message
    if echo "$OUTPUT" | grep -q "Verification failed:"; then
        echo "❌ Detected failure in run #$i"
        echo "$OUTPUT"
        exit 1
    else
        echo "✅ Run #$i completed successfully"
    fi
done

echo "🎉 All $MAX_RUNS runs completed without verification failure."