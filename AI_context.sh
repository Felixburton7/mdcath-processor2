#!/bin/bash

# Define output file (adjusted to the current project structure)
OUTPUT_FILE="/home/s_felix/mdcath-processor/AI_context.txt"

# Start writing to output file
{
    echo "Working Directory: $(pwd)"
    echo ""
    echo "File Structure:"
    tree
    echo ""
    echo "Contents of Relevant Files Below (Ignoring Binary Files):"
    echo "---------------------------------------------------------"
    
    # Search inside the mdcath subdirectories for text files (ignoring .png files)
    find src/mdcath/config src/mdcath/core src/mdcath/processing  -type f ! -name "*.png" -print0 | while IFS= read -r -d '' file; do
        if file "$file" | grep -qE "text|ASCII|UTF-8"; then
            echo "===== FILE: $file ====="
            cat "$file"
            echo ""
        fi
    done

    echo ""
    echo "======================================="
    echo "Extracting First 10 Lines from Data Directory (Ignoring Binary Files)"
    echo "======================================="
    echo ""

    # Change to the expected data directory; update this path if needed.
    if [ -d "/home/s_felix/mdcath-processor/data6" ]; then
        cd /home/s_felix/mdcath-processor/data6 || exit 1
    else
        echo "Data directory /home/s_felix/mdcath-processor/data6 does not exist."
    fi

    echo "Data Directory: $(pwd)"
    echo ""
    echo "Folder Structure in Data Directory:"
    tree
    echo ""
    echo "Extracting First 10 Lines from Each File in Data Directory (Excluding Binary & pipeline.log):"
    echo "-------------------------------------------------------------------------------------"

    find . -type f ! -name "pipeline.log" ! -name "*.png" -print0 | while IFS= read -r -d '' file; do
        if file "$file" | grep -qE "text|ASCII|UTF-8"; then
            echo "===== FILE: $file ====="
            head -n 10 "$file"
            echo ""
        fi
    done
} > "$OUTPUT_FILE"
