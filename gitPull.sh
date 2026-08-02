#!/bin/bash

# 1. Ask for directory with tab completion enabled (-e)
read -e -p "Enter project path (Press Tab to autocomplete, or type 'q' to quit): " repo_dir

# Check for quit option
if [ "$repo_dir" = "q" ] || [ "$repo_dir" = "Q" ]; then
    echo "Exiting script."
    exit 0
fi

# Move to the directory if specified
if [ ! -z "$repo_dir" ]; then
    # Expand tilde (~) manually if used
    repo_dir="${repo_dir/#\~/$HOME}"
    cd "$repo_dir" || { echo "Directory not found!"; exit 1; }
fi

echo -e "\nTarget directory: $(pwd)\n"

# 2. Menu options
echo "Select a Git update method:"
echo "1) Standard Update (Pull)"
echo "2) Safe Update (Fetch + Merge)"
echo "3) Overwrite Local Changes (Hard Reset)"
echo "q) Quit"
read -p "Enter choice [1-3 or q]: " choice

case $choice in
    1)
        echo "Running standard pull..."
        git pull origin main
        ;;
    2)
        echo "Fetching and merging..."
        git fetch origin
        git merge origin/main
        ;;
    3)
        echo "Warning: This will overwrite local changes!"
        read -p "Are you sure? (y/n): " confirm
        if [ "$confirm" = "y" ] || [ "$confirm" = "Y" ]; then
            git fetch origin
            git reset --hard origin/main
        else
            echo "Operation canceled."
        fi
        ;;
    q|Q)
        echo "Exiting script."
        exit 0
        ;;
    *)
        echo "Invalid choice. Exiting."
        exit 1
        ;;
esac

