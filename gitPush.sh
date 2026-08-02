#!/bin/bash

# 1. Ask for directory with tab completion enabled (-e)
read -e -p "Enter project path (Press Tab to autocomplete, or type 'q' to quit): " repo_dir

# Check for quit option
if [ "$repo_dir" = "q" ] || [ "$repo_dir" = "Q" ]; then
    echo "Exiting script."
    exit 0
fi

# Move to the directory if specified
if [ -n "$repo_dir" ]; then
    # Expand tilde (~) manually if used
    repo_dir="${repo_dir/#\~/$HOME}"
    cd "$repo_dir" || { echo "Directory not found!"; exit 1; }
fi

echo -e "\nTarget directory: $(pwd)"

# Verify the target is actually a git repository
if ! git rev-parse --is-inside-work-tree > /dev/null 2>&1; then
    echo "Error: Directory is not a git repository."
    exit 1
fi

# Display current status so you know what you are pushing
echo -e "\n--- Current Git Status ---"
git status -s
echo "--------------------------"

# Detect current branch automatically
current_branch=$(git branch --show-current)

# 2. Menu options
echo -e "\nSelect a Git push method (Current Branch: $current_branch):"
echo "1) Standard Push (Add all changes, prompt for commit, and push)"
echo "2) Push Only (If you have already committed locally)"
echo "3) Force Push (Overwrite remote history - Use with caution!)"
echo "q) Quit"
read -p "Enter choice [1-3 or q]: " choice

case $choice in
    1)
        echo "Staging all changes..."
        git add .
        
        # Prompt for commit message, but allow it to be blank
        read -p "Enter commit message (Press Enter for default): " commit_msg
        
        # If the user pressed Enter without typing, assign a default timestamp message
        if [ -z "$commit_msg" ]; then
            commit_msg="Automated update: $(date +'%Y-%m-%d %H:%M:%S')"
            echo "No message entered. Using default: $commit_msg"
        fi
        
        git commit -m "$commit_msg"
        echo "Pushing to origin/$current_branch..."
        git push origin "$current_branch"
        ;;

    2)
        echo "Pushing existing commits to origin/$current_branch..."
        git push origin "$current_branch"
        ;;
    3)
        echo "Warning: This will forcefully overwrite remote history on origin/$current_branch!"
        read -p "Are you sure? (y/n): " confirm
        if [ "$confirm" = "y" ] || [ "$confirm" = "Y" ]; then
            git push origin "$current_branch" --force
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
