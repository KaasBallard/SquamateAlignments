#!/bin/bash

# This script is used to start or attach to a tmux session for the SquamateAlignments project

# Check if the session exists, and attach if it does
tmux has-session -t SquamateAlignments 2>/dev/null

# $? returns 0 if the session exists, so create it if it doesn't
if [ $? != 0 ]; then
  # Create a new session
  tmux new-session -d -s SquamateAlignments -n main  # The dash -d flag is used to prevent the tmux session from immediately attaching to the terminal

  # Create a new window named 'htop'
  tmux new-window -d -n htop
fi

# Attach to the session
tmux attach-session -t SquamateAlignments:main