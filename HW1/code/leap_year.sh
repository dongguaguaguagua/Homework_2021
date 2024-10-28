#!/bin/bash

# Function to check if the year is a leap year
is_leap_year() {
	local year=$1
	if ((year % 400 == 0)); then
		echo "$year is a leap year."
	elif ((year % 100 == 0)); then
		echo "$year is not a leap year."
	elif ((year % 4 == 0)); then
		echo "$year is a leap year."
	else
		echo "$year is not a leap year."
	fi
}

# Main script logic
if [ -z "$1" ]; then
	# Interactive mode
	read -p "Enter a year: " year
	is_leap_year "$year"
else
	# Positional mode
	for year in "$@"; do
		is_leap_year "$year"
	done
fi
