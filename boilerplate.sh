#!/bin/sh
#author: Lucie Ruzickova
#script: boilerplate.sh
#date: Sept 2022
#arguments: none
#description: First shell-bash script
#tools that can be used through the terminal
#can automare usage of commands, create tools, scripts and programs
#file backups, directories and files, converting formats
#using uppercase snakecase and no spaces
#MY_VAR=value as an argument
echo -e "\nThis is a shell script! \n"
#shell variables are $0 - filename
#$n where n=integer is the position of the argument
#$# the number of arguments in a script
#$@ all the arguments are individually printed
#asigning variables: MY_VAR=myvalue
#read MY_VAR
#command substitution MY_VAR=$(command) such as MY_VAR( (ls | wc -1) )
#exit