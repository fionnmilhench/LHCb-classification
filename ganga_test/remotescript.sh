#!/usr/bin/env bash

ARGS=("$@")
SHELLCMD="${ARGS[@]:0}"
echo Running $SHELLCMD
$SHELLCMD 
