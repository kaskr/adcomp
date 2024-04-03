#!/bin/bash -e
git show $1:TMB/DESCRIPTION | grep Version | sed 's/.* //g'
