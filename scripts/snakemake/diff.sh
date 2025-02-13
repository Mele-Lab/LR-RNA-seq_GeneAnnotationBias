#!/bin/bash

diff $1 $2 | grep '^< ' | sed 's/^< //' > $3