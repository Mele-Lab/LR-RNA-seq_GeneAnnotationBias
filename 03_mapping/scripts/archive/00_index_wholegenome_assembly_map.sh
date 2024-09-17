#!/bin/bash

module load minimap
minimap2 --MD -x {params.minimap_preset} -t {threads} --secondary=no -L -a {input.genome} {input.reads} > {TMPDIR}/$uuid